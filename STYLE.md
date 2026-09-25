# House style for the seismology lecture notes (PEP, Lectures 12–23)

This file records the formatting and style conventions applied in the revised
lecture notes, problem sets and worked examples. It is written so that a bulk
decision (e.g. "no, keep hyphens in Piola-Kirchhoff") can be reversed with a
single search-and-replace.

## 1. File skeleton

```latex
\documentclass[a4paper,12pt]{paper}      % problem sets and solutions: 10pt

\usepackage{mycommands}

\title{Lecture 12: Finite elasticity}      % sentence case after the colon
\author{David Al-Attar \\ Michaelmas Term, \the\year}

\begin{document}

\maketitle

\section*{Outline and motivation}          % lower-case m
...
\section*{What you need to know and be able to do}
\begin{enumerate}
    \item[(i)] ...
\end{enumerate}

\end{document}
```

* No `\usepackage{amsmath}` in the file: `mycommands.sty` already loads it.
* Section and subsection titles in sentence case. Non-examinable material is
  flagged as `\subsection*{Title (non-examinable)}`.
* Figures:
  ```latex
  \begin{figure}[ht]
      \centering
      \includegraphics[width=0.7\textwidth]{figures/L1F1.png}
      \caption{Sentence-case caption ending with a full stop.}
      \label{fig:1}
  \end{figure}
  ```
  Labels are `fig:1`, `fig:2`, ... within each file (not `fig:fig1`,
  `fig:slowness`). Fix captions that were pasted from a PDF with hard
  hyphenation ("inter-\naction").

## 2. Equations and cross-references

* Every displayed equation that is referred to carries `\label{eq:N}` where
  `N` is its position in the file (`eq:1`, `eq:2`, ...). In practice label
  every displayed equation, so later insertions do not force relabelling of
  the text. Multi-line `align` blocks that were meant to be one equation use
  `\nonumber` on all but the last line (or `split` inside `equation`).
* **Never hard-code equation or figure numbers in the text.** Use
  `eq.~(\ref{eq:7})`, `eqs.~(\ref{eq:7}) and (\ref{eq:8})`,
  `eqs.~(\ref{eq:7})--(\ref{eq:9})`, `Fig.~\ref{fig:1}`. At the start of a
  sentence write `Eq.~(\ref{eq:7})` / `Eqs.~...`. Lectures are cited as
  `Lecture~12`; problem sets as "the first problem set" / "the second problem
  set" (there are only two).
* Displayed equations are part of the sentence: they end with a comma or full
  stop as the grammar requires, and the following text is not capitalised
  after a comma.
* Punctuation inside displayed maths uses `,\quad` between parallel relations:
  `\mu>0,\quad \kappa>0`.

## 3. Mathematical notation

| Item | Use | Not |
|---|---|---|
| norm | `\|\mathbf{v}\|` | `||\mathbf{v}||`, `\|\!\|` |
| differential in integrals | `\dd^{3}\mathbf{x}`, `\dd S`, `\dd t`, `\dd\sigma` | `d^3x`, `dS`, `dt`, `\dd\mathbf{x}` |
| ordinary derivative | `\frac{\ddns}{\ddns t}`, `\frac{\ddns x_{i}}{\ddns\sigma}` | `\frac{d}{dt}`, `d\sigma` |
| imaginary unit | `\ii` | `i` |
| exponential | `\ee^{\ii\omega t}` | `e^{i\omega t}` |
| ellipses | `i=1,\dots,n`, `\cdots` in sums/expansions | `...`, `\cdot\cdot\cdot` |
| bold Greek | `\bm{\Gamma}`, `\bm{\nu}`, `\bm{\Omega}`, `\bphi` (macro) | `\boldsymbol{\Gamma}`, `\mathbf{\Omega}`, `\boldsymbol{\varphi}` |
| roman labels | `\mathbf{u}^{\mathrm{obs}}`, `\mathbf{m}_{\mathrm{in}}`, `\mathrm{diag}` | `^{obs}`, `_{in}`, `\text{diag}` |
| position arguments | `u_{i}(\mathbf{x},t)`, `T(\mathbf{x})` | `u_i(x,t)`, `T(x)` |
| thin spaces | `s\,\delta\alpha`, `-p\,\delta_{ij}`, `2q\alpha\,\dd z` | `s~\delta\alpha`, `-p~\delta_{ij}`, `q\alpha~\dd z` |
| tall brackets | `\left(A_{ijkl}\frac{\partial u_{k}}{\partial x_{l}}\right)` | `(A_{ijkl}\frac{\partial u_{k}}{\partial x_{l}})` |
| units | `2\,\mathrm{Hz}`, `5\,\mathrm{km\,s^{-1}}`, `10^{21}\,\mathrm{Pa\,s}`, `\mathrm{m^{3}\,kg^{-1}\,s^{-2}}` | `2Hz`, `Pas`, `ms^{-1}`, `\text{km s}^{-1}` |
| ranges | `10^{5}$--$10^{8}$` in text, `(1859--1944)` | `10^{5}-10^{8}`, `(1859-1944)` |
| variables in prose | `$s$ is a perturbation parameter`, `for each $k$`, `the $i$th seismometer` | `s is`, `for each k`, `the ith`, `the $i$-th` |
| operators in bra-kets | `\langle \mathbf{w}|P|\mathbf{u}\rangle` with plain italic `P`, `W`, `H` (as in Lecture 22) | `\mathbf{P}` |
| mean radius of the Earth | `b` | `a` |

* `\partial^{2}` for second derivatives: `\frac{\partial^{2}\zeta}{\partial x_{i}\partial x_{j}}`.
* Do not reuse a symbol that already has a meaning in the same file (e.g. do
  not use `A_{ij}` for an antisymmetric matrix in a file where `A_{ijkl}` is
  the elastic tensor).

## 4. Text conventions

* British spelling: linearised, parameterise, modelled, centre, neighbourhood.
  "moduli" (not "modulii"), "phenomenon" (singular), "different from".
* Quotation marks: ``like this'' (never straight `"`), and no Unicode curly
  quotes or apostrophes (`’`, `“`, `”`, `–`) in the source.
* Dashes: spaced en-dash ` -- ` for parenthetical dashes; en-dash for name
  pairs and ranges: Piola--Kirchhoff, Cauchy--Green, Euler--Lagrange,
  Sturm--Liouville, Runge--Lenz, Backus--Gilbert.
* Footnote markers go after punctuation: `... processes.\footnote{Text of the
  footnote, ending with a full stop.}`
* Non-breaking spaces before references: `Fig.~`, `eq.~`, `Lecture~12`,
  `Dziewonski \& Anderson (1981)`.
* Abbreviations: `e.g.\ ` and `i.e.\ ` (with the control space) when followed
  by a word; expand acronyms at first use (PREM, CMB, ICB, MCMC, LLSVP, PDF).
* Compound terms: wavefront, wavefield, waveform, travel time (noun),
  travel-time curve (adjective), length scale, time scale, half-space,
  whole space, first-order (adjective) / to first order, well-posed,
  ill-posed, underdetermined, overdetermined, non-linear, self-adjoint,
  co-rotating, squared eigenfrequency (no hyphen), zeroth-order term / at
  zeroth order, frequency domain (noun) / frequency-domain response.
* Capitalisation: "the Earth" (the planet) but "earth model" (a model, as in
  Dahlen & Tromp); P-wave, S-wave (capital letters); "Poisson's equation";
  "Hamilton's principle"; "Newton's second law".
* "data" is plural: "the data are", "these data".
* Bold (`\textbf`) is used for a term only where it is first defined or where
  a key physical statement is emphasised, as in the original notes; do not add
  new bold.
* Lists in the summary section use `enumerate` with `\item[(i)]` labels.

## 5. What is *not* changed

* The lecture numbering (12–23) and the overall structure of each lecture.
* The author's voice: first person, informal asides and footnotes stay.
* Mathematical content, unless it is an outright error; every such correction
  is recorded in `EDIT_NOTES.md` so it can be checked.
* Existing line breaks in the source wherever the text is untouched, so that
  `git diff` and `latexdiff` stay readable.
