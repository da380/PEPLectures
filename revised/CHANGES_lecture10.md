# Changes to lecture10.tex ("Lecture 21: Equilibrium figures")

## A. Corrections and rewordings

**Summary section ("What you need to know and be able to do") -- PLEASE CHECK**
- Two items referred to material that is not in this lecture:
  - (ii) "That the equilibrium equations are under-determined, but they admit non-unique solutions subject to the two compatibility conditions." → "That the equilibrium equations are underdetermined, and why." (no compatibility conditions are discussed anywhere in the lecture; the footnote on divergence-free symmetric tensor fields is the only supporting material).
  - (vi) "How minimum stress fields can be defined and their relation to slow viscous flow. Again, you do not need to learn this derivation in full, but may be asked to perform related calculations given the necessary information." → **deleted** (minimum stress fields and slow viscous flow are not treated in this lecture).
  - Items (i), (iii), (iv), (v) keep their hand-written labels; nothing else renumbered. See section B: the author may prefer to restore the missing material rather than trim the summary.
- `itemize` → `enumerate` (STYLE.md §4).

**Outline and motivation**
- "the constraints imposed by equilibrium equations" → "imposed by the equilibrium equations" (missing article).

**The equilibrium equations**
- "planet\footnote{...}." → "planet.\footnote{...}" (footnote marker after punctuation).
- "To simplify notations" → "To simplify notation".
- "over long time-scales" → "over long time scales" (STYLE.md §4).
- "then all terms in the equilibrium equations are fixed" → "would be fixed" (sequence of tenses after "If we knew").
- "under-determined" → "underdetermined" (here and in the hydrostatic subsection; also "over-determined" → "overdetermined").
- "we cannot expect to uniquely determine stress tensor\footnote{...}." → "the stress tensor.\footnote{...}" (missing article; footnote moved after the full stop). Inside the footnote "(i.e. divergence of a curl is zero)" → "(i.e.\ the divergence of a curl is zero)".
- "Lecture 13" → "Lecture~13".

**Hydrostatic equilibrium**
- `-p~\delta_{ij}` → `-p\,\delta_{ij}` (eqs. 3 and 4).
- "by requiring the velocity be divergence-free" → "requiring the velocity to be divergence-free".
- Units: "$\eta\sim10^{21}$ Pas" → `\eta\sim10^{21}\,\mathrm{Pa\,s}`; "$V\sim10^{-8}ms^{-1}$" → `V\sim10^{-8}\,\mathrm{m\,s^{-1}}`; "$L\sim10^{6}$ m" → `L\sim10^{6}\,\mathrm{m}`; "$\eta V/L\sim10^{7}$ Pa" → `\ldots10^{7}\,\mathrm{Pa}`; "$10^{10}$Pa\footnote{...}," → "$10^{10}\,\mathrm{Pa}$,\footnote{...}".
- Eq. 5 ended with "." but the sentence continues "along with the boundary condition" → ",".
- "recalling $\nabla\times\nabla f=0$" → "$\nabla\times\nabla f=\mathbf{0}$" (the curl is a vector; consistent with the bold zeros in eqs. 9--11).
- "and so note that $\rho$ and $\gamma$ also have coincident level surfaces" → "and so $\rho$ and $\gamma$ also have ..." (stray "note that").
- "for an planet" (twice) → "for a planet" / "of a planet".

**Non-rotating hydrostatic planets**
- Eq. 12: `=0; ,` → `=0,` (stray semicolon).
- "for a non-rotating planets" → "for a non-rotating planet".
- "not easy\footnote{... "Ellipsoidal figures of equilibrium".}." → "not easy.\footnote{... ``Ellipsoidal figures of equilibrium''.}" (straight quotes → LaTeX quotes; footnote after the full stop).

**Slowly rotating hydrostatic planets**
- "and b the mean radius" → "and $b$ the mean radius"; "with M the planet's mass" → "with $M$ the planet's mass".
- Units: "$\overline{\rho}\sim5000~\text{kg m}^{-3}$" → `5000\,\mathrm{kg\,m^{-3}}`; "$\Omega\sim7\times10^{-5}\text{ s}^{-1}$" → `7\times10^{-5}\,\mathrm{s^{-1}}`.
- `\cdot\cdot\cdot` → `\cdots` (6×); `s~\phi_{1}`, `s~\psi_{1}` → `s\,\phi_{1}`, `s\,\psi_{1}`.
- "internal boundaries\footnote{But these can be easily included into the theory.}," → "boundaries,\footnote{But these can be easily included in the theory.}".
- "we find at first-order" → "at first order"; "Expanding this out to first-order in $s$" → "to first order in $s$" (adverbial use, STYLE.md §4).
- "and hence the boundary conditions simplify to" → "the boundary condition simplifies to" (there is a single condition, eq. 26).
- "We will show that eq.(22) and (26) allow" → "eqs.~(\ref{eq:22}) and (\ref{eq:26})"; "Eq.(38) and (39) constitute" → "Eqs.~(\ref{eq:38}) and (\ref{eq:39})".
- "surfaces\footnote{...}, and hence" → "surfaces,\footnote{...} and hence".
- "Clairaut (1713-1765)" → "(1713--1765)".
- "are found to oblate spheroids" → "are found to be oblate spheroids".
- `\epsilon~Y_{20} .` → `\epsilon\,Y_{20}.`; Clairaut's equation ended "=0 ." followed by "which is known as" → "=0,"; `\frac{\ddns ^{2}\epsilon}` → `\frac{\ddns^{2}\epsilon}`.
- Terminal punctuation added to displays that lacked it: eqs. 3, 4, 6, 7, 13, 14, 30, 36, 37, 38, 41, 43 now end with a comma (sentence continues); eqs. 8, 15, 31, 34 with a full stop.
- Figure: `[h]` → `[ht]`; `fig:fig1` → `fig:1`; Unicode apostrophes in "Earth’s", "Clairaut’s" → ASCII; "The lower figure then plots" → "The lower panel then plots" (parallel to "The upper panel").

**Equation references**
- The original had 43 numbered displays (39 `equation` + one 4-line `align`). All hard-coded references -- eq.(2), (4), (9)×3, (22)×3, (26), (27), (35)×2, (38)×2, (39) -- pointed at the correct equations; each is now `eq.~(\ref{eq:N})` with the same N. Every display now carries `\label{eq:N}`; the printed numbers in the PDF are unchanged from the original.

## B. Left unchanged / questions for the author

- **Summary items (ii) and (vi)**: as noted above, the "two compatibility conditions" and "minimum stress fields ... slow viscous flow" are not covered in this lecture (they read like leftovers from an earlier version in which the non-uniqueness of the equilibrium stress was resolved by a minimum-stress/viscous-flow argument). I have trimmed (ii) and deleted (vi) rather than invent content; if that material is meant to be examinable, it needs to be restored to the lecture (or to the problem set) instead.
- The figure (Earth's density profile and $1/\epsilon$ from Clairaut's equation) is never referred to in the text. A natural place for `Fig.~\ref{fig:1}` would be the sentence "... in later lectures we will see how this was determined for the Earth", or after "for the Earth this value is around $1/300$". Not added, since that would be new text.
- $g_{0}$ is defined as $\ddns\phi_{0}/\ddns r$ (positive, since $\phi_{0}$ increases outwards), so $g_{0}$ here is the magnitude of gravity while $g_{i}=-\partial\phi/\partial x_{i}$ is the acceleration vector. This is self-consistent (eqs. 22, 25, 31, 34, 40, 42, 43 all check out with this convention) but the reader might be helped by a word saying so.
- $\Omega$ is called "the angular velocity" on p.1 and "the angular frequency" in the slowly-rotating section; both are in common use, left as they are.
- "co-ordinates", "book-keeping" left as written (not covered by STYLE.md).
- The footnote after "stress tensor" in the Hydrostatic-equilibrium subsection sits mid-sentence ("... the stress tensor\footnote{...} takes the form"); there is no punctuation to move it after, so it was left in place.
- The double spaces in "restrictions on,  and interrelations between, its  shape" were left so the line stays untouched for latexdiff.

## C. Style pass

Applied: `\label{eq:1}`--`\label{eq:43}` on every display (align lines labelled individually); all `eq.(N)`/`Eq.(N) and (M)` → `eq.~(\ref{eq:N})`/`Eqs.~(\ref{eq:N}) and (\ref{eq:M})`; `Lecture 13` → `Lecture~13`; `~` → `\,` thin spaces in `-p\,\delta_{ij}`, `s\,\phi_{1}`, `s\,\psi_{1}`, `\epsilon\,Y_{20}`; `\cdot\cdot\cdot` → `\cdots`; units in `\,\mathrm{...}`; straight quotes → ``...''; hyphen → en-dash in the date range; Unicode apostrophes removed; six footnote markers moved after punctuation; `i.e.\ `; under-/over-determined and time-scales per §4; figure `[ht]`/`fig:1`; summary `itemize` → `enumerate`. Already conformant and untouched: `\ddns` derivatives, `\hat{\mathbf{r}}`, `\mathbf{0}`, `,\quad` inside eqs. 35, 37, 39, 43, title in sentence case, `\mathbf{u}`-free text. Compiles with no errors or warnings (8 pages, as the original; the original's one overfull box in summary item (ii) is gone); latexdiff runs cleanly.
