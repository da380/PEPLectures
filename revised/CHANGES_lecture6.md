# Changes to `lecture6.tex` (Lecture 17: The spherical Earth)

Original has 19 numbered displayed equations (all `equation` environments), now labelled
`eq:1`–`eq:19` in order. The hard-coded references eq.(5), (6), (8), (12), (13), (14), (15)
all pointed at the correct equations and were converted one-to-one; no mapping corrections
were needed.

## A. Corrections and rewordings

### Hypocentre location using high-frequency waves
* "reasonable first-approximation" → "reasonable first approximation" (noun phrase, no hyphen).
* "smooth over the length-scale of" → "length scale" (house style).
* "the ordinary differential equations described in Lecture 17" → "Lecture~16" (this file *is*
  Lecture 17; ray theory was the previous lecture).
* "at the ith seismometer which has position" → "at the $i$th seismometer, which has position".
* "If we record waves at maximum frequency of 2 Hz" → "at a maximum frequency of $2\,\mathrm{Hz}$".
* Footnote marker moved after the full stop; footnote text "isn't well-defined nor is it important"
  → "isn't well-defined, nor is it important".
* "dependent on the relative locations of the seismometers to the source" → "dependent on the
  locations of the seismometers relative to the source" (clearer).
* "in term of both azimuth and distance" → "in terms of".
* Eqs. (1) and (2) given terminal commas (they run on into "where ...").

### The spherically averaged velocity structure of the Earth
* "Studies of earthquake and Earth structure are fundamentally linked" → "Studies of earthquakes
  and of Earth structure".

### Ray theory in spherically symmetric earth models
* "integrable systems which having sufficiently many conserved quantities" → "which have".
* "right hand sides" → "right-hand sides".
* "Nöether's theorem" → "Noether's theorem".
* "$q=\|\mathbf{q}\|$ is constant along a ray, and that this can be written" → "and this can be
  written".
* "up going (+) or down going (-)" → "up-going ($+$) or down-going ($-$)"; "initially down going
  ray" → "initially down-going ray"; "the down going leg" → "the down-going leg".
* "are left with single equation" → "a single equation".
* Eq. (12): `\frac{\ddns x_{i}}{d\sigma}` → `\frac{\ddns x_{i}}{\ddns\sigma}` (typo in the
  differential).
* "Combining eq.(12) and (13)" → "Combining eqs.~(\ref{eq:12}) and (\ref{eq:13})".
* "portions of the ray in which the epicentral angle $\Delta$ or the travel time $T$ are
  single-valued functions of radius $r$" → "in which the epicentral angle $\Delta$ and the travel
  time $T$ are single-valued functions of the radius $r$" (grammar: plural verb needs "and").
* **Surface radius symbol $a$ → $b$**: "a ray starting at the surface $r=a$" → "$r=b$", and the
  upper limits of the $\Delta$ and $T$ integrals `\int_{r_{t}}^{a}` → `\int_{r_{t}}^{b}` (house
  style; $b$ is the planetary radius in Lectures 21 and 23).
* Eq. (8) `q=\frac{r\sin\theta}{\alpha}` and eq. (16) `q=\frac{r_t}{\alpha(r_t)}` had no terminal
  punctuation; given a comma and a full stop respectively.
* Eq. (17) (epicentral-angle integral) ended with ";" → ",".
* "with r monotonically increasing" → "with $r$ monotonically increasing".
* "the up and down going legs of the rays path" → "the up- and down-going legs of the ray's path".
* "Using eq.(14) which is valid throughout the down going leg, we obtain" → "Using
  eq.~(\ref{eq:14}), which is valid ..., we obtain"; "Using eq.(15) we similarly" → "Using
  eq.~(\ref{eq:15}), we similarly".
* Stray spaces before the full stops in eqs. (7) and (11) (`\mathbf{0} .`, `\hat{\bm{\Delta}} .`)
  removed.

### Travel time inversions for spherical velocity structure
* "By repeating this process for many many events, global averaged travel time curves have been
  constructed" → "for many events, globally averaged travel-time curves".
* "more-or-less ad hoc methods" → "more or less ad hoc methods" (adverbial).
* Life dates (1858-1936), (1889-1960), (1891-1989), (1888-1993) → en-dashes.

### An overview of the Earth's spherically averaged velocity structure
* "The model PREM by Dziewonski \& Anderson (1981) shows the variation of P-wave and S-wave
  velocity with depth in the earth." → "Figure~\ref{fig:5} shows the variation of P-wave and
  S-wave velocity with depth in the Earth according to the Preliminary Reference Earth Model
  (PREM) of Dziewonski \& Anderson (1981)." (the sentence was describing the figure, not the
  model; acronym expanded at first use).
* "elastic modulii" → "elastic moduli".
* "their presence must be due to sharp changes in either the composition or through a phase
  change within the same material" → "must be due either to sharp changes in composition or to a
  phase change within the same material" (parallel construction).
* "core mantle boundary (CMB) which is" → "core--mantle boundary (CMB), which is"; "inner core
  boundary (ICB) which is" → "(ICB), which is".
* Fluid outer core list: "(i) observations of solid-Earth tides, (ii) the absence of S-waves ...,
  and (iii) is also required for the generation of the geodynamo" → "(iii) the fact that a fluid
  core is required for the generation of the geodynamo" (item (iii) was not parallel with (i) and
  (ii)).
* "And lastly, some studies that have identify inner core S-waves" → "Lastly, some studies have
  identified inner core S-waves".
* "the robustness of such observations are still debated" → "is still debated".
* "upper most mantle" → "uppermost mantle".
* "410 km and 670 km" (twice) → `$410\,\mathrm{km}$ and $670\,\mathrm{km}$`.
* Caption of Fig. 5: "between the the rocky mantle" → "between the rocky mantle".

### What you need to know and be able to do
* `itemize` → `enumerate` (house style; `\item[(i)]` labels unchanged).
* Item (ii): "Within the second problem set, you will work through similar calculations in a
  horizontally stratified earth model" → "Within the first problem set" (this is Problem 7 of
  Problem Set 1).

## B. Left unchanged / questions for the author
* Eq. (7) is stated for the P-wave Hamiltonian using the eikonal equation `H = 1/2`; the
  simplification `\partial H/\partial\mathbf{p} = \alpha\hat{\mathbf{p}}` relies on
  `\alpha\|\mathbf{p}\| = 1`, which is what the text says ("where we have used the eikonal
  equations"). Correct as written; noted only because the step is compressed.
* In eq. (8) $\theta$ is "measured positively clockwise", but eq. (9)/(10) treat $\theta$ as the
  angle from $\hat{\mathbf{x}}$ towards $\hat{\bm{\Delta}}$; the sense of "clockwise" depends on
  the orientation of the ray plane. Left as is.
* Eqs. (14)–(15) carry the minus sign of the down-going leg; eqs. (17)–(18) drop it by the
  symmetry argument (integral from $r_t$ to $b$, factor 2). Consistent; left as is.
* Citation style "(Dziewonski \& Gilbert 1975)", "(e.g.\ Deuss et al. 2000)" left as in the
  original, as instructed.
* "two-fold", "first arriving waves", "vertical component seismogram", "later arriving phases",
  "higher frequency components" left unhyphenated (author's usage).
* Fig. 2 caption "in Finland following an earthquake in California" — no data source given; left.
* The last paragraph of the ray-theory subsection says the extension to more complicated velocity
  structures "takes a little effort" but Fig. 3 (shadow zone) is never referred to in the text;
  no reference added, since that would be new material.

## C. Style pass
Applied: `\documentclass[a4paper,12pt]`; `||...||` → `\|...\|` (3 places); `, \quad` → `,\quad`
inside displays; `e.g.`/`i.e.` → `e.g.\ `/`i.e.\ ` before words (4 places); straight and Unicode
quotes → ``...''; Unicode apostrophe → ASCII; units to `\,\mathrm{}` form (2 Hz ×3, km s^-1,
2.5 km, 410/670 km ×2); hyphen → en-dash in life dates and "core--mantle"; figure labels
`fig:fig1..5` → `fig:1..5`; every display labelled `eq:1`–`eq:19`; double spaces in prose
collapsed. Compiles with no errors/undefined references; the two overfull boxes (lines 45–55,
56–67 of the source) are inherited from the original. Nothing in the style guide failed to apply.

## D. Second pass (25 September 2026)

* Renumbered and retitled: "Lecture 18: Spherical Earth structure" (body-only file).
* Fig. 1 redrawn in TikZ; Figs. 3–5 regenerated in PREM with `figsrc/prem.py` and `figsrc/sphray.py`; Fig. 3 cited in the text; captions of Figs. 3 and 4 rewritten to describe what is now shown.
