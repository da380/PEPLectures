# Slides 2026: what changed from last year, for review

Delete this file once reviewed. ⚠ marks a change with physics or notation content that needs your check.
Each deck follows last year's slide order and wording (Lectures 13–14 follow your 2026 Impress decks). All maths
is typeset from the 2026 notes, and notes figures replace last year's images wherever they show the same thing.
Borrowed images were extracted unchanged from last year's PDFs/odp files and cropped; sources are in
`figures/slides/SOURCES.md`. All 12 decks build with no warnings (`./build.sh`).

| Lecture | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| last year | 15* | 19* | 22 | 25 | 27 | 29 | 29 | 31 | 25 | 26 | 24 | 36 |
| now | 17 | 21 | 23 | 27 | 28 | 32 | 38 | 43 | 26 | 30 | 25 | 42 |

\* your 2026 Impress decks. Counts include last year's admin slides, now dropped, and this year's closing slide.

**Added after the first draft, across the decks:** every deck ends with the notes' "What you need to know and be
able to do" list. Lectures 16, 17 and 20 also build as `slidesNN-live.pdf`, where clicking a wavefield still
opens its movie (`figures/anim/*.webm`); upload the plain `slidesNN.pdf`.

## Lecture 13
* ⚠ Rigid-body arrow V′(t) = V(t − τ) → V′(t) = V(t). "station point" → "stationary point".
* New slide after the first Piola–Kirchhoff slide: a sketch of a surface element in the reference and deformed
  bodies, with dF = T n̂ dS = σ n̂_t dS_t (σ the Cauchy stress, as in Lecture 22's notes).

## Lecture 14
* Transverse-isotropy slide split in two: the tensor (as in the notes), then Love's A, C, F, L, N as components
  with the axis along x̂₃ and the isotropic limit (from the notes).
* ⚠ "At equilibrium, the motion is time-dependent" → time-independent.
* "third problem set" → first problem set (angular momentum), as in the notes.
* Added from the notes: R rotation, U and V symmetric positive-definite (polar decomposition); "the Lamé
  parameters λ and μ"; κ named as the bulk modulus.
* Question: the notes call the linearised stress relation **Hooke's law**; last year's slide did, your 2026 slide
  doesn't. Put it back? Last year's "first-order expansions" slide is omitted, as in your 2026 deck.

## Lecture 15
* Pearson block diagram → `L3F1`; textbook elastic-rebound panels → `L3F2` (fence photo kept).
* ⚠ n̂ throughout; ideal-fault stress glut S̄ = A_ijkl w_k n̂_l δ_Σ as in the notes; "(w₃ = 0)" on the planar fault.
* Credit "Kozdon et al. 2013" → Kozdon & Dunham (2013).
* Question: is the fence photo G. K. Gilbert's (USGS)? Assumed, not verified.

## Lecture 16
* ⚠ Homogeneity Γ(γp) = γ²Γ(p) (was λ, which clashes with Lamé).
* ⚠ Perturbation slides: Γ = Γ⁰ + s δΓ with δΓ_ik = ρ⁻¹ δA_ijkl p̂_j p̂_l, as in the notes.
* "A non-trivial solution exists iff" → "The equations of motion are satisfied iff"; "2D eigenvalue problem" →
  "2×2 symmetric eigenvalue problem".
* Slowness-surface raster → `L4F1` (TI and olivine). Plane-wave diagram redrawn (vector).
* New slide after the slowness surface: a point force in an isotropic medium and in a strongly anisotropic (TI)
  one, side by side; quasi-P and quasi-S fronts that are not circles, and cusps on the quasi-S front
  (`pepcode/scripts/slides_anisotropy.py`, with a TI subclass added to `pepseis/fd2d_psv.py`).
* Boundary-interaction slide kept, although it is in neither year's notes. The two snapshots after it are
  regenerated with the FD code in the style of `L5F1`/`L5F2` (`pepcode/scripts/slides_boundary.py`), with a
  one-line caption each.
* Questions: sources of the rock/Grand Canyon photos? Separately, the legend in the TI panel of `L4F1` overlaps
  the curves (a notes-figure issue).

## Lecture 17
* ⚠ Last year's errors: the recap showed the Christoffel equation as det[Γ − I] = 0 (a copy of the slowness-surface
  line), now [Γ(p) − I]a = 0; Maxwell's ∇×H + ∂D/∂t = J → ∇×H − ∂D/∂t = J.
* ⚠ "3D level surface of the Hamiltonian" → "3D surface lying within a level set of the ray Hamiltonian" (the
  level set is 5D in phase space); initial data "a 2D surface within this level set".
* Caustics: rays "cross", caustic = envelope of the ray family (notes). Euler's theorem step shown, marked
  non-examinable. "Example sheet 3" → first problem set.
* Animations: slides 6 and 7 show stills (`L5F1`, `L5F2`) where the GIFs belong. Two new stills (t = 0, 0.9)
  replace last year's "canonical problem" and "initial ansatz" images; script
  `pepcode/scripts/slides_L5F2_stills.py`.
* Your phase-space sketch is redrawn in TikZ (notation **x**⁰, **p**⁰, H_k(**x**, **p**) = ½).
* Question: source of the Newton prism engraving?

## Lecture 18
* Two short slides added from the notes: α, β formulae, why speeds rise with depth, discontinuities; CMB/ICB and
  the evidence for a fluid outer core and solid inner core.
* ⚠ "Travel time can be calculated **exactly** using ray theory" → "**efficiently**" (the notes' reason).
* ⚠ Error variance ς_i² (not σ_i², which clashes with the ray parameter); turning-point integrals to r = b;
  recap Hamiltonian H_k = ½‖p‖²c_k(x, p̂)².
* Small additions: T[x(σ)] = σ; "(for a down-going ray)"; "valid only where Δ and T are single-valued".
* Typos: "Rat theoretical", "Palvis", "Plus of minus", an unfinished sentence ("… and **p**"), "20s and 30s".
* Notes figures `L6F1`–`L6F5` replace the fault block, Loma Prieta record, PKP diagram, IASP91 chart and PREM plots.
* Ray-geometry drawing redrawn in TikZ, labelled hypocentre/epicentre/seismometer as on the earlier slides
  (last year: focus/seismograph).
* The three TauP plots are regenerated with `obspy.taup` and PREM in the style of `L6F3`/`L6F4`
  (`pepcode/scripts/slides_taup.py`). The first is **PKKKKP** (four K legs), as in last year's legend; the slide
  caption had been mistyped as PKKKP in my first draft.
* The two-panel 1940s/1981 comparison is now one plot of Jeffreys–Bullen (1940) against PREM
  (`pepcode/scripts/slides_jb_prem.py`; Gutenberg's model, in last year's textbook panel, is no longer shown).
* Questions: unknown sources: pioneers composite, Dziewonski and Anderson photos, the mineralogy panel. References to check: Pavlis & Das (2000), Hosseini & Sigloch (2015).

## Lecture 19
* Eight slides added for the 2026 material: pixel-basis toy problem, null space and SVD via AᵀA, singular
  values, discrepancy principle (χ² ≈ n), damping in SVD form, spike test, Bayesian toy problem, function-space
  inversion ("you need only the idea"). Last year had no L-curve slide.
* ⚠ Null-space vector m₀ → m_null (m₀ is the preferred model later); M > N → m > n; "trade-off parameter" →
  "damping parameter".
* ⚠ "Bayesian methods … only as good as the prior" → "Any choice of regularisation amounts to a choice of prior:
  the Bayesian approach makes this choice explicit."
* Credits: "De Witt et al." → de Wit et al. (2012); "Dziewonski & Woodhouse (1989)" → Woodhouse & Dziewonski.
* Questions: keep the SVD-damping and function-space slides (37 vs 29)? Sources of the convection cartoons
  and the delay-time seismogram figure?

## Lecture 20
* Eleven slides added for the 2026 material: data set (`L8F1`); steepest descent; time reversal as the device
  turning the terminal-value problem into an initial-value one; kernels with explicit (x, t), (x, T − t) and
  the minus sign on K_ρ; windowed cross-correlation delay with the sign remark; derivation of δτ̄ and the adjoint
  traction (replacing "not one you need to know"); comparison with δT = −∫δc/c² ds; `L8F3`, `L8F4`; GLAD-M35,
  uncertainty, six-model comparison, hot-spot sections (`L8F6`–`L8F9`).
* ⚠ "Density/elastic tensor – piecewise **constant**" → "piecewise **smooth**". Not stated in the notes; judgement.
* ⚠ Recap regulariser λ/2 mᵀBm (Lecture 19 notes; was (m − m₀)); "condition for a minimum – global and then
  local" → "(local) minimum"; u′, h′ → u†, h†; equations of motion with S̄ inside the divergence (equivalent).
* Wavefront-healing snapshots regenerated with the FD code (`pepcode/scripts/slides_healing.py`, drawn as `L5F1`).
  Weather-forecast aside dropped, as you asked.
* The two Liu & Tromp kernel slides now show YSpec kernels in PREM without its ocean
  (`pepcode/scripts/slides_bdkernels.py`), each with its seismogram and window: P at 60° in two bands, and S at
  60°; a new slide shows SS. The notes still use Liu & Tromp's `L8F5`.
* New slide after the 2D kernel: a Taylor test (`pepcode/scripts/slides_taylor.py`).
* ⚠ The Taylor test exposed a bug in `pepseis/kernel2d.py`: its adjoint source used a single window and the derivative
  of the unwindowed trace, not eqs (40)–(41) of the notes (s = w u, N = ∫ ṡ² dt), so the 2D kernel was about 9% off
  the measurement it is meant to linearise. Fixed to follow the notes; `L8F3`, `L8F4` and the kernel movie are
  regenerated (the shapes are unchanged to the eye), and the notes' footnote now reads "within about two per cent"
  (it said ten, which was this bug).
* The GLAD-M15 (Bozdağ et al. 2016) slide is dropped, as you asked; GLAD-M35 carries the point.
* Questions: all four GLAD-M35 slides? Sources of the contour and local-minima images; Farnetani (1997)
  reference unchecked.

## Lecture 21
* ⚠ Last year's errors: R_kj dR_kl/dt = ε_jkl Ω_k (repeated k) → ε_jml Ω_m; γ_i = R_ij γ̄_i → γ̄_j; centrifugal
  scaling 2ρ̄Ω²U → ρ̄Ω²U (notes).
* "Example sheet 4" → second problem set; "integro-differential equations"; d³x.
* The two ratio slides (gravity/elasticity, Coriolis/inertia) are one slide with a table of both ratios against
  period, from 1 s to 12 hours, taking L = v̄T up to the Earth's radius (so the gravity ratio is 4πGρ̄T²).
* Added `L9F1` (g and pressure in PREM); linearisation slide split (zeroth / first order).
* Recap: `L6F5` and `L8F6` (GLAD-M35) replace the PREM image and the GLAD-M15 globes.
* ⚠ Last slide: the Gharti et al. Cowling-vs-self-gravitating comparison → your `L9F2` yspec comparison, which
  switches off *all* gravitational effects (caption says so). The Gharti benchmark slide is dropped, as you asked.
* Questions: GLAD-M35 instead of the GLAD-M15 globes? A line on the Tisserand frame?

## Lecture 22
* ⚠ Last year's error: deviatoric stress ηV/L ∼ 10⁷ Pa → 10⁶ Pa (notes).
* Added: `L10F1` level surfaces; Clairaut's equation and shooting; "Rapidly rotating bodies (non-examinable)"
  section; Maclaurin–Jacobi sequences (`L10F3`).
* ⚠ Haumea period 3.9 h (was 4 h), axes within 1 % of a Jacobi ellipsoid, too fast to be homogeneous so it has a
  denser core; axis plot → `L10F4`. Fission described sympathetically (notes). Earth: hydrostatic 1/299.9 vs
  observed 1/298.3 (`L10F2`). Derivation now "in an appendix to the notes".
* Small additions from the notes (general materials, dp/dr + ρ dφ/dr = 0, γ = φ + ψ, γ₁, g₀, centrifugal ∝ Y₂₀).
* Questions: is the Phobos figure Myhill et al. (2025), arXiv:2508.07910? NASA credits assumed; Haumea render
  and Darwin fission illustration uncredited.

## Lecture 23
* ⚠ "Complete basis for vector fields" (README decision); added "countably many squared eigenfrequencies,
  bounded below, accumulating only at infinity".
* ⚠ Dropped "a finite number of eigenspaces with negative squared eigenfrequencies"; added "for realistic earth
  models such as PREM, calculations find no negative squared eigenfrequencies" and "self-gravitation is
  inherently destabilising".
* Time domain: "inverse Fourier transformation – not covered" → "contour integration – not examinable, but you
  will use it in the second problem set".
* New slide: finite record length (u_obs = h u, peaks broadened by ~T⁻¹).
* Bolivia 1994 → Tohoku 2011 at BFO (`L11F1`, `L11F2`).
* ⚠ Rigid motion u = a + Bx, B antisymmetric (was A); stability constant c₀ (was k); ω_k², |km⟩. The
  self-adjointness slide now shows ⟨w|H|u⟩ = ⟨u|H|w⟩* (last year repeated the definition).
* Question: the recap shows `L10F2` (density and 1/ε); density panel only?

## Lecture 24
* ⚠ Perturbation parameter ε (was s, clashing with the eigenfunction); toroidal test function's radial part w
  (was W′); bulk sound speed v_φ (was γ).
* Notes figures replace Woodhouse & Deuss and the Benioff scan (`L12F1`–`L12F4`, `L12F6`, `L12F7`); the two
  eigenfunction slides merged.
* Six slides added: Slichter/Stoneley/inner-core modes (`L12F8`); Tohoku splitting of ₀S₂, ₀S₃ (`L12F9`); kernel
  equations; full coupling (expansion in the spherical basis); direct solution method (dense matrix,
  complex-shifted frequencies); self-coupling and density (Akbarashrafi et al. 2018) with adjoint kernels
  (Adourian et al. 2024).
* Credits: Lau et al. 2017 (was 2018; the reference list says 2017); Widmer-Schnidrig & Laske (2007).
* New slide after the adjoint-kernels slide: Adourian et al. (2024) Figs 3–4, kernels of a spectrum about ₀S₂ and
  ₂S₂ for density, vs and vp (CC BY 4.0).
* Questions: keep the Feynman slide (photo source unknown)? Are three coupling/adjoint slides the right amount?
