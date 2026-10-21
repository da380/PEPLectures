# Borrowed figures for the slides

Cropped from the publisher PDF (vector graphics, rendered at 300 dpi), captions removed.

## Lecture 13 (slides13.tex, formerly slides1.odp), introduction

Source: Tromp, J., Komatitsch, D., Hjörleifsdóttir, V., Liu, Q., Zhu, H., Peter, D., Bozdağ, E., McRitchie, D.,
Friberg, P., Trabant, C. & Hutko, A., 2010. Near real-time simulations of global CMT earthquakes,
*Geophys. J. Int.*, 183, 381–389, doi:10.1111/j.1365-246X.2010.04734.x.
PDF: https://global.shakemovie.princeton.edu/pdf/j.1365-246X.2010.04734.x.pdf

Suggested slide credit: "Tromp et al. (2010), GJI, 183, 381–389".

* `tromp2010_fig7_body_waves.png` (their Fig. 7): vertical-component record section, data (black) and 3-D
  spectral-element synthetics (red; S362ANI + Crust2.0) for the 2008 September 3, Mw 6.3 Santiago del Estero,
  Argentina earthquake (depth 571 km). Aligned on P, plotted against epicentral distance (54–144°), bandpass
  17–60 s. Major body-wave arrivals labelled.
* `tromp2010_fig8_love_waves.png` (their Fig. 8): transverse-component record section, data (black) and
  3-D synthetics (red) for the 2007 August 17, Mw 6.4 Banda Sea earthquake (depth 21 km). Aligned on the Love
  wave (4.4 km/s), plotted against source azimuth, bandpass 60–200 s.

## Lecture 13 (slides13.tex), "the overall point": convection, tomography, earthquakes, finite fault

First choice listed first in each group; the others are alternatives.

Convection
* `davies2022_convection_pacific.png` (1009×1009): Davies, D. R., Kramer, S. C., Ghelichkhan, S. & Gibson, A.,
  2022. Towards automatic finite-element methods for geodynamics via Firedrake, *Geosci. Model Dev.*, 15,
  5127–5166, doi:10.5194/gmd-15-5127-2022, Fig. 11b (panel label and colour-bar ticks removed). CC BY 4.0.
  Cut-away global simulation driven by 230 Myr of plate history: hot plumes (red), slabs (blue), plate velocities.
* `davies2022_convection_africa_preprint.png` (2973×2977): same authors, discussion paper version,
  *Geosci. Model Dev. Discuss.*, doi:10.5194/gmd-2021-367, Fig. 11a (earlier model run: fewer, fatter plumes).
  CC BY 4.0. Extracted from the PDF; CMYK converted to sRGB without a colour profile, so colours are more
  saturated than published.

Tomography
* `garnero2024_semucb-wm1_globe.jpg` (2409×2700): image by Edward Garnero, adapted from model SEMUCB-WM1
  (French & Romanowicz 2014, *Geophys. J. Int.*, 199, 1303); in Garnero, E. & Richardson, C., 2024. The
  mysterious, massive structures in Earth's deep mantle, *Physics Today*, 77(12), 36. Publisher copyright (AIP),
  no licence stated. Isosurfaces of the two large low-velocity provinces in a cut-away globe.
* `pugmire_glad-m25.jpg` (1879×950): visualisation by David Pugmire, ORNL (OLCF), of GLAD-M25 (Lei et al. 2020,
  *Geophys. J. Int.*, 223, 1–21). CC BY 2.0, via Wikimedia Commons ("A Data-Driven Journey to the Center of the
  Eath"); colour-bar label cropped off.

Earthquakes
* `seismicity_M5_1976-2025.png` (3200×1333): made by `pepcode/scripts/slides_seismicity.py`. USGS ComCat
  catalogue, M ≥ 5, 1976–2025 (83 901 events), coloured by depth, over NASA Blue Marble Next Generation
  (topography and bathymetry, December 2004; NASA Earth Observatory). Both public domain. Suggested credit:
  "Data: USGS; background: NASA Blue Marble".

Finite fault
* `gsi2011_tohoku_slip.png` (1900×2250): Geospatial Information Authority of Japan (GSI) and Japan Coast
  Guard, 2011. Coseismic slip distribution on the plate interface from GPS and seafloor geodetic observations,
  2011 Tohoku earthquake, https://www.gsi.go.jp/common/000060854.pdf (Japanese title and notes cropped).
  GSI content is under the Public Data License 1.0; credit "Source: GSI website" with the URL.
* `barbot2023_turkiye_slip.png` (3895×1530): Barbot, S. et al., 2023. Slip distribution of the February 6, 2023
  Mw 7.8 and Mw 7.6, Kahramanmaraş, Turkey earthquake sequence in the East Anatolian Fault Zone, *Seismica*,
  2(3), doi:10.26443/seismica.v2i3.502, slip-distribution figure, panel (a), Mw 7.8 mainshock (background made
  white). CC BY 4.0.

## Lecture 15 (slides15.tex), seismic sources

All four images were extracted from last year's deck (`old_slides/Lecture_14_-_Seismic_sources_-_Slides.pdf`)
with `pdfimages -all` and cropped to remove the captions and margins baked into them.

* `ikluft2007_san_andreas_carrizo.jpg` (1062×704): the San Andreas Fault in the Carrizo Plain, California, aerial
  photograph by Ikluft (Ian Kluft), 16 November 2007, Wikimedia Commons,
  "File:Kluft-photo-Carrizo-Plain-Nov-2007-Img_0327.jpg", CC BY-SA 4.0. Confirmed against the Commons original
  (3072×2048); last year's slide used a contrast-enhanced version, kept here. Credit: "Ikluft, Wikimedia Commons
  (CC BY-SA 4.0)".
* `bearphoto1906_san_francisco_street.jpg` (1056×766): fissured cobbled street with streetcar tracks after the
  1906 San Francisco earthquake. The print carries the stamp "Bear Photo S.F. 260" (Bear Photo Co., San Francisco).
  Public domain (published 1906). Archive record not identified.
* `gilbert1906_fence_offset.jpg` (368×543): fence offset by the 1906 rupture (hand-tinted print), panel (b) of the
  textbook elastic-rebound figure used last year (the textbook panels (a) are replaced by `L3F2.pdf`). Usually
  credited to G. K. Gilbert, USGS (1906, near Woodville, Marin County; USGS Bulletin 324); this attribution was
  not checked against the USGS photo library record. USGS, public domain. Credit: "Photo: G. K. Gilbert, USGS (1906)".
* `kozdon2013_tohoku_wavefield.jpg` (1208×608): horizontal and vertical particle velocity at 15 s and 110 s from a
  dynamic rupture simulation of the 2011 Tohoku earthquake. Kozdon, J. E. & Dunham, E. M., 2013. Rupture to the
  trench: dynamic rupture simulations of the 11 March 2011 Tohoku earthquake, *Bull. seism. Soc. Am.*, 103(2B),
  1275–1289, doi:10.1785/0120120136. Publisher copyright. Last year's slide credited "Kozdon et al. 2013"; the
  paper has two authors, so the credit is now "Kozdon & Dunham (2013)".

## Lecture 16 (slides16.tex), plane wave propagation

The first three extracted at original resolution from last year's `slides/slides4.odp` (git `d38c2f9`), not from the PDF.

* `vandriel_nissenmeyer_global_wavefield.png` (937×928): cut-away globe showing a simulated global seismic
  wavefield. Credited on the 2025 slide to Martin van Driel & Tarje Nissen-Meyer (presumably an AxiSEM
  visualisation; original publication not traced). The credit baked into the image was painted out and the
  white margins trimmed; the slide carries `\credit{Martin van Driel \& Tarje Nissen-Meyer}`.
* `long2009_splitting_flow.png` (1248×538): shear-wave-splitting fast directions compared with plate- and
  density-driven mantle flow, (a) Atlantic and (b) Pacific basins. Credited on the 2025 slide to Long, M. D. &
  Silver, P. G., 2009. Shear wave splitting and mantle anisotropy: measurements, interpretations, and new
  directions, *Surv. Geophys.*, 30, 407–461, doi:10.1007/s10712-009-9075-1 (the panel may be reproduced there
  from earlier work; not checked). Baked-in credit painted out, margins trimmed.
* `anisotropy_rocks_photos.png` (912×735): composite of three photographs (granite hand specimen, banded gneiss
  hand specimen, Grand Canyon layered strata) used on the "Anisotropy in the Earth" slide. No credit on the 2025
  slide; sources not recorded (probably Wikimedia Commons). Margins trimmed.
* `boundary_t1.pdf`, `boundary_t2.pdf`: P-wave from an explosive point source meeting a horizontal boundary with a
  faster, denser medium, at t = 0.36 and 0.62 (incident wave; reflected and transmitted P- and S-waves and head
  waves). Made by `pepcode/scripts/slides_boundary.py` with `pepseis.fd2d_psv`, drawn as `L5F1`/`L5F2`, on a grid
  larger than the window shown, so that no wave reaches its absorbing edges before the second snapshot. They
  replace last year's two uncredited snapshots.
* `anisotropy_iso.pdf`, `anisotropy_ti.pdf`: a vertical point force in an isotropic medium and in a strongly
  anisotropic one (transversely isotropic, axis vertical; Thomsen ε = 0.3, δ = −0.1, α₀ = 1, β₀ = 0.5) at the same
  time; quasi-P and quasi-S fronts, the latter with cusps. Made by `pepcode/scripts/slides_anisotropy.py` with
  `pepseis.fd2d_psv.TIPSVSolver`, drawn as `L5F1`/`L5F2`. New in 2026.

## Lecture 17 (slides17.tex), ray theory

* `newton_prism_engraving.jpg` (1722×965): extracted with `pdfimages` from the 2025 slides
  (`old_slides/Lecture_16_-_Ray_theory_-_Slides.pdf`, p. 3), white border trimmed. Uncredited there: a
  hand-coloured 19th-century engraving of Newton's prism experiment (the one with the cat). Original source and
  licence not recorded; presumably public domain, David to confirm. No credit on the slide (as in 2025).
* `L5F2_t00.pdf`, `L5F2_t09.pdf`: stills at t = 0 and t = 0.9 of the Lecture 17 Fig. 2 finite-difference simulation
  (identical set-up to `pepcode/scripts/L5F2.py`, which gives t = 1.2), made with `pepseis.fd2d_psv.snapshot_figure`
  (vmax = 1). Script: `pepcode/scripts/slides_L5F2_stills.py`. They replace the 2025 "canonical problem" and "initial ansatz" images.

## Lecture 18 (slides18.tex), spherical Earth structure

All extracted with `pdfimages -all` from last year's `old_slides/Lecture_17_-_Spherical_Earth_-_Slides.pdf`
(JPEGs as embedded, soft masks fully opaque, so unchanged apart from cropping); page numbers refer to that PDF.
Captions that were baked into the images were cropped off and are typeset on the slides instead.

Earthquake location
* `pavlis2000_pamir_map.jpg` (583×820) and `pavlis2000_pamir_sections.jpg` (806×647), old p. 8: map of the
  Pamir–Hindu Kush region and two competing tectonic cross-sections (a), (b). Last year's slide credited
  "Palvis & Das 2000", i.e. Pavlis, G. L. & Das, S., 2000, *Tectonics* 19 (Pamir–Hindu Kush seismic zone as a
  strain marker for flow in the upper mantle); full reference to be checked. Slide credit "Pavlis & Das (2000)".

Rays in spherical earth models
* The ray-geometry drawing (hypocentre, epicentre, seismometer, Δ) is drawn in TikZ in `slides18.tex`.

Spherical velocity structure
* `rays_pkkkkp.pdf`, `rays_pcp_pcs.pdf`, `rays_pdiff_sdiff.pdf`: ray paths and travel-time curves in PREM for a
  surface source (PKKKKP; PcP and PcS; Pdiff and Sdiff), from `obspy.taup`, drawn in the style of `L6F3`/`L6F4` by
  `pepcode/scripts/slides_taup.py`. They replace last year's TauP plots; the first is PKKKKP (four K legs), as last year's legend says.
* `iris_global_stack.jpg` (951×1000), old p. 20, left half of a composite (the right half, the IASP91 travel-time
  chart, is replaced by `figures/L6F4.pdf`): "Global stacking three-component long period & broadband
  seismograms", IRIS (logo in image). Slide credit "IRIS".
* `seismology_pioneers.jpg` (1721×1067; also on the Lecture 19 recap slide), old p. 22: composite with photographs of Oldham, Lehmann, Jeffreys and
  Gutenberg, an Earth cut-away cartoon and the Oldham (1906) quotation. Sources not recorded.
* `dziewonski_photo.jpg` (250×318) and `anderson_photo.jpg` (413×616), old p. 23: photographs of Adam Dziewonski
  and Don Anderson. Sources not recorded. (The PREM plot from that slide is replaced by `figures/L6F5.pdf`.)
* `jb_prem.pdf`: P- and S-wave speeds in the Jeffreys–Bullen model (1940; obspy's TauP `jb.nd`, without its
  later inner-core S speed) and in PREM, drawn as `L6F5` by `pepcode/scripts/slides_jb_prem.py`. It replaces last
  year's two panels (a textbook figure of the Jeffreys and Gutenberg models of the 1940s, and PREM).
* `l18_mantle_mineralogy.jpg` (1038×948), old p. 27, right panel: mineral volume fractions against depth in the
  upper mantle, transition zone and top of the lower mantle, with mineral photographs. Source not recorded.
  (The annotated PREM plot on the left of that slide is replaced by `figures/L6F5.pdf`.)
* `hosseini2015_pdiff.jpg` (1253×930), old p. 29: data, synthetics and Green's functions for core-diffracted P
  waves against epicentral distance. Credited "Hosseini & Sigloch 2015" on last year's slide; presumably
  Hosseini, K. & Sigloch, K., 2015, *Geophys. J. Int.* (multifrequency Pdiff measurements); to be checked.
  Slide credit "Hosseini & Sigloch (2015)".

## Lecture 19 (slides19.tex), delay time tomography

Extracted with `pdfimages -all` from last year's `old_slides/Lecture_18_-_Delay_time_tomography_-_Slides.pdf`
(JPEG copied byte for byte unless a crop is stated). Last year's slides gave no credit for the first four.

* `seismology_pioneers.jpg`: see Lecture 18.
* `mantle_convection_cartoon.jpg` (1671×920): whole-mantle convection cartoon (ridge, trenches, convection cells).
  Source not recorded on last year's slides.
* `mantle_convection_models.jpg` (1650×889): four cartoons of competing mantle convection models, with the
  "Th & U poor / Th & U rich" labels added on last year's slide. Source not recorded.
* `delay_time_seismograms.jpg` (1639×789): observed and reference synthetic seismograms with P- and S-wave picks,
  illustrating a delay time. Source not recorded.
* `dewit2012_nullspace_model.jpg`, `dewit2012_nullspace_robust.jpg` (about 720×845): cropped from last year's
  slide images (credit text and caption removed, re-saved as JPEG at quality 92). de Wit, R. W. L.,
  Trampert, J. & van der Hilst, R. D., 2012. Toward quantifying uncertainty in travel time tomography using the
  null-space shuttle, *J. geophys. Res.*, 117, B03301, doi:10.1029/2011JB008754. Cross-sections through the
  P-wave model MIT-P08 beneath the Americas; the second shows the same sections with features not required by the
  data greyed out. Publisher copyright (AGU). Last year's slides misspelt the credit as "De Witt et al. 2012".
  Suggested slide credit: "de Wit et al. (2012)".

Also used, already in `figures/`: `L7F1.png` (Woodhouse & Dziewonski 1989), `L12F5.png` (Ritsema et al. 2011),
`L7F2`–`L7F8` (generated).

## Lecture 20 (slides20.tex), waveform tomography

All extracted with `pdfimages -all` from last year's `old_slides/Lecture_19_-_Waveform_tomography_-_Slides.pdf`
(page given), then cropped with ImageMagick to remove text and maths that had been pasted into the images (those
labels are now typeset on the slide). No originals were downloaded.

* `healing_t1.pdf`, `healing_t2.pdf`, `healing_t3.pdf`: wavefront healing, a plane P-wave passing a small
  low-velocity anomaly, made by `pepcode/scripts/slides_healing.py` with `pepseis.fd2d_psv` (drawn as `L5F1`); they
  replace last year's three snapshots. Last year's weather-forecast aside is dropped.
* `taylor_test.pdf`: a Taylor test of the 2D kernel of `L8F3`/`L8F4`, made by `pepcode/scripts/slides_taylor.py`.
* `bdkernels_p.pdf`, `bdkernels_s.pdf`, `bdkernels_ss.pdf`: cross-correlation traveltime kernels in PREM without
  its ocean, by the adjoint method, from YSpec (`yspec-seismo` 0.1.0): P at 60° in bands about 25 s and 12.5 s
  (ln vp), S and SS at 25 s (ln vs), with each seismogram and its window, on sections along and across the path. Made by `pepcode/scripts/slides_bdkernels.py` (kernels cached in
  `pepcode/data/kernels/`). They replace Liu & Tromp (2008)'s figures on the two slides (the notes keep `L8F5`).

* `farnetani1997_plumes.png` (1193×529; old p. 6): plume snapshots in a 2900 km deep box. Old slide credit
  "Farnetani 1997", probably Farnetani, C. G., 1997. Excess temperature of mantle plumes: the role of chemical
  stratification across D″, *Geophys. Res. Lett.*, 24, 1583–1586 (not checked against the paper). Text lines
  above the plot removed.
* `gradient_descent_contours.png` (932×526; old p. 11): steepest descent vs conjugate gradient on the contours of a
  quadratic. No credit on the old slide; source unknown.
* `local_minima_surface.png` (870×435; old p. 14): misfit surface J(θ0, θ1) with descent paths A–D reaching
  different minima. No credit on the old slide; source unknown (resembles a figure from Andrew Ng's
  machine-learning course).

## Lecture 21 (slides21.tex), self-gravitation and rotation

Both extracted unchanged (pdfimages -all) from last year's deck,
`old_slides/Lecture_20_-_Self-gravitation_-_Slides.pdf`.

* `dehant2015_rotation.jpg` (1000×996), slide 13 of the old deck: schematic of the Earth's rotation axis with
  precession, nutation, the 23.5° obliquity and the tidal torque. Credited "Dehant & Mathews (2015)" (spelt
  "Matthews" on last year's slide): Dehant, V. & Mathews, P. M., 2015. *Precession, Nutation and Wobble of the
  Earth*, Cambridge University Press, ISBN 9781107092549. Publisher copyright; the exact figure number in the
  book has not been checked.

## Lecture 22 (slides22.tex), equilibrium figures

All extracted with `pdfimages -all` from last year's `old_slides/Lecture_21_-_Figure_of_the_Earth_-_Slides.pdf`
(embedded JPEGs at about 150–250 dpi) and cropped of white or black margins and of typed labels/captions
(ImageMagick `-trim`/`-crop`); no other processing. The old slides recorded only the first two credits.

* `chandrasekhar1969_newton_cassini.png` (1630×987, old slide 4): Newton's oblate and Cassini's prolate figures of the
  Earth with Voltaire's couplet, as reproduced in Chandrasekhar, S., 1969. *Ellipsoidal Figures of Equilibrium*, Yale
  University Press. Credit on the slide: "Chandrasekhar (1969)" (as last year).
* `myhill2025_phobos_potential.png` (1219×969, old slide 12): equatorial slice through the gravitational potential of
  Phobos for a constant interior density, panel (c). Credit on the slide: "Myhill et al. (2025)" (as last year).
  Probably Myhill, A. et al., 2025, "Forward and adjoint calculations of gravitational potential in heterogeneous,
  aspherical planets", arXiv:2508.07910 — David to confirm the reference.
* `nasa_newhorizons_pluto.png` (789×777, old slide 24, "A planet"): enhanced-colour image of Pluto, presumably
  NASA/JHUAPL/SwRI (New Horizons, 2015). Not credited last year; slide now says "NASA". David to confirm.
* `nasa_hirise_phobos.png` (860×485, old slide 24, "Not a planet"): Phobos, presumably the 2008 HiRISE image,
  NASA/JPL-Caltech/University of Arizona (Mars Reconnaissance Orbiter). Credit "NASA" as above.
* `haumea_render.jpg` (600×600, old slide 25, unmodified): grey rendering of Haumea as a triaxial ellipsoid.
  Source not recorded; not credited on the slide. David to supply.
* `darwin_fission_illustration.png` (587×878, old slide 26): early twentieth-century book illustration of George
  Darwin's fission sequence ("First a sphere, then an egg, afterwards pearshaped, the stalk-end breaks and forms the
  Moon"), headed "A history of the Moon's formation". Source not recorded; a web search did not identify the book.
  Not credited on the slide.

## Lecture 23 (slides23.tex)

No new borrowed images. Figures used: figures/L10F2.pdf (recap), L11F1.pdf, L11F2.pdf (notes figures);
the compression sketch on the last slide is drawn in TikZ in slides23.tex (redrawn from the old Impress shape).

## Lecture 24 (slides24.tex), mode splitting and coupling

All extracted from last year's deck `old_slides/Lecture_23_-_Mode_coupling_-_Slides.pdf` with `pdfimages -all`
(JPEG as embedded, about 160–190 ppi at the old display size); caption and credit text that Impress had baked into
the bitmaps was cropped off (or painted white), then whitespace trimmed. Original sources as credited on the old
slides; not re-fetched.

* `feynman_lectures_vol1_ch51.jpg` (1334×971, old slide 24): composite of the Benioff et al. (1961) spectrum as
  reproduced in Feynman, R. P., Leighton, R. B. & Sands, M., 1963. *The Feynman Lectures on Physics*, Vol. I,
  Lecture 51 "Waves" (Fig. 51-8 and the accompanying paragraph on the split ellipsoidal mode), plus a photograph
  of Richard Feynman at a blackboard (source not recorded). Slide credit: "The Feynman Lectures on Physics,
  Volume 1, Lecture 51: Waves in Solids" (as last year).
* `iris_travel_time_stack_iasp91.jpg` (1389×814, old slide 25): left, IRIS global stack of three-component
  long-period and broadband seismograms (time against distance); right, IASP91 travel-time curves for a surface
  focus (Kennett & Engdahl 1991) with observed high-frequency phases. IRIS educational material; no credit was on
  the old slide beyond the IRIS logo inside the image.
* `widmer2007_bfo_sumatra_spectrum.jpg` (1284×606, old slide 27) and `widmer2007_measured_eigenfrequencies.jpg`
  (1237×819, old slide 28): Widmer-Schnidrig, R. & Laske, G., 2007. Theory and observations: normal mode and
  surface wave observations, in *Treatise on Geophysics*, Vol. 1, ed. G. Schubert, Elsevier, 67–125. Strain
  spectrum at BFO after the 2004 Sumatra–Andaman earthquake; measured eigenfrequencies against angular order by
  method. Slide credit "Widmer-Schnidrig & Laske (2007)" (last year's slides said "Widmer & Laske").
* `ishii1999_density_model.jpg` (1229×807, old slide 34): Ishii, M. & Tromp, J., 1999. *Science*, 285,
  1231–1236; shear velocity, bulk sound velocity and density at the base of the mantle, total model and degree 2.
  The "Ishii & Tromp (1999)" credit baked into the bitmap was painted out; the panel labels were kept.
* `lau2017_tidal_tomography.jpg` (1183×590, old slide 35): Lau, H. C. P. et al., 2017. Tidal tomography
  constrains Earth's deep-mantle buoyancy, *Nature*, 551, 321–326 (map of excess shear wave speed with GPS
  tidal residuals). Last year's slide credited it as "Lau et al. (2018)"; the reference list gives 2017.
* `garnero2008_lower_mantle_cartoon.jpg` (1306×518, old slide 36): Garnero, E. J. & McNamara, A. K., 2008.
  Structure and dynamics of Earth's lower mantle, *Science*, 320, 626–628 (cartoon cross-section of the lower
  mantle: LLSVP, ULVZ, post-perovskite, D'').
* `adourian2024_spectrum_kernels.png` (1416×1117): Figs 3 and 4 of Adourian, S., Dursun, M. S., Lau, H. C. P. &
  Al-Attar, D., 2024. Adjoint sensitivity kernels for free oscillation spectra, *Geophys. J. Int.*, 238, 257–271,
  doi:10.1093/gji/ggae136, rendered at 300 dpi from the publisher PDF (vector graphics) and stacked, captions
  removed: kernels of a spectrum in a narrow band about ₀S₂ (top) and ₂S₂ (bottom) for density, vs and vp.
  CC BY 4.0. Credit "Adourian et al. (2024)".

## Movies for the live slides (`figures/anim/*.webm`)

`L5F1.webm`, `L5F2.webm` (Lecture 17), `boundary.webm`, `anisotropy_iso.webm`, `anisotropy_ti.webm` (Lecture 16) and
`healing.webm` and `L8F3.webm` (Lecture 20, the adjoint construction of a 2-D kernel): the finite-difference runs behind the stills, written by the same scripts with
`--webm` (`pepseis.fd2d_psv.movie`). Opened from `slidesNN-live.pdf` by clicking the still.
