"""Slides, Lecture 20: banana-doughnut kernels in PREM (without its ocean, which adds noise near the surface)
by the adjoint method, from YSpec (yspec-seismo).
Cross-correlation traveltimes at 60 degrees, measured in a window from one dominant period before the arrival
(TauP) to one and a half after, with cosine tapers of half a period: P on the vertical at periods of 25 s and
12.5 s (kernel of ln vp), and S and SS on the transverse at 25 s (ln vs).  Each row shows the seismogram with the
window shaded, and the kernel on the section along the path and across it at its middle, each on its own colour
scale, clipped at the 98th percentile away from the source and the receiver, with blue negative as in L8F3/L8F4.
Kernels and seismograms are cached in data/kernels/; the calculation uses OMP_NUM_THREADS threads, or 4 (the
12.5 s band then takes about twenty minutes).
Written as PDF to ../figures/slides.
Usage:  poetry run python scripts/slides_bdkernels.py"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis.paths import FIG

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "kernels")
os.makedirs(DATA, exist_ok=True)
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 13})
R, R_CMB = 6368.0, 3480.0                        # surface of PREM without its ocean
DISTANCE, AZIMUTH, DEPTH, RECEIVER_DEPTH = 60.0, 70.0, 29.5, 0.0
MOMENT = dict(rr=4.22e17, rt=3.08e17, rp=6.99e17, tt=-2.976e18, tp=6.54e17, pp=2.55e18)   # 2023 Kahramanmaras

# band: period (s), [(phase, component)], calculation, section spacing (degrees, km)
BANDS = {
    "25s": (25.0, [("P", "Z"), ("S", "T"), ("SS", "T")],
            dict(degrees=450, frequencies_mhz=(7, 40), length_min=30, time_step_s=1.2, filter_mhz=(8, 12, 35, 40)),
            (0.2, 25.0)),
    "12s": (12.5, [("P", "Z")],
            dict(degrees=900, frequencies_mhz=(14, 80), length_min=30, time_step_s=0.6, filter_mhz=(15, 25, 70, 80)),
            (0.12, 15.0)),
}
COMPONENT = {"P": "Z", "S": "T", "SS": "T"}


def windows(band):
    """{phase: (start, end, taper)} in s, about the TauP arrival."""
    from obspy.taup import TauPyModel
    period, phases = BANDS[band][:2]
    taup = TauPyModel("prem")
    out = {}
    for phase, _ in phases:
        t = min(a.time for a in taup.get_travel_times(DEPTH, DISTANCE, [phase]))
        out[phase] = (round(t - period), round(t + 1.5 * period), period / 2)
    return out


def _setup(band):
    import yspec
    yspec.set_threads(int(os.environ.get("OMP_NUM_THREADS", 4)))
    model = yspec.Model.prem("prem-no-ocean")
    source = yspec.Source(latitude=36.87, longitude=35.58, depth_km=DEPTH, moment_tensor=MOMENT)
    receiver = yspec.great_circle(source, azimuth_deg=AZIMUTH, distances_deg=[DISTANCE], depth_km=RECEIVER_DEPTH)[0]
    calculation = yspec.Calculation(quantity="velocity", **BANDS[band][2])
    return yspec, model, source, receiver, calculation


def kernels(band):
    """Kernels on the two sections, cached: {phase_section: values, section_coordinate: ...}."""
    phases = BANDS[band][1]
    path = os.path.join(DATA, f"bdkernels_{band}_{'-'.join(p for p, _ in phases)}.npz")
    if os.path.exists(path):
        return dict(np.load(path))
    yspec, model, source, receiver, calculation = _setup(band)
    spacing, radial = BANDS[band][3]
    w = windows(band)
    measured = [yspec.Traveltime(c, w[p][:2], w[p][2], p) for p, c in phases]
    measurements = yspec.traveltimes(model, calculation, source, receiver, measured)
    options = dict(margin_deg=8.0, half_width_deg=10.0, spacing_deg=spacing, radial_spacing_km=radial)
    sections = [yspec.path_section(model, source, receiver, **options),
                yspec.cross_section(model, source, receiver, **options)]
    result = yspec.sampled_kernels(model, calculation, measurements, sections, parameters="velocities",
                                   scaling="relative", progress=True)
    out = {}
    for (phase, _), row in zip(phases, result):
        for name, k in zip(("path", "across"), row):
            out[f"{phase}_{name}"] = k["vp" if phase.startswith("P") else "vs"]
            out[f"{name}_radius"], out[f"{name}_layers"] = k.radius_km, k.layers
            out[f"{name}_angle"] = k.along_deg if name == "path" else k.across_deg
    np.savez_compressed(path, **out)
    return out


def traces(band):
    """The vertical and transverse velocity seismograms of the band, cached."""
    path = os.path.join(DATA, f"traces_{band}.npz")
    if os.path.exists(path):
        return dict(np.load(path))
    yspec, model, source, receiver, calculation = _setup(band)
    s = yspec.seismograms(model, calculation, source, receiver)
    out = dict(time=s.time_s, Z=s.trace("Z"), T=s.trace("T"))
    np.savez_compressed(path, **out)
    return out


def trace(ax, band, phase, text):
    """The seismogram measured, four minutes either side of the window, scaled to the phase in the window, with
    the window and its tapers shaded and the part inside it drawn in black."""
    d = traces(band)
    start, end, taper = windows(band)[phase]
    t, u = d["time"], d[COMPONENT[phase]]
    w = np.clip(np.minimum(t - start, end - t) / taper, 0, 1)
    w = np.where(w < 1, 0.5 * (1 - np.cos(np.pi * w)), 1.0)
    u = 0.8 * u / np.abs(u[(t > start) & (t < end)]).max()
    ax.fill_between(t / 60, -w, w, color="0.88", lw=0, zorder=0)
    ax.plot(t / 60, u, color="0.6", lw=0.6)
    inside = w > 0
    ax.plot(t[inside] / 60, u[inside], color="k", lw=1.0)
    ax.set_xlim(start / 60 - 4, end / 60 + 4); ax.set_ylim(-1.05, 1.05); ax.set_yticks([])
    ax.set_xlabel("time / minutes")
    ax.text(0.03, 0.9, text, transform=ax.transAxes, va="top")


def section(ax, d, phase, name):
    """One section, drawn in its own plane, mantle only, with the source (star) and receiver (triangle)."""
    r, layers, angle = d[f"{name}_radius"], d[f"{name}_layers"], d[f"{name}_angle"]
    K = d[f"{phase}_{name}"]
    centre = 0.5 * DISTANCE if name == "path" else 0.0
    theta = np.radians(angle - centre)
    far = (np.abs(angle) > 10) & (np.abs(angle - DISTANCE) > 10) if name == "path" else np.ones(angle.size, bool)
    limit = np.percentile(np.abs(K[(r > R_CMB) & (R - r > 150)][:, far]), 98)
    for layer in np.unique(layers):
        i = (layers == layer) & (r >= R_CMB)
        if i.sum() < 2:
            continue
        o = np.argsort(r[i], kind="stable")
        rr, KK = r[i][o], K[i][o]
        ax.pcolormesh(rr[:, None] * np.sin(theta), rr[:, None] * np.cos(theta), KK / limit, cmap="RdBu_r",
                      vmin=-1, vmax=1, shading="gouraud", rasterized=True)
    arc = np.linspace(theta[0], theta[-1], 400)
    for rb in (R, R_CMB):
        ax.plot(rb * np.sin(arc), rb * np.cos(arc), "k", lw=0.7)
    for side in (theta[0], theta[-1]):
        ax.plot([R_CMB * np.sin(side), R * np.sin(side)], [R_CMB * np.cos(side), R * np.cos(side)], "k", lw=0.7)
    if name == "path":
        for a, depth, mark in ((0.0, DEPTH, "*"), (DISTANCE, 0.0, "v")):
            t = np.radians(a - centre)
            ax.plot((R - depth) * np.sin(t), (R - depth) * np.cos(t), mark, ms=13 if mark == "*" else 10,
                    mfc="gold" if mark == "*" else "w", mec="k", zorder=5, clip_on=False)
    ax.set_aspect("equal"); ax.axis("off")


def figure(rows, label, fname):
    """rows: [(band, phase, label on the trace)]."""
    fig, axes = plt.subplots(len(rows), 3, figsize=(15, 3.4 * len(rows) + 0.4), squeeze=False,
                             gridspec_kw={"width_ratios": [1.7, 3.2, 1]})
    for (band, phase, text), (a0, a1, a2) in zip(rows, axes):
        trace(a0, band, phase, text)
        d = kernels(band)
        section(a1, d, phase, "path"); section(a2, d, phase, "across")
    axes[0, 1].set_title("along the path", fontsize=13); axes[0, 2].set_title("across it, at the middle", fontsize=13)
    sm = plt.cm.ScalarMappable(norm=plt.Normalize(-1, 1), cmap="RdBu_r")
    cb = fig.colorbar(sm, ax=axes[:, 1:].ravel().tolist(), fraction=0.02, pad=0.02, shrink=0.8, ticks=[-1, 0, 1])
    cb.set_label(label)
    fig.savefig(str(FIG / "slides" / fname), bbox_inches="tight")
    plt.close(fig)


figure([("25s", "P", "P on Z, 25 s"), ("12s", "P", "P on Z, 12.5 s")], r"$K_{\alpha}$, scaled", "bdkernels_p.pdf")
figure([("25s", "S", "S on T, 25 s")], r"$K_{\beta}$, scaled", "bdkernels_s.pdf")
figure([("25s", "SS", "SS on T, 25 s")], r"$K_{\beta}$, scaled", "bdkernels_ss.pdf")
