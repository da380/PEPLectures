"""Slides, Lecture 18: ray paths (left) and travel-time curves (right) in PREM for a surface source, for the
core phases PKKKKP, PcP and PcS, and the core-diffracted Pdiff and Sdiff.  Uses obspy.taup (which handles the
diffracted phases); drawn in the style of L6F3 and L6F4.  Written as PDF to ../figures/slides.
Usage:  poetry run python scripts/slides_taup.py"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
from pepseis import prem
from pepseis.paths import FIG

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 13})   # larger than the notes figures: read at slide scale
MODEL = TauPyModel("prem")
R = 6371.0

# figure name: [(phase, colour, epicentral angles of the rays drawn)]
FIGURES = {
    "rays_pkkkkp": [("PKKKKP", "tab:red", np.arange(0, 181, 6))],
    "rays_pcp_pcs": [("PcP", "tab:red", np.arange(10, 91, 10)), ("PcS", "tab:blue", np.arange(10, 61, 10))],
    "rays_pdiff_sdiff": [("Pdiff", "tab:red", np.arange(100, 161, 10)), ("Sdiff", "tab:blue", np.arange(100, 161, 10))],
}


def globe(ax):
    for r, col in [(R, "0.85"), (prem.R_CMB, "0.72"), (prem.R_ICB, "0.6")]:
        ax.add_patch(plt.Circle((0, 0), r, color=col, zorder=0))
        ax.add_patch(plt.Circle((0, 0), r, fill=False, color="k", lw=0.7, zorder=1))
    for d in range(0, 181, 30):
        t = np.radians(d)
        ax.text(1.08 * R * np.sin(t), 1.08 * R * np.cos(t), r"$%d^\circ$" % d, ha="center", va="center", fontsize=11)
    ax.plot([0], [R], "k*", ms=11, zorder=5)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-1.15 * R, 1.15 * R); ax.set_ylim(-1.15 * R, 1.15 * R)


def curve(phase):
    """(Delta, T) along the travel-time curve, ordered by ray parameter so that branches join up."""
    rows = []
    for d in np.arange(0.0, 180.01, 0.5):
        for a in MODEL.get_travel_times(0.0, d, phase_list=[phase]):
            rows.append((a.ray_param, d, a.time))
    rows.sort()
    return np.array(rows)[:, 1:].T if rows else np.zeros((2, 0))


for name, phases in FIGURES.items():
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9.5, 4.5), gridspec_kw={"width_ratios": [1, 1.05]})
    globe(ax0)
    for phase, col, dists in phases:
        for d in dists:
            for a in MODEL.get_ray_paths(0.0, d, phase_list=[phase]):
                r = R - a.path["depth"]; D = a.path["dist"]
                for s in (1, -1):                      # both sides of the source
                    ax0.plot(s * r * np.sin(D), r * np.cos(D), color=col, lw=0.6, zorder=2)
        d, t = curve(phase)
        jumps = np.where(np.abs(np.diff(d)) > 2)[0]
        for k, seg in enumerate(np.split(np.arange(d.size), jumps + 1)):
            ax1.plot(d[seg], t[seg] / 60, color=col, lw=1.5, label=phase if k == 0 else None)
    ax1.set_xlim(0, 180); ax1.set_xticks(range(0, 181, 30))
    ax1.set_xlabel(r"epicentral angle $\Delta$ / degrees"); ax1.set_ylabel("travel time / minutes")
    ax1.grid(True, ls=":", lw=0.5)
    ax1.legend(loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(str(FIG / "slides" / (name + ".pdf")))
    plt.close(fig)
