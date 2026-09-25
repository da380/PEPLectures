"""Lecture 16, Fig. 1: plane sections through the slowness surfaces of two anisotropic media.
Left: a transversely isotropic medium (Love moduli A, C, F, L, N; section containing the
symmetry axis).  Right: single-crystal olivine (orthorhombic; section in the [100]-[001] plane).
The three sheets (qP, qS1, qS2) are shown in different colours."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})
from pepseis.elastic import voigt_to_tensor, ti_voigt, olivine_voigt, slowness_section  # noqa: E402
from pepseis.paths import FIG

# --- media ------------------------------------------------------------------
# (i) a transversely isotropic medium: a strongly anisotropic, PREM-like upper mantle
#     rock with the Love moduli chosen so that alpha_h > alpha_v, beta_h > beta_v.
rho_ti = 3300.0                                             # kg m^-3
A_, C_, F_, L_, N_ = 240.0, 200.0, 70.0, 60.0, 72.0        # GPa
A_ti = voigt_to_tensor(ti_voigt(A_, C_, F_, L_, N_))
# (ii) olivine
rho_ol = 3355.0
A_ol = voigt_to_tensor(olivine_voigt())

e1, e3 = np.eye(3)[0], np.eye(3)[2]
fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.1))
colours = ["tab:blue", "tab:red", "tab:green"]
labels = ["fastest sheet", "intermediate sheet", "slowest sheet"]
for ax, (A, rho, title) in zip(axes, [(A_ti, rho_ti, "transversely isotropic"), (A_ol, rho_ol, "olivine")]):
    phi, s = slowness_section(A, rho, e1, e3)
    for k, (col, lab) in enumerate(zip(colours, labels)):
        # sheet index: s[2] slowest speed = outermost slowness
        r = s[2 - k]
        ax.plot(r * np.cos(phi), r * np.sin(phi), color=col, lw=1.3, label=lab)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$p_1$ / s km$^{-1}$")
    ax.set_ylabel(r"$p_3$ / s km$^{-1}$")
    ax.set_xlim(-0.26, 0.26); ax.set_ylim(-0.26, 0.26)
    ax.set_xticks([-0.2, -0.1, 0, 0.1, 0.2]); ax.set_yticks([-0.2, -0.1, 0, 0.1, 0.2])
    ax.grid(True, ls=":", lw=0.5)
    ax.set_title(title, fontsize=11)
axes[0].legend(loc="upper right", fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig(str(FIG / "L4F1.pdf"))
fig.savefig(str(FIG / "L4F1.png"), dpi=200)
