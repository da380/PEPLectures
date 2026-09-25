"""Lecture 19, Fig. 5: the effect of the damping parameter, with the trade-off curve between
data misfit and model norm."""
import numpy as np
from pepseis.toytomo import ToyTomography, figure_style
from pepseis.paths import FIG
plt = figure_style()
T = ToyTomography(); lam_star = T.lambda_for_chi2()
lams = lam_star * np.logspace(-4, 4, 41)
misfit = []; norm = []
for lam in lams:
    m = T.solve(lam); misfit.append(T.chi2_per_datum(m)); norm.append(np.linalg.norm(m))
fig, axes = plt.subplots(2, 2, figsize=(8.0, 7.6), constrained_layout=True)
for ax, f, lab in [(axes[0, 0], 0.01, r"$\lambda=\lambda_{*}/100$"), (axes[0, 1], 1.0, r"$\lambda=\lambda_{*}$"), (axes[1, 0], 100.0, r"$\lambda=100\lambda_{*}$")]:
    im = T.draw(ax, T.solve(f * lam_star), vmax=0.004, rays=False, title=lab)
ax = axes[1, 1]
ax.loglog(misfit, norm, "k-", lw=1)
for f, mk in [(0.01, "s"), (1.0, "o"), (100.0, "^")]:
    m = T.solve(f * lam_star); ax.loglog(T.chi2_per_datum(m), np.linalg.norm(m), "r" + mk, ms=7)
ax.axvline(1.0, color="0.6", lw=0.6, ls="--")
ax.set_xlabel(r"data misfit $\chi^{2}/n$"); ax.set_ylabel(r"model norm $\|\mathbf{m}\|$"); ax.set_title("trade-off curve", fontsize=11)
cb = fig.colorbar(im, ax=axes[1, 0], fraction=0.05, pad=0.03, shrink=0.8, location="bottom"); cb.set_label(r"$\delta p/p$")
fig.savefig(str(FIG / "L7F5.pdf"), dpi=200); fig.savefig(str(FIG / "L7F5.png"), dpi=200)
