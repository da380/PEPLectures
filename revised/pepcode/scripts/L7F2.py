"""Lecture 19, Fig. 2: input model with rays, and the regularised least-squares recovery."""
from pepseis.toytomo import ToyTomography, figure_style
from pepseis.paths import FIG
plt = figure_style()
T = ToyTomography(); lam = T.lambda_for_chi2()
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.9), constrained_layout=True)
im = T.draw(axes[0], T.m_true, vmax=0.004, title="input model")
T.draw(axes[1], T.solve(lam), vmax=0.004, title="recovered model")
cb = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.03, shrink=0.85); cb.set_label(r"slowness perturbation $\delta p/p$")
fig.savefig(str(FIG / "L7F2.pdf"), dpi=200); fig.savefig(str(FIG / "L7F2.png"), dpi=200)
