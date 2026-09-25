"""Lecture 19, Fig. 7: Bayesian solution of the toy problem with a Gaussian prior having a
correlation length of 0.15 (box units): the input model, the posterior mean and the
posterior standard deviation."""
import numpy as np
from pepseis.toytomo import ToyTomography, figure_style
from pepseis.paths import FIG
plt = figure_style()
T = ToyTomography()
Cm = T.prior_covariance(std=0.003, corr_len=0.15)
m_p, sd = T.posterior(Cm)
fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.7), constrained_layout=True)
im = T.draw(axes[0], T.m_true, vmax=0.004, title="input model")
T.draw(axes[1], m_p, vmax=0.004, rays=False, title="posterior mean")
cb = fig.colorbar(im, ax=axes[:2], fraction=0.03, pad=0.02, shrink=0.85); cb.set_label(r"$\delta p/p$")
im2 = T.draw(axes[2], sd, vmin=0, vmax=0.003, rays=False, cmap="viridis", title="posterior standard deviation")
cb2 = fig.colorbar(im2, ax=axes[2], fraction=0.05, pad=0.02, shrink=0.85); cb2.set_label(r"$\sigma(\delta p/p)$")
fig.savefig(str(FIG / "L7F7.pdf"), dpi=200); fig.savefig(str(FIG / "L7F7.png"), dpi=200)
print("prior std 0.003; posterior std min %.5f max %.5f" % (sd.min(), sd.max()))
