"""Lecture 19, Fig. 6: spike tests.  The regularised least-squares recovery of a model that is
a single spike (noise-free data) is the corresponding column of the resolution matrix.  Input
and output are shown side by side for a spike in a well-covered block and in a poorly covered one."""
import numpy as np
from pepseis.toytomo import ToyTomography, figure_style
from pepseis.paths import FIG
plt = figure_style()
T = ToyTomography(); lam = T.lambda_for_chi2(); n = T.n
hits = (T.A != 0).sum(axis=0).reshape(n, n)             # number of rays through each block
good = np.unravel_index(np.argmax(hits), hits.shape)
poor_candidates = np.argwhere((hits >= 3) & (hits <= 4))
poor = tuple(poor_candidates[np.argmax(np.abs(poor_candidates - np.array(good)).sum(axis=1))])
print("well-covered block", good, "rays:", hits[good], "; poorly covered block", poor, "rays:", hits[poor])
fig, axes = plt.subplots(2, 2, figsize=(8.0, 7.6), constrained_layout=True)
for row, (cell, name) in enumerate([(good, "well-covered block"), (poor, "poorly covered block")]):
    m_in = np.zeros((n, n)); m_in[cell] = 1.0                # unit spike
    m_out = T.solve(lam, d=T.A @ m_in.ravel())                # noise-free synthetic data
    print(name, ": peak recovered / input = %.2f" % m_out.max())
    for col, (m, what) in enumerate([(m_in, "input"), (m_out, "recovered")]):
        ax = axes[row, col]
        im = T.draw(ax, m, vmax=1.0, rays=False, title="%s: %s" % (name, what))
        for src, rec in T.rays:
            ax.plot([src[0], rec[0]], [src[1], rec[1]], "k", lw=0.2, alpha=0.25)
cb = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.03, shrink=0.8); cb.set_label("amplitude relative to the input spike")
fig.savefig(str(FIG / "L7F6.pdf"), dpi=200); fig.savefig(str(FIG / "L7F6.png"), dpi=200)
