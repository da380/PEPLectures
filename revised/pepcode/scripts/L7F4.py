"""Lecture 19, Fig. 4: singular values of the toy problem, and two singular vectors: one that
is well constrained by the data and one (chosen to be localised in a few cells) that is barely constrained."""
import numpy as np
from pepseis.toytomo import ToyTomography, figure_style
from pepseis.paths import FIG
plt = figure_style()
T = ToyTomography()
s2, V = np.linalg.eigh(T.AtA * T.sigma ** 2)          # eigen-decomposition of A^T A
order = np.argsort(s2)[::-1]; s = np.sqrt(np.clip(s2[order], 0, None)); V = V[:, order]
n_data = T.A.shape[0]
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), constrained_layout=True)
ax = axes[0]
ax.semilogy(np.arange(1, n_data + 1), s[:n_data] / s[0], "k.", ms=3)
ax.set_xlabel(r"index $j$"); ax.set_ylabel(r"$s_{j}/s_{1}$"); ax.set_title("singular values", fontsize=11)
ax.set_xlim(0, n_data + 5); ax.grid(True, ls=":", lw=0.5)
for ax, j in [(axes[1], 0), (axes[2], 277)]:   # the 278th of at most 280 non-zero singular values: s_j/s_1 ~ 7e-4
    v = V[:, j].reshape(T.n, T.n); v /= np.abs(v).max()
    T.draw(ax, v, vmax=1, rays=False, title=r"singular vector $j=%d$" % (j + 1))
fig.savefig(str(FIG / "L7F4.pdf"), dpi=200); fig.savefig(str(FIG / "L7F4.png"), dpi=200)
print("rank (s > 1e-10 s1):", np.sum(s > 1e-10 * s[0]), "of", T.n ** 2, "; number of data", n_data)
