"""Lecture 19, Fig. 3: a model in the null space of the toy problem, and the input model with
this null-space model added; the two models produce identical data."""
import numpy as np
from pepseis.toytomo import ToyTomography, figure_style
from pepseis.paths import FIG
plt = figure_style()
T = ToyTomography()
rng = np.random.default_rng(11)
n = T.n
kx = 2 * np.pi * np.fft.fftfreq(n, 1 / n); KX, KY = np.meshgrid(kx, kx, indexing="ij")
field = np.real(np.fft.ifft2(np.fft.fft2(rng.standard_normal((n, n))) * np.exp(-(KX ** 2 + KY ** 2) * 0.06 ** 2 / 4)))
m_null = T.null_space_component(0.004 * field / np.abs(field).max())
m_null *= 0.004 / np.abs(m_null).max()
print("max |A m_null| / sigma =", np.abs(T.A @ m_null.ravel()).max() / T.sigma)
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.9), constrained_layout=True)
im = T.draw(axes[0], m_null, vmax=0.004, title="a model in the null space")
T.draw(axes[1], T.m_true + m_null, vmax=0.004, title="input model plus null-space model")
cb = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.03, shrink=0.85); cb.set_label(r"slowness perturbation $\delta p/p$")
fig.savefig(str(FIG / "L7F3.pdf"), dpi=200); fig.savefig(str(FIG / "L7F3.png"), dpi=200)
