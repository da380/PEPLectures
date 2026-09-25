"""Lecture 19, Fig. 8: a function-space Bayesian tomography example using pygeoinf (David's
package; this follows its Tutorial 10 for the planar geometry).  The model is a function in a
Sobolev space on a plane with a heat-kernel Gaussian prior, the data are path averages along
straight rays with Gaussian errors, and the posterior is computed with a preconditioned
conjugate-gradient solver.  Run with the pygeoinf virtual environment, e.g.
    <venv>/bin/python L7F8.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pygeoinf as inf
from pygeoinf.symmetric_space.plane import Sobolev, plot, plot_geodesic_network
from pepseis.paths import FIG

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})
np.random.seed(42)
order, scale, prior_scale = 2, 0.2, 0.1
n_sources, n_receivers = 12, 12

model_space = Sobolev.from_heat_kernel_prior(prior_scale, order, scale, ax=0.0, bx=6.0, cx=0.5,
                                             ay=0.0, by=6.0, cy=0.5, power_of_two=True, min_degree=32)
receivers = model_space.random_points(n_receivers); sources = model_space.random_points(n_sources)
paths = [(src, rec) for src in sources for rec in receivers]
forward_operator = model_space.path_average_operator(paths)
data_error_measure = inf.GaussianMeasure.from_standard_deviation(forward_operator.codomain, 0.05)
forward_problem = inf.LinearForwardProblem(forward_operator, data_error_measure=data_error_measure)
model_prior = model_space.point_value_scaled_heat_kernel_gaussian_measure(prior_scale)
model_true, data_obs = forward_problem.joint_measure(model_prior).sample()

inverse_problem = inf.LinearBayesianInversion(forward_problem, model_prior)
surrogate_space = model_space.with_degree(model_space.degree // 4)
raw_surrogate_prior = surrogate_space.point_value_scaled_heat_kernel_gaussian_measure(prior_scale)
woodbury_solver = inf.CholeskySolver(galerkin=True)
damped_surrogate_prior = raw_surrogate_prior.with_regularized_inverse(woodbury_solver, damping=1e-6)
precon = inverse_problem.surrogate_inversion(
    alternate_forward_operator=surrogate_space.path_average_operator(paths),
    alternate_prior_measure=damped_surrogate_prior).woodbury_data_preconditioner(woodbury_solver)
solver = inf.CGMatrixSolver()
model_posterior = inverse_problem.model_posterior_measure(data_obs, solver, preconditioner=precon)
print("CG iterations:", solver.iterations, "; model space dimension:", model_space.dim)
posterior_std = model_posterior.sample_pointwise_std(200)

fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.4), constrained_layout=True)
vmax = float(max(np.max(np.abs(model_true.data)), np.max(np.abs(model_posterior.expectation.data))))
src_style = {"marker": "*", "color": "black", "s": 60, "edgecolor": "black"}      # as in the toy-problem figures
rec_style = {"marker": "v", "color": "black", "s": 35, "edgecolor": "black"}
for ax, field, title in [(axes[0], model_true, "true model"), (axes[1], model_posterior.expectation, "posterior expectation")]:
    _, im = plot(model_space, field, ax=ax, colorbar=False, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    plot_geodesic_network(paths, ax=ax, alpha=0.3, linewidth=0.4, color="black", source_kwargs=src_style, receiver_kwargs=rec_style)
    ax.set_title(title, fontsize=11); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
cb = fig.colorbar(im, ax=axes[:2], orientation="horizontal", shrink=0.5, pad=0.04); cb.set_label("model value")
_, im2 = plot(model_space, posterior_std, ax=axes[2], colorbar=False, cmap="viridis")
cb2 = fig.colorbar(im2, ax=axes[2], orientation="horizontal", shrink=0.8, pad=0.04); cb2.set_label("standard deviation")
axes[2].set_title("posterior standard deviation", fontsize=11); axes[2].set_aspect("equal"); axes[2].set_xticks([]); axes[2].set_yticks([])
for ax in axes:                                   # avoid hairline seams between cells in PDF viewers
    for coll in ax.collections:
        coll.set_rasterized(True); coll.set_edgecolor("face")
fig.savefig(str(FIG / "L7F8.pdf"), dpi=200); fig.savefig(str(FIG / "L7F8.png"), dpi=200)
