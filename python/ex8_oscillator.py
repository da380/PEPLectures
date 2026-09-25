"""Figure for Example 2 of examples8.tex: forward and adjoint fields for the damped
oscillator x'' + 2 gamma x' + omega0^2 x = f(t) with f(t) = exp[-(t-2)^2],
omega0 = 2, gamma = 0.1, T = 10; x_obs generated with omega0 = 2.2, gamma = 0.15.

Saves examples/figures_ex/ex8_oscillator.png.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

here = os.path.dirname(os.path.abspath(__file__))
outdir = os.path.join(here, "..", "examples", "figures_ex")
os.makedirs(outdir, exist_ok=True)

T = 10.0
omega0, gamma = 2.0, 0.1
omega0_true, gamma_true = 2.2, 0.15
opts = dict(method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True)


def force(t):
    return np.exp(-(t - 2.0) ** 2)


def forward(om2, gam):
    def rhs(t, y):
        return [y[1], force(t) - 2 * gam * y[1] - om2 * y[0]]
    return solve_ivp(rhs, (0.0, T), [0.0, 0.0], **opts)


sol_obs = forward(omega0_true ** 2, gamma_true)
sol = forward(omega0 ** 2, gamma)
x = lambda t: sol.sol(t)[0]
x_obs = lambda t: sol_obs.sol(t)[0]


def rhs_adj(t, y):
    return [y[1], (x(t) - x_obs(t)) + 2 * gamma * y[1] - omega0 ** 2 * y[0]]


sol_adj = solve_ivp(rhs_adj, (T, 0.0), [0.0, 0.0], **opts)

t = np.linspace(0, T, 1001)
fig, axes = plt.subplots(2, 1, figsize=(7.5, 5.2), sharex=True, facecolor="white")
axes[0].plot(t, force(t), color="0.6", lw=1.0, label=r"$f(t)$")
axes[0].plot(t, x(t), color="C0", lw=1.4, label=r"$x(t)$")
axes[0].plot(t, x_obs(t), color="k", lw=1.0, ls="--", label=r"$x^{\mathrm{obs}}(t)$")
axes[0].set_ylabel("displacement (dimensionless)")
axes[0].legend(frameon=False, fontsize=9, loc="upper right", ncol=3)
axes[1].plot(t, sol_adj.sol(t)[0], color="C3", lw=1.4, label=r"$x'(t)$")
axes[1].plot(t, x(t) - x_obs(t), color="0.5", lw=1.0, ls=":", label=r"$x-x^{\mathrm{obs}}$")
axes[1].axhline(0, color="0.8", lw=0.6)
axes[1].set_xlabel(r"$t$ (dimensionless)")
axes[1].set_ylabel("adjoint field (dimensionless)")
axes[1].legend(frameon=False, fontsize=9, loc="upper right", ncol=2)
axes[1].set_xlim(0, T)
fig.tight_layout()
fig.savefig(os.path.join(outdir, "ex8_oscillator.png"), dpi=200, facecolor="white")
print("saved", os.path.join(outdir, "ex8_oscillator.png"))
