"""
Slowness curves in the (p_1, p_3)-plane for the sample transversely isotropic
medium of Example 3 (worked examples for Lecture 15).  The three sheets are
obtained from a numerical eigen-decomposition of the Christoffel matrix; the
analytic ellipse for the SH sheet, mu p_1^2 + (mu - zeta) p_3^2 = rho, is
overlaid as a dashed curve.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d3 = np.eye(3)


def isotropic(lam, mu):
    return (lam * np.einsum("ij,kl->ijkl", d3, d3)
            + mu * (np.einsum("ik,jl->ijkl", d3, d3) + np.einsum("il,jk->ijkl", d3, d3)))


def transversely_isotropic(lam, mu, gam, xi, zeta, nu):
    A = isotropic(lam, mu)
    A += 8 * gam * np.einsum("i,j,k,l->ijkl", nu, nu, nu, nu)
    A += 4 * xi * (np.einsum("i,j,kl->ijkl", nu, nu, d3) + np.einsum("ij,k,l->ijkl", d3, nu, nu))
    A -= zeta * (np.einsum("i,k,jl->ijkl", nu, nu, d3) + np.einsum("j,k,il->ijkl", nu, nu, d3)
                 + np.einsum("j,l,ik->ijkl", nu, nu, d3) + np.einsum("i,l,jk->ijkl", nu, nu, d3))
    return A


lam, mu, gam, xi, zeta, rho = 80e9, 70e9, 2e9, 1e9, 8e9, 3300.0
nu = np.array([0.0, 0.0, 1.0])
A = transversely_isotropic(lam, mu, gam, xi, zeta, nu)

thetas = np.linspace(0, 2 * np.pi, 721)
# branches identified by polarisation: SH (a parallel to e_2), quasi-SV, quasi-P
e2 = np.array([0.0, 1.0, 0.0])
branches = {"SH": [], "qSV": [], "qP": []}
for th in thetas:
    ph = np.array([np.sin(th), 0.0, np.cos(th)])
    G = np.einsum("ijkl,j,l->ik", A, ph, ph) / rho
    w, v = np.linalg.eigh(G)
    c = np.sqrt(w)
    ish = int(np.argmax(np.abs(v.T @ e2)))          # the e_2-polarised eigenvector
    others = [k for k in range(3) if k != ish]
    iqp = max(others, key=lambda k: w[k])           # faster of the remaining two
    iqsv = min(others, key=lambda k: w[k])
    for name, k in (("SH", ish), ("qSV", iqsv), ("qP", iqp)):
        branches[name].append(ph[[0, 2]] / c[k] * 1e3)   # slowness in s/km
branches = {k: np.array(v) for k, v in branches.items()}

fig, ax = plt.subplots(figsize=(5.5, 6.0), facecolor="white")
style = {"qP": ("tab:green", "quasi-P"), "qSV": ("tab:orange", "quasi-SV"), "SH": ("tab:blue", "SH (numerical)")}
for name in ("qP", "qSV", "SH"):
    col, lab = style[name]
    ax.plot(branches[name][:, 0], branches[name][:, 1], color=col, lw=1.8, label=lab)
# analytic SH ellipse
t = np.linspace(0, 2 * np.pi, 400)
a1 = np.sqrt(rho / mu) * 1e3
a3 = np.sqrt(rho / (mu - zeta)) * 1e3
ax.plot(a1 * np.cos(t), a3 * np.sin(t), "k--", lw=1.0,
        label=r"SH ellipse $\mu p_1^2+(\mu-\zeta)p_3^2=\rho$")
ax.set_aspect("equal")
ax.set_xlabel(r"$p_1$ (s km$^{-1}$)")
ax.set_ylabel(r"$p_3$ (s km$^{-1}$)")
ax.axhline(0, color="0.8", lw=0.5, zorder=0)
ax.axvline(0, color="0.8", lw=0.5, zorder=0)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=8, frameon=False)
lim = 0.25
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "examples", "figures_ex")
os.makedirs(out, exist_ok=True)
fig.savefig(os.path.join(out, "ex4_slowness_TI.png"), dpi=200, facecolor="white")
print("saved", os.path.join(out, "ex4_slowness_TI.png"))
# sanity check: the numerical SH branch lies on the analytic ellipse
p = branches["SH"] / 1e3
print("max |mu p1^2 + (mu - zeta) p3^2 - rho|/rho on the SH branch:",
      np.max(np.abs(mu * p[:, 0]**2 + (mu - zeta) * p[:, 1]**2 - rho)) / rho)
