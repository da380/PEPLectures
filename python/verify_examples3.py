"""Numerical checks for examples/examples3.tex (Worked examples for Lecture 14: Seismic sources).

Every numerical or symbolic-by-numbers claim in the solutions is checked below and reported
as an OK/FAIL line.  Only numpy and the Python standard library are used.
"""
import numpy as np
from fractions import Fraction
from itertools import permutations

rng = np.random.default_rng(20260903)
failures = 0


def report(name, ok, detail=""):
    global failures
    if not ok:
        failures += 1
    print(f"{'OK  ' if ok else 'FAIL'} {name}" + (f"  [{detail}]" if detail else ""))


def isotropic_A(lam, mu):
    d = np.eye(3)
    return (lam * np.einsum("ij,kl->ijkl", d, d)
            + mu * (np.einsum("ik,jl->ijkl", d, d) + np.einsum("il,jk->ijkl", d, d)))


def levi_civita():
    eps = np.zeros((3, 3, 3))
    for p in permutations(range(3)):
        i, j, k = p
        eps[i, j, k] = np.linalg.det(np.eye(3)[list(p)])
    return eps


eps = levi_civita()
e1, e2, e3 = np.eye(3)

# ----------------------------------------------------------------------------------------
print("\nExample 1: stress glut in an isotropic body, moment tensor, double couple")
lam, mu = rng.uniform(1, 5, 2)
A = isotropic_A(lam, mu)
w = np.array([rng.normal(), rng.normal(), 0.0])          # tangential slip, w3 = 0
lhs = np.einsum("ijk,k->ij", A[:, :, :, 2], w)          # A_{ijk3} w_k
rhs = mu * (np.outer(w, e3) + np.outer(e3, w))          # mu (w_i d_j3 + w_j d_i3)
report("A_ijk3 w_k = mu(w_i d_j3 + w_j d_i3) for w3=0", np.allclose(lhs, rhs),
       f"max diff {np.abs(lhs - rhs).max():.2e}")
report("glut is symmetric", np.allclose(lhs, lhs.T))
report("glut is traceless", abs(np.trace(lhs)) < 1e-12, f"trace {np.trace(lhs):.2e}")
# with an opening component the lambda term appears: A_ijk3 w_k = lam w3 d_ij + mu(...)
w_open = rng.normal(size=3)
lhs_o = np.einsum("ijk,k->ij", A[:, :, :, 2], w_open)
rhs_o = lam * w_open[2] * np.eye(3) + mu * (np.outer(w_open, e3) + np.outer(e3, w_open))
report("general w: A_ijk3 w_k = lam w3 d_ij + mu(w_i d_j3 + w_j d_i3)", np.allclose(lhs_o, rhs_o))
# general normal n, tangential w: A_ijkl w_k n_l = mu (w_i n_j + w_j n_i)
n = rng.normal(size=3); n /= np.linalg.norm(n)
wt = rng.normal(size=3); wt -= (wt @ n) * n
lhs_n = np.einsum("ijkl,k,l->ij", A, wt, n)
report("general normal: A_ijkl w_k n_l = mu(w_i n_j + w_j n_i)",
       np.allclose(lhs_n, mu * (np.outer(wt, n) + np.outer(n, wt))))

# double couple
M0 = mu * 3.7 * 1.3                                     # mu * area * D (arbitrary numbers)
M = M0 * (np.outer(e1, e3) + np.outer(e3, e1))
evals, evecs = np.linalg.eigh(M)
report("double-couple eigenvalues (-M0, 0, M0)", np.allclose(np.sort(evals), [-M0, 0, M0]),
       f"eigenvalues/M0 = {np.sort(evals)/M0}")
t_hat = (e1 + e3) / np.sqrt(2); p_hat = (e1 - e3) / np.sqrt(2); b_hat = e2
report("M t = +M0 t", np.allclose(M @ t_hat, M0 * t_hat))
report("M p = -M0 p", np.allclose(M @ p_hat, -M0 * p_hat))
report("M b = 0", np.allclose(M @ b_hat, 0))
report("M = M0 (t t^T - p p^T)", np.allclose(M, M0 * (np.outer(t_hat, t_hat) - np.outer(p_hat, p_hat))))
report("trace M = 0", abs(np.trace(M)) < 1e-12)
report("t, p at 45 degrees to fault plane",
       np.isclose(np.degrees(np.arcsin(abs(t_hat @ e3))), 45) and np.isclose(np.degrees(np.arcsin(abs(p_hat @ e3))), 45))
# auxiliary plane: slip along e3 on a fault with normal e1 gives the same M
M_aux = mu * 3.7 * 1.3 * (np.outer(e3, e1) + np.outer(e1, e3))
report("auxiliary plane gives identical moment tensor", np.allclose(M, M_aux))
# characteristic polynomial -sigma (sigma^2 - M0^2)
sig = rng.normal() * M0
report("det(M - sigma I) = -sigma(sigma^2 - M0^2)",
       np.isclose(np.linalg.det(M - sig * np.eye(3)), -sig * (sig**2 - M0**2)))

# ----------------------------------------------------------------------------------------
print("\nExample 2: equivalent body force, net force and torque, force couples")
# Smooth (Gaussian) delta on a grid; f_i = -M_ij d_j delta.  Quadrature of net force and torque.
Lbox, N = 6.0, 121
x = np.linspace(-Lbox, Lbox, N); dx = x[1] - x[0]
X = np.stack(np.meshgrid(x, x, x, indexing="ij"))
s = 0.5
delta_s = np.exp(-np.sum(X**2, axis=0) / (2 * s**2)) / (2 * np.pi * s**2) ** 1.5
grad_delta = -X / s**2 * delta_s                        # analytic gradient of the Gaussian
report("smoothed delta integrates to 1", np.isclose(delta_s.sum() * dx**3, 1.0, rtol=1e-6),
       f"{delta_s.sum()*dx**3:.8f}")
M_sym = rng.normal(size=(3, 3)); M_sym = M_sym + M_sym.T
f = -np.einsum("ij,jxyz->ixyz", M_sym, grad_delta)
net_force = f.sum(axis=(1, 2, 3)) * dx**3
torque = np.einsum("ijk,jxyz,kxyz->i", eps, X, f) * dx**3
report("net force of point source vanishes", np.allclose(net_force, 0, atol=1e-8), f"{net_force}")
report("net torque vanishes for symmetric M", np.allclose(torque, 0, atol=1e-8), f"{torque}")
# with an antisymmetric part the torque is eps_ijk M_kj (check formula, not zero)
M_gen = rng.normal(size=(3, 3))
f_gen = -np.einsum("ij,jxyz->ixyz", M_gen, grad_delta)
torque_gen = np.einsum("ijk,jxyz,kxyz->i", eps, X, f_gen) * dx**3
report("torque = eps_ijk M_kj for general M", np.allclose(torque_gen, np.einsum("ijk,kj->i", eps, M_gen), atol=1e-6),
       f"quadrature {torque_gen}, formula {np.einsum('ijk,kj->i', eps, M_gen)}")
# torque equals minus twice the axial vector of the antisymmetric part
axial = 0.5 * np.einsum("ijk,jk->i", eps, M_gen)         # a_i = (1/2) eps_ijk M_jk
report("torque = -2 x axial vector of antisymmetric part", np.allclose(torque_gen, -2 * axial, atol=1e-6))

# Force couples: torque of each couple in the double couple
h = 1e-3
F = M0 / h
tau1 = np.cross(0.5 * h * e3, F * e1) + np.cross(-0.5 * h * e3, -F * e1)
tau2 = np.cross(0.5 * h * e1, F * e3) + np.cross(-0.5 * h * e1, -F * e3)
report("couple 1 torque = +M0 e2", np.allclose(tau1, M0 * e2), f"{tau1/M0}")
report("couple 2 torque = -M0 e2", np.allclose(tau2, -M0 * e2), f"{tau2/M0}")
report("torques cancel", np.allclose(tau1 + tau2, 0))
# sign of the finite-difference representation of d_3 delta (eq. 16): pair both sides with a
# smooth phi.  <phi, d_3 delta> = -d_3 phi(0), while <phi, delta(x + a)> = phi(-a).
phi = lambda p: np.cos(p[0] + 2 * p[1]) * np.exp(-p[2])
d3phi0 = -1.0                                             # d_3 phi at the origin
hh = 1e-4
pairing = (phi(-hh / 2 * e3) - phi(hh / 2 * e3)) / hh     # <phi, [delta(x+h/2 e3)-delta(x-h/2 e3)]/h>
report("d_3 delta = lim [delta(x+h/2 e3)-delta(x-h/2 e3)]/h (sign check)", np.isclose(pairing, -d3phi0, rtol=1e-6),
       f"pairing {pairing:.6f}, -d_3 phi(0) {-d3phi0}")

# ----------------------------------------------------------------------------------------
print("\nExample 3: isotropic operator, curl of the equation of motion")
# Smooth analytic field and finite differences on a grid
Lb, Ng = 1.0, 41
g = np.linspace(-Lb, Lb, Ng); dg = g[1] - g[0]
G = np.meshgrid(g, g, g, indexing="ij")
U = np.stack([np.sin(G[0]) * np.cos(2 * G[1]) * G[2] ** 2,
              np.cos(G[0] * G[2]) + G[1] ** 3,
              np.exp(-G[0] ** 2 - G[1] ** 2) * np.sin(3 * G[2])])


def grad(fld):  # returns d_j fld
    return np.stack(np.gradient(fld, dg, edge_order=2))


def div(V):
    return sum(np.gradient(V[i], dg, axis=i, edge_order=2) for i in range(3))


def lap(fld):
    return sum(np.gradient(np.gradient(fld, dg, axis=i, edge_order=2), dg, axis=i, edge_order=2) for i in range(3))


def curl(V):
    dV = np.stack([grad(V[i]) for i in range(3)])        # dV[i, j] = d_j V_i
    return np.einsum("ijk,kj...->i...", eps, dV)


lam, mu = 2.3, 1.7
A = isotropic_A(lam, mu)
dU = np.stack([grad(U[k]) for k in range(3)])            # dU[k, l] = d_l u_k
stress = np.einsum("ijkl,kl...->ij...", A, dU)
op_full = np.stack([div(stress[i]) for i in range(3)])   # d_j (A_ijkl d_l u_k)
op_navier = (lam + mu) * grad(div(U)) + mu * np.stack([lap(U[i]) for i in range(3)])
inner = (slice(4, -4),) * 3
err = np.abs(op_full[(slice(None),) + inner] - op_navier[(slice(None),) + inner]).max()
scale = np.abs(op_navier[(slice(None),) + inner]).max()
report("d_j(A_ijkl d_l u_k) = (lam+mu) grad div u + mu lap u (isotropic, homogeneous)", err / scale < 1e-6,
       f"rel err {err/scale:.2e}")
# curl of the operator equals mu lap curl u  (i.e. curl kills the grad-div term)
lhs_c = curl(op_navier)
omega = curl(U)
rhs_c = mu * np.stack([lap(omega[i]) for i in range(3)])
err = np.abs(lhs_c[(slice(None),) + inner] - rhs_c[(slice(None),) + inner]).max()
scale = np.abs(rhs_c[(slice(None),) + inner]).max()
report("curl[(lam+mu) grad div u + mu lap u] = mu lap (curl u)", err / scale < 2e-2,
       f"rel err {err/scale:.2e} (finite differences of fourth order)")
# curl of a gradient vanishes (body force is a gradient)
cg = curl(grad(np.exp(-(G[0] ** 2 + G[1] ** 2 + G[2] ** 2) / 0.3)))
report("curl grad delta_smooth = 0", np.abs(cg[(slice(None),) + inner]).max() < 1e-10,
       f"max {np.abs(cg[(slice(None),) + inner]).max():.2e}")
# P-wave speed relation and explosion moment tensor invariance
rho = 3.0
alpha = np.sqrt((lam + 2 * mu) / rho)
Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
report("isotropic moment tensor invariant under rotation", np.allclose(Q @ (M0 * np.eye(3)) @ Q.T, M0 * np.eye(3)))
report("alpha^2 = (lam + 2 mu)/rho", np.isclose(alpha**2 * rho, lam + 2 * mu))

# ----------------------------------------------------------------------------------------
print("\nExample 4: finite rupture, moment rate")
mu4, L, Wf, D, vr = 3e10, 100e3, 20e3, 2.0, 3e3
Tr = L / vr
M0_tot = mu4 * L * Wf * D
report("rupture duration L/v_r = 33.3 s (quoted 33 s)", np.isclose(Tr, 33.3333, rtol=1e-3), f"{Tr:.3f} s")
report("total moment 1.2e20 N m", np.isclose(M0_tot, 1.2e20), f"{M0_tot:.3e}")
report("moment rate mu W_f D v_r = 3.6e18 N m/s", np.isclose(mu4 * Wf * D * vr, 3.6e18), f"{mu4*Wf*D*vr:.3e}")
Mw = 2 / 3 * (np.log10(M0_tot) - 9.1)
report("moment magnitude 7.3", np.isclose(Mw, 7.3, atol=0.05), f"{Mw:.3f}")


# M0(t) by quadrature of the slip over the fault, using a smoothed Heaviside
def M0_of_t(t, tau=0.2):
    xs = np.linspace(0, L, 200001)
    Hs = 0.5 * (1 + np.tanh((t - xs / vr) / tau))
    return mu4 * Wf * D * np.trapz(Hs, xs)


ts = np.array([-5.0, 10.0, 20.0, 30.0, 40.0, 60.0])
exact = np.where(ts < 0, 0.0, np.where(ts <= Tr, mu4 * Wf * D * vr * ts, M0_tot))
num = np.array([M0_of_t(t) for t in ts])
report("M0(t) piecewise-linear ramp (quadrature)", np.allclose(num, exact, rtol=1e-4, atol=1e-6 * M0_tot),
       f"max rel diff {np.abs(num-exact).max()/M0_tot:.2e}")
# moment rate by numerical differentiation is the boxcar of height mu Wf D vr inside (0, Tr)
dt = 1e-3
for t in [5.0, 15.0, 25.0]:
    rate = (M0_of_t(t + dt) - M0_of_t(t - dt)) / (2 * dt)
    report(f"dM0/dt at t={t:.0f} s equals boxcar height", np.isclose(rate, mu4 * Wf * D * vr, rtol=1e-4), f"{rate:.4e}")
rate_out = (M0_of_t(50 + dt) - M0_of_t(50 - dt)) / (2 * dt)
report("dM0/dt at t=50 s is zero", abs(rate_out) < 1e-6 * mu4 * Wf * D * vr, f"{rate_out:.2e}")
# area of the boxcar is the total moment
report("boxcar area = total moment", np.isclose(mu4 * Wf * D * vr * Tr, M0_tot))

# ----------------------------------------------------------------------------------------
print("\nExample 5: polynomial interpolation as a linear system")
xs3 = np.sort(rng.uniform(-2, 2, 3))
V3 = np.vander(xs3, 3, increasing=True)
det_formula = (xs3[1] - xs3[0]) * (xs3[2] - xs3[0]) * (xs3[2] - xs3[1])
report("(a) det V = (x2-x1)(x3-x1)(x3-x2)", np.isclose(np.linalg.det(V3), det_formula),
       f"{np.linalg.det(V3):.6f} vs {det_formula:.6f}")
y3 = rng.normal(size=3)
c3 = np.linalg.solve(V3, y3)
report("(a) solution interpolates the data", np.allclose(np.polyval(c3[::-1], xs3), y3))
# Lagrange form gives the same polynomial
ell = lambda i, x: np.prod([(x - xs3[j]) / (xs3[i] - xs3[j]) for j in range(3) if j != i])
xt = rng.uniform(-2, 2, 5)
p_lag = np.array([sum(y3[i] * ell(i, xx) for i in range(3)) for xx in xt])
report("(a) Lagrange form agrees with monomial solution", np.allclose(p_lag, np.polyval(c3[::-1], xt)))
# relative perturbation bound with the 2-norm condition number
cond3 = np.linalg.cond(V3, 2)
for _ in range(200):
    dy = rng.normal(size=3) * 1e-6
    dc = np.linalg.solve(V3, dy)
    assert np.linalg.norm(dc) / np.linalg.norm(c3) <= cond3 * np.linalg.norm(dy) / np.linalg.norm(y3) * (1 + 1e-9)
report("(a) relative error bound with cond(V) holds for random perturbations", True)
report("(a) cond(V) = sigma_max/sigma_min", np.isclose(cond3, np.linalg.svd(V3, compute_uv=False)[0] / np.linalg.svd(V3, compute_uv=False)[-1]))

# (b) n=3, m=1: left null vector and compatibility
V31 = np.vander(xs3, 2, increasing=True)
a = np.array([xs3[1] - xs3[2], xs3[2] - xs3[0], xs3[0] - xs3[1]])
report("(b) a^T V = 0", np.allclose(a @ V31, 0), f"{a @ V31}")
y_gen = rng.normal(size=3)
res_gen = np.linalg.lstsq(V31, y_gen, rcond=None)[1]
report("(b) generic data: no exact solution (nonzero residual)", res_gen.size > 0 and res_gen[0] > 1e-6)
c_line = rng.normal(size=2)
y_col = V31 @ c_line
report("(b) collinear data satisfy a . y = 0", abs(a @ y_col) < 1e-12, f"{a @ y_col:.2e}")
res_col = np.linalg.norm(V31 @ np.linalg.lstsq(V31, y_col, rcond=None)[0] - y_col)
report("(b) collinear data: exact solution exists", res_col < 1e-12)
# equivalence of a.y=0 with the slope form and the second-difference form
slope_form = (y_gen[2] - y_gen[0]) * (xs3[1] - xs3[0]) - (y_gen[1] - y_gen[0]) * (xs3[2] - xs3[0])
report("(b) a . y = -(slope form)", np.isclose(a @ y_gen, -slope_form))
xe = np.array([0.3, 0.3 + 0.7, 0.3 + 1.4]); ae = np.array([xe[1] - xe[2], xe[2] - xe[0], xe[0] - xe[1]])
report("(b) equally spaced: a proportional to (1,-2,1)", np.allclose(ae / ae[0], [1, -2, 1]))
# explicit solution when compatible
c1 = (y_col[1] - y_col[0]) / (xs3[1] - xs3[0]); c0 = y_col[0] - c1 * xs3[0]
report("(b) c1, c0 formulas reproduce the line", np.allclose([c0, c1], c_line))

# (c) n=2, m=2: null vector and family
xs2 = xs3[:2]; y2 = rng.normal(size=2)
V22 = np.vander(xs2, 3, increasing=True)
c_null = np.array([xs2[0] * xs2[1], -(xs2[0] + xs2[1]), 1.0])
report("(c) V c0 = 0 for c0 = (x1 x2, -(x1+x2), 1)", np.allclose(V22 @ c_null, 0))
report("(c) null space is one-dimensional", np.linalg.matrix_rank(V22) == 2)
sl = (y2[1] - y2[0]) / (xs2[1] - xs2[0])
for alpha in rng.normal(size=5):
    c_fam = np.array([y2[0] - sl * xs2[0] + alpha * xs2[0] * xs2[1], sl - alpha * (xs2[0] + xs2[1]), alpha])
    assert np.allclose(V22 @ c_fam, y2)
report("(c) one-parameter family c(alpha) fits the data for all alpha", True)


# (d) condition numbers: floating point and exact rational inverse
def exact_inverse_vandermonde(n):
    xs = [Fraction(i, n - 1) for i in range(n)]
    Aug = [[xi ** k for k in range(n)] + [Fraction(int(i == j)) for j in range(n)] for i, xi in enumerate(xs)]
    for c in range(n):
        p = next(r for r in range(c, n) if Aug[r][c] != 0)
        Aug[c], Aug[p] = Aug[p], Aug[c]
        piv = Aug[c][c]
        Aug[c] = [v / piv for v in Aug[c]]
        for r in range(n):
            if r != c and Aug[r][c] != 0:
                fac = Aug[r][c]
                Aug[r] = [u - fac * v for u, v in zip(Aug[r], Aug[c])]
    return np.array([[float(v) for v in row[n:]] for row in Aug])


quoted = {5: 6.9e2, 10: 1.5e7, 15: 4.0e11, 20: 1.1e16}
for n in [5, 10, 15, 20]:
    xn = np.linspace(0, 1, n)
    Vn = np.vander(xn, increasing=True)
    cond_fp = np.linalg.cond(Vn, 2)
    Vinv = exact_inverse_vandermonde(n)
    cond_exact = np.linalg.norm(Vn, 2) * np.linalg.norm(Vinv, 2)
    ok = abs(cond_exact - quoted[n]) / quoted[n] < 0.05            # agreement to 2 significant figures
    report(f"(d) cond(V), n={n}: quoted {quoted[n]:.1e}", ok,
           f"exact-inverse {cond_exact:.4e}, numpy {cond_fp:.4e}")
report("(d) numpy value for n=20 is 1.2e16 as stated in footnote",
       abs(np.linalg.cond(np.vander(np.linspace(0, 1, 20), increasing=True), 2) - 1.2e16) / 1.2e16 < 0.05)
# amplification demonstration: n=10, data perturbed at 1e-7 gives O(1) relative coefficient error
xn = np.linspace(0, 1, 10); Vn = np.vander(xn, increasing=True)
u, sv, vt = np.linalg.svd(Vn)
cn = rng.normal(size=10); yn = Vn @ cn
dy = u[:, -1] * 1e-7 * np.linalg.norm(yn)           # worst-case direction
dc = np.linalg.solve(Vn, dy)
rel = np.linalg.norm(dc) / np.linalg.norm(cn)
report("(d) n=10: 1e-7 relative data error can give O(1) coefficient error", rel > 0.1, f"relative error {rel:.2f}")

print("\nAll checks passed." if failures == 0 else f"\n{failures} CHECK(S) FAILED.")
