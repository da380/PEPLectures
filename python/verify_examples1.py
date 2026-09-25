"""Numerical checks for the worked examples for Lecture 12 (examples/examples1.tex).

Every numerical or closed-form claim made in the solutions is checked here.
Prints OK/FAIL lines; exits with a non-zero status if anything fails.
"""
import numpy as np
from scipy.linalg import sqrtm

rng = np.random.default_rng(1)
failures = 0


def check(name, ok, val=""):
    global failures
    print(f"{'OK  ' if ok else 'FAIL'} {name}  {val}")
    if not ok:
        failures += 1


def rot3(axis, ang):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


print("\n=== Example 1: simple shear ===")
gamma, t = 0.8, 1.25
a = gamma * t
F = np.array([[1, a, 0], [0, 1, 0], [0, 0, 1.0]])
check("J = 1", abs(np.linalg.det(F) - 1) < 1e-14, np.linalg.det(F))
# inverse motion
x = rng.normal(size=3)
y = np.array([x[0] + a * x[1], x[1], x[2]])
xinv = np.array([y[0] - a * y[1], y[1], y[2]])
check("inverse motion phi^{-1}(phi(x)) = x", np.allclose(xinv, x))
# spatial velocity field v(phi^{-1}(y,t),t) = (gamma y2, 0, 0)
v_ref = np.array([gamma * xinv[1], 0, 0])
check("spatial velocity = (gamma y2,0,0)", np.allclose(v_ref, [gamma * y[1], 0, 0]))
# spatial density for rho = rho0 (1 + alpha x1)
rho0, alpha = 3.3, 0.4
varrho = rho0 * (1 + alpha * (y[0] - a * y[1]))
check("spatial density rho0(1+alpha(y1 - gamma t y2))", abs(varrho - rho0 * (1 + alpha * x[0])) < 1e-12)
C = F.T @ F
lam = np.sort(np.linalg.eigvalsh(C))
lam_p = 1 + a**2 / 2 + a * np.sqrt(1 + a**2 / 4)
lam_m = 1 + a**2 / 2 - a * np.sqrt(1 + a**2 / 4)
check("eigenvalues of C = {lam-, 1, lam+}", np.allclose(lam, [lam_m, 1, lam_p]), lam)
check("lam+ lam- = 1", abs(lam_p * lam_m - 1) < 1e-12)
sp, sm = np.sqrt(1 + a**2 / 4) + a / 2, np.sqrt(1 + a**2 / 4) - a / 2
check("principal stretches sqrt(1+a^2/4) +/- a/2", np.allclose([sp**2, sm**2], [lam_p, lam_m]))
up = np.array([1, sp, 0])
um = np.array([1, -sm, 0])
check("C u+ = lam+ u+", np.allclose(C @ up, lam_p * up))
check("C u- = lam- u-", np.allclose(C @ um, lam_m * um))
check("u+ . u- = 0", abs(up @ um) < 1e-12)
# small-a limit: eigenvectors -> (1, +/-1)
a_small = 1e-4
check("small a: lam+- = 1 +- a + O(a^2)",
      abs((1 + a_small**2 / 2 + a_small * np.sqrt(1 + a_small**2 / 4)) - (1 + a_small)) < 1e-8)

print("\n=== Example 2: polar decomposition ===")
F = np.array([[2.0, 1.0], [0.0, 1.0]])
C = F.T @ F
check("C = [[4,2],[2,2]]", np.allclose(C, [[4, 2], [2, 2]]))
lam, Q = np.linalg.eigh(C)
s5 = np.sqrt(5.0)
check("eigenvalues 3 -/+ sqrt5", np.allclose(lam, [3 - s5, 3 + s5]), lam)
up, um = np.array([2, -1 + s5]), np.array([2, -1 - s5])
check("C u+ = (3+sqrt5) u+", np.allclose(C @ up, (3 + s5) * up))
check("C u- = (3-sqrt5) u-", np.allclose(C @ um, (3 - s5) * um))
check("u+ . u- = 0", abs(up @ um) < 1e-12)
check("|u+|^2 = 10 - 2 sqrt5", abs(up @ up - (10 - 2 * s5)) < 1e-12)
check("|u-|^2 = 10 + 2 sqrt5", abs(um @ um - (10 + 2 * s5)) < 1e-12)
alpha = np.degrees(np.arctan2(up[1], up[0]))
check("angle of q+ = arctan((sqrt5-1)/2) ~ 31.7 deg", abs(alpha - 31.717474) < 1e-4, alpha)
check("sqrt(3+sqrt5) = (sqrt5+1)/sqrt2", abs(np.sqrt(3 + s5) - (s5 + 1) / np.sqrt(2)) < 1e-14)
check("sqrt(3-sqrt5) = (sqrt5-1)/sqrt2", abs(np.sqrt(3 - s5) - (s5 - 1) / np.sqrt(2)) < 1e-14)
# projectors
Pp = np.outer(up, up) / (up @ up)
check("q+ q+^T = [[(5+sqrt5)/10, sqrt5/5],[sqrt5/5,(5-sqrt5)/10]]",
      np.allclose(Pp, [[(5 + s5) / 10, s5 / 5], [s5 / 5, (5 - s5) / 10]]))
U = Q @ np.diag(np.sqrt(lam)) @ Q.T
U_cf = np.array([[6, 2], [2, 4]]) / np.sqrt(10)
check("U = (1/sqrt10)[[6,2],[2,4]] (eigendecomposition)", np.allclose(U, U_cf))
check("U = sqrt(2/5)[[3,1],[1,2]]", np.allclose(U_cf, np.sqrt(0.4) * np.array([[3, 1], [1, 2]])))
check("U agrees with scipy sqrtm", np.allclose(U, np.real(sqrtm(C))))
check("U^2 = C", np.allclose(U @ U, C))
check("tr U = sqrt10, det U = 2", abs(np.trace(U) - np.sqrt(10)) < 1e-12 and abs(np.linalg.det(U) - 2) < 1e-12)
# Cayley-Hamilton shortcut
U_ch = (C + 2 * np.eye(2)) / np.sqrt(10)
check("Cayley-Hamilton: U = (C + 2 I)/sqrt10", np.allclose(U_ch, U))
Uinv = np.linalg.inv(U)
check("U^{-1} = (1/sqrt10)[[2,-1],[-1,3]]", np.allclose(Uinv, np.array([[2, -1], [-1, 3]]) / np.sqrt(10)))
R = F @ Uinv
R_cf = np.array([[3, 1], [-1, 3]]) / np.sqrt(10)
check("R = (1/sqrt10)[[3,1],[-1,3]]", np.allclose(R, R_cf))
check("R^T R = 1", np.allclose(R.T @ R, np.eye(2)))
check("det R = 1", abs(np.linalg.det(R) - 1) < 1e-12)
theta = np.arctan2(R[1, 0], R[0, 0])
check("tan theta = -1/3, theta ~ -18.4 deg", abs(np.tan(theta) + 1 / 3) < 1e-12, np.degrees(theta))
V = R @ U @ R.T
V_cf = np.array([[7, 1], [1, 3]]) / np.sqrt(10)
check("V = R U R^T = (1/sqrt10)[[7,1],[1,3]]", np.allclose(V, V_cf))
check("V = F R^T", np.allclose(F @ R.T, V))
check("V^2 = F F^T = [[5,1],[1,1]]", np.allclose(V @ V, F @ F.T) and np.allclose(F @ F.T, [[5, 1], [1, 1]]))
check("F = R U", np.allclose(R @ U, F))
check("F = V R", np.allclose(V @ R, F))
wp = np.array([1, s5 - 2])
wm = np.array([1, -s5 - 2])
check("V w+ = sqrt(lam+) w+", np.allclose(V @ wp, np.sqrt(3 + s5) * wp))
check("V w- = sqrt(lam-) w-", np.allclose(V @ wm, np.sqrt(3 - s5) * wm))
Rup = R @ up
check("R q+ parallel to w+", abs(Rup[0] * wp[1] - Rup[1] * wp[0]) < 1e-12)
beta = np.degrees(np.arctan2(wp[1], wp[0]))
check("angle of w+ ~ 13.3 deg = 31.7 - 18.4", abs(beta - (alpha + np.degrees(theta))) < 1e-9, beta)
semi = np.sqrt(lam[::-1])
check("semi-axes of ellipse 2.288, 0.874; product = det F = 2",
      np.allclose(semi, [2.2882456, 0.8740320], atol=1e-6) and abs(semi[0] * semi[1] - 2) < 1e-12, semi)

print("\n=== Example 3: rigid motions ===")
Qr = rot3([0.3, -1.0, 0.7], 1.1)
one = np.eye(3)
check("F = Q: J = 1", abs(np.linalg.det(Qr) - 1) < 1e-12)
check("C = Q^T Q = 1", np.allclose(Qr.T @ Qr, one))
mu = 1.7
W1 = lambda F: 0.5 * mu * np.trace(F.T @ F - one)
W2 = lambda F: 0.5 * mu * np.trace(F - one)
W3 = lambda F: 0.5 * mu * np.sum((F - one) ** 2)
Fr = rng.normal(size=(3, 3)) + 2 * one
check("W1(QF) = W1(F)", abs(W1(Qr @ Fr) - W1(Fr)) < 1e-12)
Qz = rot3([0, 0, 1], np.pi / 2)
check("Q (pi/2 about x3) has tr Q = 1", abs(np.trace(Qz) - 1) < 1e-12)
check("W2(Q) = -mu != W2(1) = 0", abs(W2(Qz) + mu) < 1e-12 and W2(one) == 0, W2(Qz))
check("W3(Q) = 2 mu != W3(1) = 0", abs(W3(Qz) - 2 * mu) < 1e-12 and W3(one) == 0, W3(Qz))
psi = 0.9
Qpsi = rot3([1, 2, 0.5], psi)
check("W3(Q) = 2 mu (1 - cos psi)", abs(W3(Qpsi) - 2 * mu * (1 - np.cos(psi))) < 1e-12)
# kinetic energy of a rigid motion of the unit cube centred at the origin, uniform rho
rho = 2.0
n = 40
g = (np.arange(n) + 0.5) / n - 0.5
X = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
dV = 1.0 / n**3
def Qt(tt):
    return rot3([1, 2, 0.5], 0.7 * tt) @ rot3([0, 1, -1], 0.3 * tt**2)
adot = np.array([0.2, -0.5, 0.1])
tt, h = 0.6, 1e-5
Qd = (Qt(tt + h) - Qt(tt - h)) / (2 * h)
Q0 = Qt(tt)
vel = adot + X @ Qd.T
T_direct = 0.5 * rho * np.sum(vel**2) * dV
Om = Q0.T @ Qd
check("Omega = Q^T Qdot antisymmetric", np.allclose(Om, -Om.T, atol=1e-6))
omega = np.array([Om[2, 1], Om[0, 2], Om[1, 0]])
M = rho
I0 = rho * dV * (np.sum(np.sum(X**2, 1)) * one - X.T @ X)
check("inertia tensor of unit cube = M/6 1", np.allclose(I0, M / 6 * one, rtol=1e-3))
T_formula = 0.5 * M * adot @ adot + 0.5 * omega @ I0 @ omega
check("T = M|adot|^2/2 + omega.I omega/2", abs(T_direct - T_formula) < 1e-9, (T_direct, T_formula))
check("omega x x = Omega x", np.allclose(np.cross(omega, X[5]), Om @ X[5], atol=1e-6))

print("\n=== Example 4: uniform stretch of the unit cube ===")
lam = np.array([1.2, 0.8, 1.5])
lamdot = np.array([0.3, -0.2, 0.1])
rho = 2.7
J = np.prod(lam)
check("volume of M_t = J", abs(J - 1.44) < 1e-12, J)
T_ref = rho / 6 * np.sum(lamdot**2)
n = 60
g = (np.arange(n) + 0.5) / n
X = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
T_ref_num = 0.5 * rho * np.sum((X * lamdot) ** 2) / n**3
check("referential kinetic energy = rho/6 sum lamdot^2", abs(T_ref_num - T_ref) < 1e-4 * T_ref, (T_ref_num, T_ref))
Y = X * lam                           # midpoint grid on M_t, cell volume J/n^3
varrho = rho / J
vspat = (lamdot / lam) * Y
T_sp_num = 0.5 * varrho * np.sum(vspat**2) * J / n**3
check("spatial kinetic energy agrees", abs(T_sp_num - T_ref) < 1e-4 * T_ref, (T_sp_num, T_ref))
check("mass from spatial density", abs(varrho * J - rho) < 1e-12)
check("int_{M_t} y1^2 d3y = J lam1^2/3", abs(np.sum(Y[:, 0] ** 2) * J / n**3 - J * lam[0] ** 2 / 3) < 1e-4)

print("\n=== Example 5: one-dimensional rod ===")
E, rho, L = 4.0, 1.6, 2.0
c = np.sqrt(E / rho)
W = lambda Fv: 0.5 * E * (Fv - 1) ** 2
Fv, h = 1.3, 1e-6
check("T = dW/dF = E(F-1)", abs((W(Fv + h) - W(Fv - h)) / (2 * h) - E * (Fv - 1)) < 1e-6)
# normal mode phi = x + A cos(k x) cos(omega t) with k = n pi / L satisfies PDE and BC
nmode, A = 3, 0.05
k = nmode * np.pi / L
om = c * k
xs = np.linspace(0, L, 201)
t0 = 0.37
phi = lambda xx, tt: xx + A * np.cos(k * xx) * np.cos(om * tt)
hh = 1e-4
phi_tt = (phi(xs, t0 + hh) - 2 * phi(xs, t0) + phi(xs, t0 - hh)) / hh**2
phi_xx = (phi(xs + hh, t0) - 2 * phi(xs, t0) + phi(xs - hh, t0)) / hh**2
check("rho phi_tt = E phi_xx for the normal mode", np.max(np.abs(rho * phi_tt - E * phi_xx)) < 1e-4)
phi_x = (phi(xs + hh, t0) - phi(xs - hh, t0)) / (2 * hh)
check("phi_x = 1 (T = 0) at x = 0, L", abs(phi_x[0] - 1) < 1e-8 and abs(phi_x[-1] - 1) < 1e-8)
check("omega_n = n pi c / L", abs(om - nmode * np.pi * np.sqrt(E / rho) / L) < 1e-14)

print("\n=== Example 6: traction ===")
# Jacobi's formula by finite differences on a random F
F = rng.normal(size=(3, 3)) + 2 * np.eye(3)
J = np.linalg.det(F)
Finv = np.linalg.inv(F)
dJ = np.zeros((3, 3))
h = 1e-6
for i in range(3):
    for j in range(3):
        Fp, Fm = F.copy(), F.copy()
        Fp[i, j] += h
        Fm[i, j] -= h
        dJ[i, j] = (np.linalg.det(Fp) - np.linalg.det(Fm)) / (2 * h)
check("dJ/dF_ij = J F^{-1}_ji (random F)", np.allclose(dJ, J * Finv.T, atol=1e-6))
check("det(1 + eps A) = 1 + eps tr A + O(eps^2)",
      abs(np.linalg.det(np.eye(3) + 1e-4 * F) - (1 + 1e-4 * np.trace(F))) < 1e-6)
kappa = 5.0
lam = np.array([1.2, 0.8, 1.5])
F = np.diag(lam)
J = np.prod(lam)
T = kappa * (J - 1) * J * np.linalg.inv(F).T
check("T = kappa (J-1) J diag(1/lam)", np.allclose(T, kappa * (J - 1) * J * np.diag(1 / lam)))
# finite-difference check of T = dW/dF for this diagonal F
Wf = lambda FF: 0.5 * kappa * (np.linalg.det(FF) - 1) ** 2
Tfd = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        Fp, Fm = F.copy(), F.copy()
        Fp[i, j] += h
        Fm[i, j] -= h
        Tfd[i, j] = (Wf(Fp) - Wf(Fm)) / (2 * h)
check("T = dW/dF by finite differences", np.allclose(Tfd, T, atol=1e-6))
n1 = np.array([1.0, 0, 0])
t_ref = T @ n1
check("reference traction on x1=1 face = kappa (J-1) lam2 lam3 e1",
      np.allclose(t_ref, [kappa * (J - 1) * lam[1] * lam[2], 0, 0]), t_ref)
force = t_ref * 1.0                    # reference area of the face is 1
area_now = lam[1] * lam[2]
t_cauchy = force / area_now
check("Cauchy traction = kappa (J-1)", np.allclose(t_cauchy, [kappa * (J - 1), 0, 0]), t_cauchy)
check("ratio = lam2 lam3", abs(t_ref[0] / t_cauchy[0] - lam[1] * lam[2]) < 1e-12)
sigma = T @ F.T / J
check("sigma = J^{-1} T F^T = kappa (J-1) 1", np.allclose(sigma, kappa * (J - 1) * np.eye(3)))
check("W(QF) = W(F) for W = kappa (J-1)^2 / 2", abs(Wf(Qr @ F) - Wf(F)) < 1e-12)

print()
if failures:
    print(f"{failures} check(s) FAILED")
    raise SystemExit(1)
print("All checks passed.")
