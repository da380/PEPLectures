"""Numerical checks for examples/examples9.tex (Lecture 20: Self-gravitation and rotation).

Every numerical or symbolic-by-numbers claim in the solutions is checked here and reported
as an OK/FAIL line.  Run with  python3 verify_examples9.py ; exit status is non-zero on any FAIL.
"""
import numpy as np
from scipy.integrate import quad

G = 6.674e-11
M = 5.97e24
b = 6.371e6
Om = 7.29e-5

fails = 0


def check(name, value, expected, rtol=1e-6, atol=0.0):
    global fails
    ok = np.allclose(value, expected, rtol=rtol, atol=atol)
    fails += 0 if ok else 1
    print(f"{'OK  ' if ok else 'FAIL'} {name}: got {value!r}, expected {expected!r}")


eps = np.zeros((3, 3, 3))
eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1
eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1

# ---------------------------------------------------------------------------
print("\n=== Example 1: homogeneous sphere ===")
rho = 3 * M / (4 * np.pi * b**3)


def phi(r):
    return np.where(r <= b, -G * M / (2 * b**3) * (3 * b**2 - r**2), -G * M / r)


def gr(r):
    return np.where(r <= b, -G * M * r / b**3, -G * M / r**2)


# (i) phi solves Poisson's equation inside and Laplace's outside (finite differences).
for r0, rhs in [(0.3 * b, 4 * np.pi * G * rho), (0.8 * b, 4 * np.pi * G * rho),
                (1.5 * b, 0.0), (3.0 * b, 0.0)]:
    h = 1e-3 * b
    r = np.array([r0 - h, r0, r0 + h])
    p = phi(r)
    lap = (p[2] - 2 * p[1] + p[0]) / h**2 + (2 / r0) * (p[2] - p[0]) / (2 * h)
    check(f"Poisson at r/b={r0/b:.1f}: (1/r^2)(r^2 phi')'", lap, rhs,
          rtol=1e-5, atol=1e-6 * 4 * np.pi * G * rho)

# (ii) g_r = -dphi/dr.
for r0 in [0.5 * b, 2.0 * b]:
    h = 1e-4 * b
    check(f"g_r = -dphi/dr at r/b={r0/b:.1f}", -(phi(r0 + h) - phi(r0 - h)) / (2 * h), gr(r0), rtol=1e-6)

# (iii) continuity at r=b and the stated central value.
check("phi continuous at b", phi(b * (1 - 1e-12)), phi(b * (1 + 1e-12)))
check("g_r continuous at b", gr(b * (1 - 1e-12)), gr(b * (1 + 1e-12)))
check("phi(0) = -3GM/(2b)", float(phi(np.array(0.0))), -1.5 * G * M / b)
check("phi(0)/phi(b) = 3/2", float(phi(np.array(0.0)) / phi(np.array(b))), 1.5)
check("g_r inside = -(4/3) pi G rho r", float(gr(np.array(0.4 * b))), -4 / 3 * np.pi * G * rho * 0.4 * b)

# (iv) binding energy by the volume form 1/2 int rho phi.
Vg_vol = 0.5 * rho * quad(lambda r: phi(np.array(r)) * 4 * np.pi * r**2, 0, b)[0]
Vg_exact = -0.6 * G * M**2 / b
check("V_g (1/2 int rho phi) = -3GM^2/(5b)", Vg_vol, Vg_exact)

# (v) binding energy by the field form -(1/8 pi G) int |grad phi|^2, interior/exterior split.
I_in = quad(lambda r: gr(np.array(r))**2 * 4 * np.pi * r**2, 0, b)[0]
I_out = quad(lambda r: gr(np.array(r))**2 * 4 * np.pi * r**2, b, np.inf)[0]
check("int |grad phi|^2 = (24 pi/5) G^2 M^2 / b", I_in + I_out, 24 * np.pi / 5 * G**2 * M**2 / b)
check("V_g (field form) = -3GM^2/(5b)", -(I_in + I_out) / (8 * np.pi * G), Vg_exact)
check("exterior share of field energy = 5/6", I_out / (I_in + I_out), 5 / 6)

# (vi) numbers quoted.
check("GM^2/b = 3.73e32 J", G * M**2 / b, 3.73e32, rtol=2e-3)
check("V_g(Earth) = -2.24e32 J", Vg_exact, -2.24e32, rtol=2e-3)
check("mean density = 5510 kg/m^3", rho, 5510, rtol=1e-3)
check("surface gravity = 9.82 m/s^2", G * M / b**2, 9.82, rtol=1e-3)
ratio = -Vg_exact / np.array([1e19, 1e18])
check("|V_g| / (1e18..1e19 J) within 1e13..1e14", (ratio[0] >= 1e13) and (ratio[1] <= 1e15), True)
print(f"     |V_g|/1e19 J = {ratio[0]:.2e}, |V_g|/1e18 J = {ratio[1]:.2e}")

# ---------------------------------------------------------------------------
print("\n=== Example 2: rigid motion of a point-mass cloud ===")
rng = np.random.default_rng(1)
N = 200
X = rng.normal(size=(N, 3))            # reference positions
m = rng.uniform(0.5, 1.5, size=N)      # masses (rho d^3x)


def deform(X):
    # some smooth non-rigid motion phi(x)
    return X + 0.2 * np.sin(X[:, [1, 2, 0]]) + 0.1 * X**2


def zeta_gamma(Y, i0):
    d = Y[i0] - np.delete(Y, i0, axis=0)
    mm = np.delete(m, i0)
    r = np.linalg.norm(d, axis=1)
    return -G * np.sum(mm / r), -G * np.sum((mm / r**3)[:, None] * d, axis=0)


def random_rotation(rng):
    A = rng.normal(size=(3, 3))
    Q, _ = np.linalg.qr(A)
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


Q = random_rotation(rng)
a = rng.normal(size=3) * 5
Y = deform(X)
Yt = a + Y @ Q.T
for i0 in [0, 17, 123]:
    z, g = zeta_gamma(Y, i0)
    zt, gt = zeta_gamma(Yt, i0)
    check(f"zeta unchanged at particle {i0}", zt, z, rtol=1e-12)
    check(f"gamma -> Q gamma at particle {i0}", gt, Q @ g, rtol=1e-10)
check("Q orthogonal", Q.T @ Q, np.eye(3), atol=1e-12)

# ---------------------------------------------------------------------------
print("\n=== Example 3: steadily rotating frame ===")


def Rz(t):
    c, s = np.cos(Om * t), np.sin(Om * t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


t0, h = 1234.5, 1e-3
Rdot = (Rz(t0 + h) - Rz(t0 - h)) / (2 * h)
A = Rz(t0).T @ Rdot                       # A_{jl} = R_{kj} dR_{kl}/dt
check("R^T dR/dt (rotation about e3)", A, np.array([[0, -Om, 0], [Om, 0, 0], [0, 0, 0]]), atol=1e-12)
check("R^T dR/dt anti-symmetric", A + A.T, np.zeros((3, 3)), atol=1e-12)
Omega_vec = 0.5 * np.einsum('jnl,jl->n', eps, A)      # eq. (15)
check("Omega read off = Omega e3", Omega_vec, np.array([0, 0, Om]), atol=1e-12)
check("epsilon_{jml} Omega_m reproduces R^T dR/dt", np.einsum('jml,m->jl', eps, Omega_vec), A, atol=1e-12)

# general axis: rotation about unit vector n by angle Om t via Rodrigues; check eq. (15) recovers Om n.
n = np.array([0.3, -0.5, 0.8]); n /= np.linalg.norm(n)
K = np.einsum('ijk,j->ik', eps, n)        # K_ik = eps_ijk n_j, so that K v = n x v


def Rn(t):
    th = Om * t
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K


Rdot_n = (Rn(t0 + h) - Rn(t0 - h)) / (2 * h)
A_n = Rn(t0).T @ Rdot_n
check("general axis: R^T dR/dt anti-symmetric", A_n + A_n.T, np.zeros((3, 3)), atol=1e-12)
check("general axis: Omega read off = Omega n", 0.5 * np.einsum('jnl,jl->n', eps, A_n), Om * n, atol=1e-12)

# Coriolis components and centrifugal components for Omega = Om e3.
v = rng.normal(size=3)
x = rng.normal(size=3)
Ov = np.array([0, 0, Om])
check("Coriolis 2 eps Omega v = 2 Om(-v2, v1, 0)", 2 * np.einsum('ijk,j,k->i', eps, Ov, v),
      2 * Om * np.array([-v[1], v[0], 0]))
cent = np.einsum('ijk,klm,j,l,m->i', eps, eps, Ov, Ov, x)
check("centrifugal = -Om^2 (x1, x2, 0)", cent, -Om**2 * np.array([x[0], x[1], 0]))
check("centrifugal = Omega(Omega.x) - |Omega|^2 x", cent, Ov * (Ov @ x) - (Ov @ Ov) * x)


def psi_special(x):
    return -0.5 * Om**2 * (x[0]**2 + x[1]**2)


def psi_general(x, Ov):
    Mmat = np.outer(Ov, Ov) - (Ov @ Ov) * np.eye(3)
    return 0.5 * x @ Mmat @ x


def grad(f, x, h=1e-6):
    g = np.zeros(3)
    for i in range(3):
        e = np.zeros(3); e[i] = h
        g[i] = (f(x + e) - f(x - e)) / (2 * h)
    return g


check("grad psi (special) = centrifugal acceleration", grad(psi_special, x), cent, atol=1e-12)
check("psi general = psi special for Omega = Om e3", psi_general(x, Ov), psi_special(x))
Ogen = Om * n
cent_gen = np.einsum('ijk,klm,j,l,m->i', eps, eps, Ogen, Ogen, x)
check("grad psi (general axis) = centrifugal acceleration", grad(lambda y: psi_general(y, Ogen), x),
      cent_gen, atol=1e-12)
# Laplacian of psi = -2 Omega^2
h = 1e-3
lap = sum((psi_special(x + h * e) - 2 * psi_special(x) + psi_special(x - h * e)) / h**2 for e in np.eye(3))
check("Laplacian psi = -2 Om^2", lap, -2 * Om**2, rtol=1e-6)
check("Om^2 b = 3.39e-2 m/s^2", Om**2 * b, 3.39e-2, rtol=2e-3)
check("Om^2 b^3/(GM) = 3.4e-3", Om**2 * b**3 / (G * M), 3.4e-3, rtol=2e-2)
check("1/(Om^2 b^3/GM) ~ 290", 1 / (Om**2 * b**3 / (G * M)), 290, rtol=1e-2)

# ---------------------------------------------------------------------------
print("\n=== Example 4: scaling table ===")
rb, vb = 5000.0, 8000.0
check("4 pi G rhobar = 4.19e-6 s^-2", 4 * np.pi * G * rb, 4.19e-6, rtol=2e-3)
check("L_g = 3900 km", vb / np.sqrt(4 * np.pi * G * rb), 3.9e6, rtol=3e-3)
check("1/Om = 1.37e4 s", 1 / Om, 1.37e4, rtol=2e-3)
check("1/Om = 3.8 h", 1 / Om / 3600, 3.8, rtol=5e-3)
table = {  # (T, L): (grav/elastic, Om T, 2 (Om T)^2)
    "(a) 1 Hz body wave": ((1.0, 1e4), (6.55e-6, 7.29e-5, 1.06e-8)),
    "(b) 100 s surface wave": ((100.0, 4e5), (1.05e-2, 7.29e-3, 1.06e-4)),
    "(c) 0S2": ((54 * 60.0, b), (2.66, 0.236, 0.112)),
}
for name, ((T, L), exp) in table.items():
    vals = (4 * np.pi * G * rb * L**2 / vb**2, Om * T, 2 * (Om * T)**2)
    check(f"{name}: grav/elastic", vals[0], exp[0], rtol=5e-3)
    check(f"{name}: Om T", vals[1], exp[1], rtol=5e-3)
    check(f"{name}: 2 (Om T)^2", vals[2], exp[2], rtol=5e-3)
check("(c) grav/elastic ~ 3 (lecture)", 4 * np.pi * G * rb * b**2 / vb**2, 3.0, rtol=0.15)
check("(c) grav/elastic = (L/L_g)^2", 4 * np.pi * G * rb * b**2 / vb**2,
      (b / (vb / np.sqrt(4 * np.pi * G * rb)))**2)

# ---------------------------------------------------------------------------
print("\n=== Example 5: gamma^1 for rigid displacements (and the lecture formula itself) ===")


def gamma_exact(Y, i0):
    return zeta_gamma(Y, i0)[1]


def gamma1_formula(X, U, i0):
    r = X[i0] - np.delete(X, i0, axis=0)
    du = U[i0] - np.delete(U, i0, axis=0)
    mm = np.delete(m, i0)
    rn = np.linalg.norm(r, axis=1)
    kern = np.eye(3)[None] / rn[:, None, None]**3 - 3 * r[:, :, None] * r[:, None, :] / rn[:, None, None]**5
    return -G * np.einsum('n,nij,nj->i', mm, kern, du)


# (0) validate the lecture's formula for gamma^1 against a finite difference of the exact gamma
#     for a smooth non-rigid displacement.
U = 0.3 * np.sin(X[:, [2, 0, 1]]) + 0.1 * X * X[:, [1, 2, 0]]
s = 1e-5
for i0 in [3, 60]:
    fd = (gamma_exact(X + s * U, i0) - gamma_exact(X - s * U, i0)) / (2 * s)
    check(f"lecture formula for gamma^1 vs finite difference, particle {i0}", gamma1_formula(X, U, i0), fd, rtol=1e-6)

# (a) translation
u0 = np.array([0.7, -1.3, 2.1])
U0 = np.tile(u0, (N, 1))
for i0 in [3, 60]:
    check(f"translation: gamma^1 = 0 at particle {i0}", gamma1_formula(X, U0, i0), np.zeros(3), atol=1e-30)
    fd = (gamma_exact(X + s * U0, i0) - gamma_exact(X - s * U0, i0)) / (2 * s)
    check(f"translation: finite-difference gamma^1 = 0 at particle {i0}", fd, np.zeros(3),
          atol=1e-8 * np.linalg.norm(gamma_exact(X, i0)))

# (b) rotation u = theta x x
theta = np.array([0.4, 0.9, -0.2])
Urot = np.cross(theta, X)
for i0 in [3, 60]:
    g0 = gamma_exact(X, i0)
    check(f"rotation: gamma^1 = theta x gamma^0 at particle {i0}", gamma1_formula(X, Urot, i0),
          np.cross(theta, g0), rtol=1e-10)
    fd = (gamma_exact(X + s * Urot, i0) - gamma_exact(X - s * Urot, i0)) / (2 * s)
    check(f"rotation: finite difference = theta x gamma^0 at particle {i0}", fd, np.cross(theta, g0), rtol=1e-6)

# infinitesimal rotation matrix Q = I + s eps_{ikj} theta_k, i.e. Q x = x + s theta x x
Qs = np.eye(3) + s * np.einsum('ikj,k->ij', eps, theta)
check("Q_ij = delta_ij + s eps_ikj theta_k gives Q x = x + s theta x x", Qs @ x, x + s * np.cross(theta, x), rtol=1e-12)

print("\n" + ("ALL CHECKS PASSED" if fails == 0 else f"{fails} CHECK(S) FAILED"))
raise SystemExit(1 if fails else 0)
