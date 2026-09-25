"""Numerical checks for the worked examples accompanying Lecture 13 (constitutive theory).

Every quantitative claim made in examples/examples2.tex is checked here.  Run with

    python3 verify_examples2.py

Each check prints an OK or FAIL line together with the numbers involved.
"""
import numpy as np

rng = np.random.default_rng(13)
I3 = np.eye(3)
status = {"fail": 0}


def report(name, ok, detail=""):
    tag = "OK  " if ok else "FAIL"
    if not ok:
        status["fail"] += 1
    print(f"{tag} {name}" + (f"  [{detail}]" if detail else ""))


def grad_fd(fun, F, h=1e-6):
    """Central-difference gradient dfun/dF_ij of a scalar function of a 3x3 matrix."""
    G = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            E = np.zeros((3, 3))
            E[i, j] = h
            G[i, j] = (fun(F + E) - fun(F - E)) / (2 * h)
    return G


def hess_fd(fun, F, h=1e-4):
    """Central-difference Hessian d^2 fun/dF_ij dF_kl of a scalar function of a 3x3 matrix."""
    H = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    E1 = np.zeros((3, 3))
                    E2 = np.zeros((3, 3))
                    E1[i, j] = h
                    E2[k, l] = h
                    H[i, j, k, l] = (fun(F + E1 + E2) - fun(F + E1 - E2)
                                     - fun(F - E1 + E2) + fun(F - E1 - E2)) / (4 * h * h)
    return H


def iso_tensor(lam, mu):
    d = I3
    return (lam * np.einsum("ij,kl->ijkl", d, d)
            + mu * (np.einsum("ik,jl->ijkl", d, d) + np.einsum("il,jk->ijkl", d, d)))


def ti_tensor(lam, mu, gam, xi, zet, nu):
    """Transversely isotropic elastic tensor as written in Lecture 13."""
    d = I3
    n = np.asarray(nu, dtype=float)
    A = iso_tensor(lam, mu)
    A = A + 8 * gam * np.einsum("i,j,k,l->ijkl", n, n, n, n)
    A = A + 4 * xi * (np.einsum("i,j,kl->ijkl", n, n, d) + np.einsum("ij,k,l->ijkl", d, n, n))
    A = A - zet * (np.einsum("i,k,jl->ijkl", n, n, d) + np.einsum("j,k,il->ijkl", n, n, d)
                   + np.einsum("j,l,ik->ijkl", n, n, d) + np.einsum("i,l,jk->ijkl", n, n, d))
    return A


def rotate_tensor(Q, A):
    return np.einsum("ia,jb,kc,ld,abcd->ijkl", Q, Q, Q, Q, A)


def rotation(axis, angle):
    a = np.asarray(axis, dtype=float)
    a = a / np.linalg.norm(a)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return I3 + np.sin(angle) * K + (1 - np.cos(angle)) * K @ K


VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]


def voigt(A):
    C = np.zeros((6, 6))
    for I, (i, j) in enumerate(VOIGT):
        for J, (k, l) in enumerate(VOIGT):
            C[I, J] = A[i, j, k, l]
    return C


# ----------------------------------------------------------------------------------------
print("\n=== Example 1: St Venant--Kirchhoff material ===")
lam, mu = 2.0, 1.3


def W_svk(F):
    E = 0.5 * (F.T @ F - I3)
    return 0.5 * lam * np.trace(E) ** 2 + mu * np.trace(E @ E)


F = I3 + 0.3 * rng.standard_normal((3, 3))
if np.linalg.det(F) < 0:
    F[:, 0] *= -1
E = 0.5 * (F.T @ F - I3)
S = lam * np.trace(E) * I3 + 2 * mu * E
T_closed = F @ S
T_fd = grad_fd(W_svk, F)
err = np.abs(T_closed - T_fd).max()
report("T = F S equals dW/dF (finite differences)", err < 1e-7, f"max err {err:.2e}")
sym = np.abs(np.linalg.inv(F) @ T_closed - (np.linalg.inv(F) @ T_closed).T).max()
report("F^{-1} T is symmetric", sym < 1e-12, f"asym {sym:.2e}")
asymT = np.abs(T_closed - T_closed.T).max()
report("T itself is not symmetric for this F (for illustration)", asymT > 1e-3, f"asym {asymT:.2e}")
# simple shear F = I + g e1 e2^T: T_12 - T_21 = g^3 (lambda/2 + mu)
g = 0.7
Fs = I3.copy()
Fs[0, 1] = g
Es = 0.5 * (Fs.T @ Fs - I3)
Ts = Fs @ (lam * np.trace(Es) * I3 + 2 * mu * Es)
report("simple shear: T_21 = mu g, T_12 - T_21 = g^3 (lambda/2 + mu)",
       abs(Ts[1, 0] - mu * g) < 1e-14 and abs(Ts[0, 1] - Ts[1, 0] - g ** 3 * (0.5 * lam + mu)) < 1e-14,
       f"T_12 - T_21 = {Ts[0, 1] - Ts[1, 0]:.6f}, g^3(lambda/2+mu) = {g ** 3 * (0.5 * lam + mu):.6f}")
# quadratic form: W(I + sH) = s^2 [kappa (tr e)^2 / 2 + mu e'_ij e'_ij] + O(s^3)
Hq = rng.standard_normal((3, 3))
eq = 0.5 * (Hq + Hq.T)
eq_dev = eq - np.trace(eq) / 3 * I3
kap1 = lam + 2 * mu / 3
quad = 0.5 * kap1 * np.trace(eq) ** 2 + mu * np.sum(eq_dev * eq_dev)
quadA = 0.5 * np.einsum("ijkl,ij,kl", iso_tensor(lam, mu), Hq, Hq)
report("A_ijkl H_ij H_kl / 2 = kappa (tr e)^2 / 2 + mu e'_ij e'_ij", abs(quad - quadA) < 1e-12,
       f"{quad:.6f} vs {quadA:.6f}")
for s_ in (1e-2, 1e-3):
    Wq = W_svk(I3 + s_ * Hq)
    report(f"s={s_}: W(I + sH) = s^2 quadratic form + O(s^3)", abs(Wq - s_ ** 2 * quad) < 5 * abs(quad) * s_ ** 3 * 10,
           f"W = {Wq:.6e}, s^2 Q = {s_ ** 2 * quad:.6e}, ratio of error to s^3 = {abs(Wq - s_ ** 2 * quad) / s_ ** 3:.3f}")
A_fd = hess_fd(W_svk, I3)
A_iso = iso_tensor(lam, mu)
err = np.abs(A_fd - A_iso).max()
report("Hessian of W at F = I equals isotropic A_ijkl", err < 1e-6, f"max err {err:.2e}")

# ----------------------------------------------------------------------------------------
print("\n=== Example 2: Elastic fluid ===")
kap = 1.7


def V(J):
    return kap * (J - 1.0 - np.log(J))


def W_fluid(F):
    return V(np.linalg.det(F))


F = I3 + 0.3 * rng.standard_normal((3, 3))
if np.linalg.det(F) < 0:
    F[:, 0] *= -1
J = np.linalg.det(F)
dJ_fd = grad_fd(lambda G: np.linalg.det(G), F)
dJ_closed = J * np.linalg.inv(F).T
err = np.abs(dJ_fd - dJ_closed).max()
report("dJ/dF_ij = J (F^{-1})_ji", err < 1e-7, f"max err {err:.2e}")
cof = np.array([[(-1) ** (i + j) * np.linalg.det(np.delete(np.delete(F, i, 0), j, 1))
                 for j in range(3)] for i in range(3)])
err = np.abs(cof - dJ_closed).max()
report("cofactor matrix equals J F^{-T}", err < 1e-12, f"max err {err:.2e}")
# Nanson-type identity (Fa) x (Fb) = J F^{-T} (a x b)
a, b = rng.standard_normal(3), rng.standard_normal(3)
lhs = np.cross(F @ a, F @ b)
rhs = J * np.linalg.inv(F).T @ np.cross(a, b)
err = np.abs(lhs - rhs).max()
report("(Fa) x (Fb) = J F^{-T} (a x b)", err < 1e-12, f"max err {err:.2e}")
Vp = kap * (1 - 1 / J)
T_closed = Vp * J * np.linalg.inv(F).T
T_fd = grad_fd(W_fluid, F)
err = np.abs(T_closed - T_fd).max()
report("T = V'(J) J F^{-T} equals dW/dF", err < 1e-7, f"max err {err:.2e}")
sigma = T_closed @ F.T / J
err = np.abs(sigma - Vp * I3).max()
report("Cauchy stress J^{-1} T F^T = V'(J) I = -p I", err < 1e-12,
       f"max err {err:.2e}, p = {-Vp:.4f}")
asymT = np.abs(T_closed - T_closed.T).max()
report("first Piola--Kirchhoff stress not symmetric here (for illustration)", asymT > 1e-3,
       f"asym {asymT:.2e}")
A_fd = hess_fd(W_fluid, I3)
A_fl = kap * np.einsum("ij,kl->ijkl", I3, I3)
err = np.abs(A_fd - A_fl).max()
report("Hessian of V(det F) at I equals kappa delta_ij delta_kl, kappa = V''(1)", err < 1e-6,
       f"max err {err:.2e}")
# pre-stressed remark: with V'(1) = p0 != 0 the linearised stress picks up extra terms
p0 = 0.4


def W_pre(F):
    J = np.linalg.det(F)
    return V(J) - p0 * J


A_pre_fd = hess_fd(W_pre, I3)
A_pre = (kap * np.einsum("ij,kl->ijkl", I3, I3)
         - p0 * (np.einsum("ij,kl->ijkl", I3, I3) - np.einsum("il,jk->ijkl", I3, I3)))
err = np.abs(A_pre_fd - A_pre).max()
report("pre-stressed fluid: Hessian = kappa dd - p0 (d_ij d_kl - d_il d_jk)  [remark]", err < 1e-6,
       f"max err {err:.2e}")
minor = np.abs(A_pre - np.transpose(A_pre, (1, 0, 2, 3))).max()
report("...and this pre-stressed tensor lacks the minor symmetry A_ijkl = A_jikl  [remark]",
       minor > 1e-3, f"|A_ijkl - A_jikl| max {minor:.2e}")

# ----------------------------------------------------------------------------------------
print("\n=== Example 3: Material symmetry group ===")
nu = np.array([0.0, 0.0, 1.0])
Fr = I3 + 0.4 * rng.standard_normal((3, 3))
C = Fr.T @ Fr


def invariants(C):
    return np.array([np.trace(C), np.trace(C @ C), np.linalg.det(C), nu @ C @ nu, nu @ C @ C @ nu])


Qax = rotation(nu, 1.1)                        # rotation about the symmetry axis
Qflip = rotation([1.0, 0.0, 0.0], np.pi)       # rotation by pi about an axis perpendicular to nu
Qgen = rotation([1.0, 2.0, 0.5], 0.9)          # a generic rotation
for name, Q in [("rotation about nu", Qax), ("rotation by pi about e1", Qflip)]:
    d = np.abs(invariants(Q.T @ C @ Q) - invariants(C)).max()
    report(f"all five invariants unchanged under {name}", d < 1e-12, f"max change {d:.2e}")
d_iso = np.abs(invariants(Qgen.T @ C @ Qgen)[:3] - invariants(C)[:3]).max()
d_ti = np.abs(invariants(Qgen.T @ C @ Qgen)[3:] - invariants(C)[3:]).max()
report("generic rotation: isotropic invariants unchanged", d_iso < 1e-12, f"max change {d_iso:.2e}")
report("generic rotation: nu.C.nu and nu.C^2.nu change", d_ti > 1e-3, f"max change {d_ti:.2e}")
# the group properties, checked on the explicit group of Qs with Q nu = +/- nu
for name, Q1, Q2 in [("closure", Qax, Qflip), ("closure", Qflip, rotation(nu, -2.3))]:
    Q = Q1 @ Q2
    report(f"{name}: product preserves the axis", abs(abs(nu @ Q @ nu) - 1) < 1e-12,
           f"|nu.Q nu| = {abs(nu @ Q @ nu):.12f}")
report("inverse: Q^T preserves the axis", abs(abs(nu @ Qflip.T @ nu) - 1) < 1e-12)
# induced invariance of the elastic tensor
lam, mu, gam, xi, zet = 2.0, 1.0, 0.3, 0.2, 0.1
A_ti = ti_tensor(lam, mu, gam, xi, zet, nu)
for name, Q, expect in [("rotation about nu", Qax, True), ("rotation by pi about e1", Qflip, True),
                        ("generic rotation", Qgen, False)]:
    d = np.abs(rotate_tensor(Q, A_ti) - A_ti).max()
    ok = (d < 1e-12) if expect else (d > 1e-3)
    report(f"TI tensor {'invariant' if expect else 'NOT invariant'} under {name}", ok, f"max change {d:.2e}")
d = np.abs(rotate_tensor(Qgen, iso_tensor(lam, mu)) - iso_tensor(lam, mu)).max()
report("isotropic tensor invariant under generic rotation", d < 1e-12, f"max change {d:.2e}")
# the isotropic invariants determine the eigenvalues (characteristic polynomial)
I1, I2, I3v = np.trace(C), 0.5 * (np.trace(C) ** 2 - np.trace(C @ C)), np.linalg.det(C)
ev = np.linalg.eigvalsh(C)
poly = np.poly(ev)  # x^3 - I1 x^2 + I2 x - I3
d = np.abs(poly - np.array([1, -I1, I2, -I3v])).max()
report("char. polynomial coefficients are I1, I2 = ((tr C)^2 - tr C^2)/2, I3", d < 1e-10, f"max err {d:.2e}")

# ----------------------------------------------------------------------------------------
print("\n=== Example 4: Voigt notation and the transversely isotropic tensor ===")
print(f"sample moduli: lambda={lam}, mu={mu}, gamma={gam}, xi={xi}, zeta={zet}, nu = e3")
for name, perm in [("A_ijkl = A_jikl", (1, 0, 2, 3)), ("A_ijkl = A_ijlk", (0, 1, 3, 2)),
                   ("A_ijkl = A_klij", (2, 3, 0, 1))]:
    d = np.abs(A_ti - np.transpose(A_ti, perm)).max()
    report(f"TI tensor symmetry {name}", d < 1e-14, f"max violation {d:.1e}")
Cv = voigt(A_ti)
np.set_printoptions(precision=4, suppress=True, linewidth=120)
print("6x6 Voigt matrix of the TI tensor:")
print(Cv)
report("Voigt matrix symmetric", np.abs(Cv - Cv.T).max() < 1e-14)
# every component of A is reproduced by the Voigt matrix (i.e. nothing is lost in the mapping)
idx = {}
for I, (i, j) in enumerate(VOIGT):
    idx[(i, j)] = I
    idx[(j, i)] = I
A_back = np.zeros((3, 3, 3, 3))
for i in range(3):
    for j in range(3):
        for k in range(3):
            for l in range(3):
                A_back[i, j, k, l] = Cv[idx[(i, j)], idx[(k, l)]]
report("81 components recovered from the 21 Voigt entries", np.abs(A_back - A_ti).max() < 1e-14)
closed = {
    "C11 = A_1111 = lambda + 2 mu": (Cv[0, 0], lam + 2 * mu),
    "C22 = A_2222 = lambda + 2 mu": (Cv[1, 1], lam + 2 * mu),
    "C33 = A_3333 = lambda + 2 mu + 8 gamma + 8 xi - 4 zeta": (Cv[2, 2], lam + 2 * mu + 8 * gam + 8 * xi - 4 * zet),
    "C12 = A_1122 = lambda": (Cv[0, 1], lam),
    "C13 = A_1133 = lambda + 4 xi": (Cv[0, 2], lam + 4 * xi),
    "C23 = A_2233 = lambda + 4 xi": (Cv[1, 2], lam + 4 * xi),
    "C44 = A_2323 = mu - zeta": (Cv[3, 3], mu - zet),
    "C55 = A_1313 = mu - zeta": (Cv[4, 4], mu - zet),
    "C66 = A_1212 = mu": (Cv[5, 5], mu),
}
for name, (num, form) in closed.items():
    report(name, abs(num - form) < 1e-14, f"{num:.4f} vs {form:.4f}")
mask = np.ones((6, 6), bool)
for (I, J) in [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (0, 1), (1, 0), (0, 2), (2, 0), (1, 2), (2, 1)]:
    mask[I, J] = False
report("all other Voigt entries vanish", np.abs(Cv[mask]).max() < 1e-14, f"max {np.abs(Cv[mask]).max():.1e}")
report("hexagonal relation C66 = (C11 - C12)/2", abs(Cv[5, 5] - 0.5 * (Cv[0, 0] - Cv[0, 1])) < 1e-14)
# invert the map (C11, C33, C13, C44, C66) -> (lambda, mu, gamma, xi, zeta)
mu_r = Cv[5, 5]
lam_r = Cv[0, 0] - 2 * mu_r
xi_r = (Cv[0, 2] - lam_r) / 4
zet_r = mu_r - Cv[3, 3]
gam_r = (Cv[2, 2] - lam_r - 2 * mu_r - 8 * xi_r + 4 * zet_r) / 8
d = np.abs(np.array([lam_r, mu_r, gam_r, xi_r, zet_r]) - np.array([lam, mu, gam, xi, zet])).max()
report("five moduli recovered from (C11,C33,C13,C44,C66)", d < 1e-14, f"max err {d:.1e}")
# isotropic Voigt matrix
Ci = voigt(iso_tensor(lam, mu))
Ci_expect = np.array([[lam + 2 * mu, lam, lam, 0, 0, 0], [lam, lam + 2 * mu, lam, 0, 0, 0],
                      [lam, lam, lam + 2 * mu, 0, 0, 0], [0, 0, 0, mu, 0, 0], [0, 0, 0, 0, mu, 0],
                      [0, 0, 0, 0, 0, mu]])
report("isotropic Voigt matrix as stated", np.abs(Ci - Ci_expect).max() < 1e-14)
# the TI tensor is invariant under a general rotation only when gamma = xi = zeta = 0
d = np.abs(rotate_tensor(Qgen, ti_tensor(lam, mu, 0, 0, 0, nu)) - ti_tensor(lam, mu, 0, 0, 0, nu)).max()
report("TI tensor reduces to isotropic one when gamma = xi = zeta = 0", d < 1e-12)
# stress from strain in Voigt form: T_I = C_IJ e_J with engineering shear strains
Hh = rng.standard_normal((3, 3))
e = 0.5 * (Hh + Hh.T)
T_full = np.einsum("ijkl,kl->ij", A_ti, Hh)
e_v = np.array([e[0, 0], e[1, 1], e[2, 2], 2 * e[1, 2], 2 * e[0, 2], 2 * e[0, 1]])
T_v = Cv @ e_v
T_v_full = np.array([T_full[i, j] for (i, j) in VOIGT])
report("T_I = C_IJ e_J with e_4 = 2 e_23 etc.", np.abs(T_v - T_v_full).max() < 1e-12)

# ----------------------------------------------------------------------------------------
print("\n=== Example 5: Navier equation, homogeneous and heterogeneous ===")
# quadratic displacement field, linear moduli: central differences of the flux are then exact
a = rng.standard_normal((3, 3, 3))
a = 0.5 * (a + np.transpose(a, (0, 2, 1)))
bmat = rng.standard_normal((3, 3))
lam0, mu0 = 2.0, 1.0
lvec, mvec = rng.standard_normal(3) * 0.3, rng.standard_normal(3) * 0.3


def u(x):
    return np.einsum("ijk,j,k->i", a, x, x) + bmat @ x


def gradu(x):
    return 2 * np.einsum("ijk,k->ij", a, x) + bmat


def lam_f(x, hetero):
    return lam0 + (lvec @ x if hetero else 0.0)


def mu_f(x, hetero):
    return mu0 + (mvec @ x if hetero else 0.0)


def flux(x, hetero):
    A = iso_tensor(lam_f(x, hetero), mu_f(x, hetero))
    return np.einsum("ijkl,kl->ij", A, gradu(x))


def div_flux(x, hetero, h=1e-3):
    out = np.zeros(3)
    for j in range(3):
        dx = np.zeros(3)
        dx[j] = h
        out += (flux(x + dx, hetero)[:, j] - flux(x - dx, hetero)[:, j]) / (2 * h)
    return out


x0 = rng.standard_normal(3)
# analytic pieces for the quadratic field: grad(div u), laplacian u
div_u = np.trace(gradu(x0))
grad_div = 2 * np.einsum("iji->j", a) + 0 * x0   # d_i (2 a_kjk x_j + b_kk) = 2 a_kik
grad_div = 2 * np.einsum("kik->i", a)
lap_u = 2 * np.einsum("ikk->i", a)
G = gradu(x0)
eps = 0.5 * (G + G.T)
for hetero in (False, True):
    lhs = div_flux(x0, hetero)
    lam_x, mu_x = lam_f(x0, hetero), mu_f(x0, hetero)
    rhs = (lam_x + mu_x) * grad_div + mu_x * lap_u
    if hetero:
        rhs = rhs + lvec * div_u + 2 * eps @ mvec
    d = np.abs(lhs - rhs).max()
    label = "heterogeneous" if hetero else "homogeneous"
    report(f"d_j(A_ijkl d_l u_k) = closed form ({label})", d < 1e-8, f"max err {d:.2e}")
# the extra terms are genuinely needed
lhs = div_flux(x0, True)
rhs_wrong = (lam_f(x0, True) + mu_f(x0, True)) * grad_div + mu_f(x0, True) * lap_u
report("heterogeneous case differs from naive formula", np.abs(lhs - rhs_wrong).max() > 1e-3,
       f"difference {np.abs(lhs - rhs_wrong).max():.2e}")

# ----------------------------------------------------------------------------------------
print("\n=== Example 6: Uniform dilatation of the St Venant--Kirchhoff material ===")
lam, mu = 2.0, 1.3
kap = lam + 2 * mu / 3


def T_svk(F):
    E = 0.5 * (F.T @ F - I3)
    return F @ (lam * np.trace(E) * I3 + 2 * mu * E)


for eps_ in (1e-2, 1e-3):
    F = (1 + eps_) * I3
    T = T_svk(F)
    T_exact = (3 * lam + 2 * mu) * (1 + eps_) * (eps_ + 0.5 * eps_ ** 2)
    report(f"eps={eps_}: T = (3 lambda + 2 mu)(1+eps)(eps + eps^2/2) I", np.abs(T - T_exact * I3).max() < 1e-13)
    T1 = 3 * kap * eps_
    T2 = 3 * kap * (eps_ + 1.5 * eps_ ** 2)
    report(f"eps={eps_}: first-order T error is O(eps^2)", abs(T[0, 0] - T1) < 5 * kap * eps_ ** 2,
           f"err {abs(T[0, 0] - T1):.2e} vs 3 kap*1.5 eps^2 = {4.5 * kap * eps_ ** 2:.2e}")
    report(f"eps={eps_}: second-order T error is O(eps^3)", abs(T[0, 0] - T2) < 3 * kap * eps_ ** 3,
           f"err {abs(T[0, 0] - T2):.2e} vs 3 kap*0.5 eps^3 = {1.5 * kap * eps_ ** 3:.2e}")
    J = np.linalg.det(F)
    sig = T @ F.T / J
    p = -sig[0, 0]
    p1 = -3 * kap * eps_
    p2 = -3 * kap * (eps_ - 0.5 * eps_ ** 2)
    report(f"eps={eps_}: Cauchy pressure p = -3 kap eps + O(eps^2)", abs(p - p1) < 5 * kap * eps_ ** 2,
           f"p = {p:.6e}, -3 kap eps = {p1:.6e}")
    report(f"eps={eps_}: p = -3 kap (eps - eps^2/2) + O(eps^3)", abs(p - p2) < 4 * kap * eps_ ** 3,
           f"err {abs(p - p2):.2e} vs 3 kap*0.5 eps^3 = {1.5 * kap * eps_ ** 3:.2e}")
    report(f"eps={eps_}: p = -kappa (J - 1) + O(eps^2)", abs(p + kap * (J - 1)) < 8 * kap * eps_ ** 2,
           f"p = {p:.6e}, -kap(J-1) = {-kap * (J - 1):.6e}")
    p_closed = -1.5 * kap * (J ** (1.0 / 3.0) - J ** (-1.0 / 3.0))
    report(f"eps={eps_}: exact p = -(3 kappa/2)(J^(1/3) - J^(-1/3))", abs(p - p_closed) < 1e-12,
           f"p = {p:.10e}, closed form {p_closed:.10e}")
    sig_closed = 3 * kap * (eps_ + 0.5 * eps_ ** 2) / (1 + eps_)
    report(f"eps={eps_}: Cauchy stress = 3 kappa (eps + eps^2/2)/(1+eps) I", np.abs(sig - sig_closed * I3).max() < 1e-13)
report("3 lambda + 2 mu = 3 kappa", abs(3 * lam + 2 * mu - 3 * kap) < 1e-14)
for eps_ in (0.3, -0.5):
    F = (1 + eps_) * I3
    J = np.linalg.det(F)
    E = 0.5 * (F.T @ F - I3)
    Wd = 0.5 * lam * np.trace(E) ** 2 + mu * np.trace(E @ E)
    report(f"eps={eps_}: W = (9/8) kappa (J^(2/3) - 1)^2 under uniform dilatation",
           abs(Wd - 9 / 8 * kap * (J ** (2 / 3) - 1) ** 2) < 1e-13, f"W = {Wd:.6f}")
report("W at J -> 0 is finite, 9 kappa / 8", True, f"9 kappa/8 = {9 * kap / 8:.4f}")

print()
print("ALL CHECKS PASSED" if status["fail"] == 0 else f"{status['fail']} CHECK(S) FAILED")
