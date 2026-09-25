"""
Numerical checks for the worked examples accompanying Lecture 18
(Delay time tomography): examples/examples7.tex.

Every numerical or closed-form claim made in the solutions is checked here.
Run with:  python verify_examples7.py
"""
import numpy as np
from scipy.optimize import brentq, minimize

np.set_printoptions(precision=6, suppress=True)
status = {"fail": 0}


def check(name, val, ref, tol=1e-9, rel=False):
    val = np.asarray(val, dtype=float)
    ref = np.asarray(ref, dtype=float)
    err = np.max(np.abs(val - ref))
    if rel:
        err = err / max(np.max(np.abs(ref)), 1e-300)
    ok = err <= tol
    if not ok:
        status["fail"] += 1
    print(f"{'OK  ' if ok else 'FAIL'} {name}: value={val.tolist()} ref={ref.tolist()} err={err:.3e}")
    return ok


# ---------------------------------------------------------------------------
# Common parameters
# ---------------------------------------------------------------------------
ell = 100.0      # block length, km
alpha = 8.0      # reference P velocity, km/s

print("\n=== Example 1: toy tomography with two blocks ===")
# case (a): one ray through both blocks
Aa = np.array([[ell, ell]])
null = np.array([1.0, -1.0])
check("E1(a) A m0 = 0 for m0=(1,-1)", Aa @ null, [0.0])
# the null space is one-dimensional
check("E1(a) rank A = 1", np.linalg.matrix_rank(Aa), 1)

# case (b): add ray through block 1 only
Ab = np.array([[ell, ell], [ell, 0.0]])
check("E1(b) det A = -ell^2", np.linalg.det(Ab), -ell**2, tol=1e-6)
Ab_inv = np.array([[0.0, 1.0], [1.0, -1.0]]) / ell
check("E1(b) A^{-1} formula", Ab_inv, np.linalg.inv(Ab), tol=1e-12)
d_b = np.array([1.2, 0.5])
m_b = np.linalg.solve(Ab, d_b)
check("E1(b) m = (d2/ell, (d1-d2)/ell) = (0.005, 0.007)", m_b, [0.005, 0.007], tol=1e-12)
dalpha_b = -alpha**2 * m_b
check("E1(b) delta alpha = (-0.32, -0.448) km/s", dalpha_b, [-0.32, -0.448], tol=1e-12)
check("E1(b) delta alpha/alpha = (-4%, -5.6%)", dalpha_b / alpha, [-0.04, -0.056], tol=1e-12)
# model covariance for unit-variance data: A^{-1} A^{-T} = (1/ell^2) [[1,-1],[-1,2]]
check("E1(b) A^{-1}A^{-T} = (1/ell^2)[[1,-1],[-1,2]]", Ab_inv @ Ab_inv.T,
      np.array([[1, -1], [-1, 2]]) / ell**2, tol=1e-15)

# case (c): third ray through block 2 only, inconsistent data
Ac = np.array([[ell, ell], [ell, 0.0], [0.0, ell]])
d_c = np.array([1.2, 0.5, 0.4])
AtA = Ac.T @ Ac
check("E1(c) A^T A = ell^2 [[2,1],[1,2]]", AtA, ell**2 * np.array([[2, 1], [1, 2]]))
AtA_inv = np.array([[2, -1], [-1, 2]]) / (3 * ell**2)
check("E1(c) (A^T A)^{-1} formula", AtA_inv, np.linalg.inv(AtA), tol=1e-15)
Atd = Ac.T @ d_c
check("E1(c) A^T d = ell (d1+d2, d1+d3)", Atd, ell * np.array([1.7, 1.6]))
m_c = AtA_inv @ Atd
m_c_formula = np.array([d_c[0] + 2 * d_c[1] - d_c[2], d_c[0] - d_c[1] + 2 * d_c[2]]) / (3 * ell)
check("E1(c) closed form for m_LS", m_c, m_c_formula, tol=1e-15)
check("E1(c) m_LS = (0.006, 0.005) s/km", m_c, [0.006, 0.005], tol=1e-12)
m_lstsq = np.linalg.lstsq(Ac, d_c, rcond=None)[0]
check("E1(c) agrees with numpy lstsq", m_c, m_lstsq, tol=1e-12)
r_c = Ac @ m_c - d_c
check("E1(c) predicted data (1.1, 0.6, 0.5)", Ac @ m_c, [1.1, 0.6, 0.5], tol=1e-12)
check("E1(c) residual = -(d1-d2-d3)/3 (1,-1,-1) = (-0.1,0.1,0.1)", r_c, [-0.1, 0.1, 0.1], tol=1e-12)
check("E1(c) residual orthogonal to columns of A", Ac.T @ r_c, [0, 0], tol=1e-12)
check("E1(c) delta alpha/alpha = (-4.8%, -4.0%)", -alpha * m_c, [-0.048, -0.040], tol=1e-12)

# ---------------------------------------------------------------------------
print("\n=== Example 2: regularised least squares (case (a), B = I, C = I) ===")
d = 1.2
A = Aa
for lam_over_ell2, mj_ref, mis_ref in [(0.0, 0.006, 0.0), (0.5, 0.0048, 0.24),
                                        (2.0, 0.003, 0.6), (8.0, 0.0012, 0.96)]:
    lam = lam_over_ell2 * ell**2
    M = A.T @ A + lam * np.eye(2)
    if lam > 0:
        m_lam = np.linalg.solve(M, A.T @ np.array([d]))
    else:
        m_lam = np.linalg.pinv(A) @ np.array([d])   # minimum-norm solution
    m_formula = ell * d / (2 * ell**2 + lam) * np.array([1.0, 1.0])
    check(f"E2 m(lambda) closed form, lambda/ell^2={lam_over_ell2}", m_lam, m_formula, tol=1e-12)
    check(f"E2 m_j = {mj_ref}, lambda/ell^2={lam_over_ell2}", m_lam, [mj_ref, mj_ref], tol=1e-12)
    misfit = abs(float((A @ m_lam)[0]) - d)
    check(f"E2 |Am-d| = lambda d/(2ell^2+lambda) = {mis_ref}", misfit, mis_ref, tol=1e-12)
    check(f"E2 |Am-d| = {mis_ref} (closed form)", lam * d / (2 * ell**2 + lam), mis_ref, tol=1e-12)
    norm_m = np.linalg.norm(m_lam)
    check(f"E2 ||m|| = sqrt2 ell d/(2ell^2+lambda), lambda/ell^2={lam_over_ell2}",
          norm_m, np.sqrt(2) * ell * d / (2 * ell**2 + lam), tol=1e-12)
    # straight-line trade-off curve
    check(f"E2 trade-off line |Am-d| = d - sqrt2 ell ||m||, lambda/ell^2={lam_over_ell2}",
          misfit, d - np.sqrt(2) * ell * norm_m, tol=1e-12)
# determinant and inverse of A^T A + lambda I
lam = 0.37 * ell**2
M = A.T @ A + lam * np.eye(2)
check("E2 det(A^T A + lambda I) = lambda(2 ell^2 + lambda)", np.linalg.det(M), lam * (2 * ell**2 + lam), tol=1e-6, rel=True)
Minv = np.array([[ell**2 + lam, -ell**2], [-ell**2, ell**2 + lam]]) / (lam * (2 * ell**2 + lam))
check("E2 inverse formula", Minv, np.linalg.inv(M), tol=1e-12, rel=True)
# eigen-decomposition of A^T A
w, V = np.linalg.eigh(A.T @ A)
check("E2 eigenvalues of A^T A are 0 and 2 ell^2", np.sort(w), [0.0, 2 * ell**2], tol=1e-9)
# large-lambda behaviour
lam_big = 1e8 * ell**2
m_big = np.linalg.solve(A.T @ A + lam_big * np.eye(2), A.T @ np.array([d]))
check("E2 m -> 0 as lambda -> infinity (lambda = 1e8 ell^2)", m_big, [0, 0], tol=1e-9)

# ---------------------------------------------------------------------------
print("\n=== Example 3: Fermat's principle checked exactly (slab, Snell's law) ===")
h = 100.0     # slab thickness, km
a = 200.0     # distance of source and receiver from the slab faces, km
b = 250.0     # half the total transverse offset, km
theta0 = np.arctan2(2 * b, 2 * a + h)   # unperturbed take-off angle from the slab normal
check("E3 unperturbed angle is 45 degrees", np.degrees(theta0), 45.0, tol=1e-12)
ell_in = h / np.cos(theta0)
T0 = (2 * a + h) / (alpha * np.cos(theta0))
check("E3 ell_in = h sqrt2 = 141.42 km", ell_in, 100 * np.sqrt(2), tol=1e-9)
check("E3 T0 = 500 sqrt2/8 = 88.388 s", T0, 500 * np.sqrt(2) / 8, tol=1e-9)
check("E3 ell_in/alpha = 17.678 s", ell_in / alpha, 17.677669529663688, tol=1e-9)


def snell_travel_time(eps):
    """Exact travel time of the ray from (-a,-b) to (h+a,b) when the slab velocity is alpha(1+eps)."""
    v2 = alpha * (1 + eps)

    def offset(p):
        s1, s2 = p * alpha, p * v2
        return 2 * a * s1 / np.sqrt(1 - s1**2) + h * s2 / np.sqrt(1 - s2**2) - 2 * b

    pmax = min(1 / alpha, 1 / v2) * (1 - 1e-12)
    p = brentq(offset, 0.0, pmax, xtol=1e-15, rtol=1e-15, maxiter=500)
    s1, s2 = p * alpha, p * v2
    return 2 * a / (alpha * np.sqrt(1 - s1**2)) + h / (v2 * np.sqrt(1 - s2**2)), p


def trial_path_time(y, eps):
    """Travel time of the three-segment path with entry point (0,y1) and exit point (h,y2)."""
    y1, y2 = y
    Tout = (np.hypot(a, y1 + b) + np.hypot(a, b - y2)) / alpha
    Tin = np.hypot(h, y2 - y1) / (alpha * (1 + eps))
    return Tout + Tin


check("E3 Snell solver reproduces T0 at eps=0", snell_travel_time(0.0)[0], T0, tol=1e-9)
y0 = np.array([-b + a * np.tan(theta0), b - a * np.tan(theta0)])   # unperturbed entry/exit points
check("E3 unperturbed entry/exit points (-50, 50) km", y0, [-50.0, 50.0], tol=1e-9)
check("E3 trial-path time at y0, eps=0 equals T0", trial_path_time(y0, 0.0), T0, tol=1e-9)

# second-order coefficients derived in the text
c_nonlin = ell_in / alpha                                   # from 1/(1+eps) = 1 - eps + eps^2 ...
c_bend = -(ell_in / alpha) * np.tan(theta0)**2 * a / (h + 2 * a)   # -1/2 g^T H^{-1} g
c_total = c_nonlin + c_bend
check("E3 bending coefficient = -7.071 s", c_bend, -7.0710678118654755, tol=1e-9)
check("E3 total second-order coefficient = 10.607 s", c_total, 10.606601717798213, tol=1e-9)

# Hessian / gradient check of the closed form -1/2 g^T H^{-1} g
def hess_T0(y, dy=0.1):
    H = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            e_i = np.zeros(2); e_i[i] = dy
            e_j = np.zeros(2); e_j[j] = dy
            H[i, j] = (trial_path_time(y + e_i + e_j, 0) - trial_path_time(y + e_i - e_j, 0)
                       - trial_path_time(y - e_i + e_j, 0) + trial_path_time(y - e_i - e_j, 0)) / (4 * dy * dy)
    return H


H_fd = hess_T0(y0)
H_formula = (np.cos(theta0)**3 / alpha) * np.array([[1 / a + 1 / h, -1 / h], [-1 / h, 1 / a + 1 / h]])
check("E3 Hessian of T0 w.r.t. (y1,y2) closed form (finite differences)", H_fd, H_formula, tol=1e-4, rel=True)
g = (np.sin(theta0) / alpha) * np.array([-1.0, 1.0])
check("E3 -1/2 g^T H^{-1} g equals bending coefficient", -0.5 * g @ np.linalg.solve(H_formula, g), c_bend, tol=1e-12)

print("  eps      T_exact        T_fixed        T_lin        (Tex-Tlin)/eps^2  (Tex-Tfix)/eps^2")
for eps in [0.001, 0.01, 0.03, 0.1]:
    Tex, p = snell_travel_time(eps)
    # cross-check: direct minimisation over the trial paths
    res = minimize(trial_path_time, y0, args=(eps,), method="Nelder-Mead",
                   options={"xatol": 1e-10, "fatol": 1e-13, "maxiter": 20000})
    check(f"E3 Snell and direct minimisation agree, eps={eps}", Tex, res.fun, tol=1e-9)
    Tfix = trial_path_time(y0, eps)
    check(f"E3 fixed-path time = T0 + (1/(1+eps)-1) ell_in/alpha, eps={eps}",
          Tfix, T0 + (1 / (1 + eps) - 1) * ell_in / alpha, tol=1e-12)
    Tlin = T0 - eps * ell_in / alpha
    print(f"  {eps:<7} {Tex:<14.6f} {Tfix:<14.6f} {Tlin:<12.6f} {(Tex-Tlin)/eps**2:<17.4f} {(Tex-Tfix)/eps**2:.4f}")
    # the exact ray always beats the fixed path (Fermat: minimum)
    check(f"E3 T_exact <= T_fixed, eps={eps}", float(Tex <= Tfix + 1e-15), 1.0)
    # second-order coefficients approached as eps -> 0 (error O(eps))
    check(f"E3 (T_exact - T_lin)/eps^2 -> 10.607 within 30 eps, eps={eps}",
          (Tex - Tlin) / eps**2, c_total, tol=30 * eps)
    check(f"E3 (T_exact - T_fixed)/eps^2 -> -7.071 within 30 eps, eps={eps}",
          (Tex - Tfix) / eps**2, c_bend, tol=30 * eps)

# table values quoted in the text (ms precision)
table = {}
for eps in [0.01, 0.03, 0.1]:
    Tex, _ = snell_travel_time(eps)
    table[eps] = (Tex - T0, -eps * ell_in / alpha, Tex - T0 + eps * ell_in / alpha, Tex - trial_path_time(y0, eps))
    print(f"  eps={eps}: dT_exact={table[eps][0]:.5f} s, dT_lin={table[eps][1]:.5f} s, "
          f"error={table[eps][2]*1e3:.3f} ms, bending part={table[eps][3]*1e3:.3f} ms")
# values quoted in the table of the text
check("E3 table eps=0.01: dT_lin = -0.17678 s", table[0.01][1], -0.17678, tol=6e-6)
check("E3 table eps=0.01: dT_exact = -0.17573 s", table[0.01][0], -0.17573, tol=6e-6)
check("E3 table eps=0.01: error = 1.05 ms, bending part = -0.70 ms",
      [table[0.01][2] * 1e3, table[0.01][3] * 1e3], [1.05, -0.70], tol=6e-3)
check("E3 table eps=0.03: dT_lin = -0.53033 s", table[0.03][1], -0.53033, tol=6e-6)
check("E3 table eps=0.03: dT_exact = -0.52114 s", table[0.03][0], -0.52114, tol=6e-6)
check("E3 table eps=0.03: error = 9.19 ms, bending part = -6.25 ms",
      [table[0.03][2] * 1e3, table[0.03][3] * 1e3], [9.19, -6.25], tol=6e-3)
check("E3 table eps=0.1: dT_lin = -1.76777 s", table[0.1][1], -1.76777, tol=6e-6)
check("E3 table eps=0.1: dT_exact = -1.67393 s", table[0.1][0], -1.67393, tol=6e-6)
check("E3 table eps=0.1: error = 93.84 ms, bending part = -66.87 ms",
      [table[0.1][2] * 1e3, table[0.1][3] * 1e3], [93.84, -66.87], tol=6e-3)
# quadratic scaling: error(0.03)/error(0.01) ~ 9, error(0.1)/error(0.01) ~ 100 (up to O(eps) corrections)
check("E3 error(0.03)/error(0.01) = 8.78 (quadratic scaling gives 9)", table[0.03][2] / table[0.01][2], 8.78, tol=0.01)
check("E3 error(0.1)/error(0.01) = 89.6 (quadratic scaling gives 100)", table[0.1][2] / table[0.01][2], 89.6, tol=0.1)
# the second-order formula itself predicts the errors to within a few percent
for eps in [0.01, 0.03, 0.1]:
    check(f"E3 second-order formula c_total eps^2 within 15% of the true error, eps={eps}",
          c_total * eps**2 / table[eps][2], 1.0, tol=0.15)
check("E3 relative error of linearised delay at eps=0.1 is 5.6% (quoted as about 6%)",
      table[0.1][2] / abs(table[0.1][0]) * 100, 5.6, tol=0.05)
check("E3 eps=0.1: |bending part| < |non-linearity part| = error - bending",
      float(abs(table[0.1][3]) < abs(table[0.1][2] - table[0.1][3])), 1.0)
check("E3 second-order formula predicts 1.06, 9.5, 106 ms",
      [c_total * e**2 * 1e3 for e in (0.01, 0.03, 0.1)], [1.06, 9.55, 106.1], tol=0.06)
# a negative perturbation too (slow slab)
Tex_neg, _ = snell_travel_time(-0.03)
check("E3 eps=-0.03: error still positive and ~ c_total eps^2",
      (Tex_neg - (T0 + 0.03 * ell_in / alpha)) / 0.03**2, c_total, tol=1.0)

# ---------------------------------------------------------------------------
print("\n=== Example 4: Bayesian view of case (a) ===")
sig_d = 0.1        # s
sig_m = 0.002      # s/km
d = 1.2
A = Aa
C_inv = np.eye(1) / sig_d**2
Cm_inv = np.eye(2) / sig_m**2
Cp = np.linalg.inv(A.T @ C_inv @ A + Cm_inv)
P_par = 0.5 * np.array([[1, 1], [1, 1]])
P_perp = 0.5 * np.array([[1, -1], [-1, 1]])
c = sig_m**2 * sig_d**2 / (sig_d**2 + 2 * ell**2 * sig_m**2)
check("E4 c = sigma_m^2/9", c, sig_m**2 / 9, tol=1e-18)
check("E4 C_p = sigma_m^2 P_perp + c P_par", Cp, sig_m**2 * P_perp + c * P_par, tol=1e-12, rel=True)
check("E4 C_p explicit [[2.222e-6,-1.778e-6],[-1.778e-6,2.222e-6]]", Cp,
      [[2.2222222e-6, -1.7777778e-6], [-1.7777778e-6, 2.2222222e-6]], tol=1e-6, rel=True)
e_par = np.array([1, 1]) / np.sqrt(2)
e_perp = np.array([1, -1]) / np.sqrt(2)
check("E4 posterior variance along (1,1)/sqrt2 = sigma_m^2/9", e_par @ Cp @ e_par, sig_m**2 / 9, tol=1e-18)
check("E4 posterior variance along (1,-1)/sqrt2 = sigma_m^2 (unchanged)", e_perp @ Cp @ e_perp, sig_m**2, tol=1e-18)
check("E4 posterior std along (1,1)/sqrt2 = sigma_m/3", np.sqrt(e_par @ Cp @ e_par), sig_m / 3, tol=1e-12)
rho = Cp[0, 1] / np.sqrt(Cp[0, 0] * Cp[1, 1])
check("E4 posterior correlation coefficient = -0.8", rho, -0.8, tol=1e-12)
check("E4 rho = (c - sigma_m^2)/(c + sigma_m^2)", rho, (c - sig_m**2) / (c + sig_m**2), tol=1e-12)
check("E4 marginal posterior std of each block = 1.49e-3 s/km", np.sqrt(Cp[0, 0]), 1.4907e-3, tol=1e-6)
# posterior mean
mp = Cp @ (A.T @ C_inv @ np.array([d]))
lam_eq = sig_d**2 / sig_m**2
check("E4 lambda = sigma_d^2/sigma_m^2 = 2500 km^2", lam_eq, 2500.0)
check("E4 m_p = ell d/(2 ell^2 + sigma_d^2/sigma_m^2) (1,1)", mp, ell * d / (2 * ell**2 + lam_eq) * np.array([1, 1]), tol=1e-15)
check("E4 m_p = (0.005333, 0.005333) s/km", mp, [0.0053333333, 0.0053333333], tol=1e-9)
m_reg = np.linalg.solve(A.T @ A + lam_eq * np.eye(2), A.T @ np.array([d]))
check("E4 m_p equals regularised LS solution with lambda = sigma_d^2/sigma_m^2", mp, m_reg, tol=1e-15)
check("E4 delta alpha/alpha = -4.3% in each block", -alpha * mp, [-0.042667, -0.042667], tol=1e-6)
check("E4 posterior std of (m1+m2)/sqrt2 is 3 times smaller than prior", sig_m / np.sqrt(e_par @ Cp @ e_par), 3.0, tol=1e-9)
check("E4 marginal std reduced from 2.0e-3 to 1.5e-3", np.sqrt(Cp[0, 0]) * 1e3, 1.5, tol=0.01)
# resolution matrix
R = Cp @ A.T @ C_inv @ A
check("E4 resolution matrix R = (8/9) P_par = (4/9)[[1,1],[1,1]]", R, (4 / 9) * np.array([[1, 1], [1, 1]]), tol=1e-12)
check("E4 R (1,1) = (8/9)(1,1)", R @ np.array([1, 1]), [8 / 9, 8 / 9], tol=1e-12)
check("E4 R (1,-1) = 0", R @ np.array([1, -1]), [0, 0], tol=1e-12)
# limit of perfect data
Cp0 = np.linalg.inv(A.T @ A / 1e-4**2 + Cm_inv)
check("E4 sigma_d -> 0 (1e-4 s): C_p -> sigma_m^2 P_perp", Cp0, sig_m**2 * P_perp, tol=1e-6, rel=True)

# ---------------------------------------------------------------------------
print("\n=== Example 5: sensitivity along the ray ===")
v = 8.0
R_s = 50.0
delta = -0.08
dT = -2 * R_s * delta / v**2
check("E5 dT = -2 R delta / v^2 = 0.125 s", dT, 0.125, tol=1e-12)


def dT_sphere(offset, nquad=200001):
    """Ray-theoretic delay for a sphere of radius R_s centred a distance `offset` from a straight ray."""
    s = np.linspace(-2 * R_s, 2 * R_s, nquad)
    inside = (s**2 + offset**2) < R_s**2
    dv = np.where(inside, delta, 0.0)
    return -np.trapz(dv / v**2, s)


check("E5 quadrature along the ray reproduces 2R delta/v^2 (offset 0)", dT_sphere(0.0), dT, tol=1e-4, rel=True)
for off, chord in [(30.0, 2 * np.sqrt(R_s**2 - 30.0**2)), (50.0 + 1e-9, 0.0), (80.0, 0.0)]:
    check(f"E5 offset {off:.0f} km: chord length {chord:.3f} km", dT_sphere(off), -chord * delta / v**2, tol=1e-4)
check("E5 offset 30 km: dT = 0.1 s", dT_sphere(30.0), 0.1, tol=1e-4)
# volume of the anomaly scales as R^3 while dT scales as R
for Rr in [50.0, 25.0, 12.5]:
    dTr = -2 * Rr * delta / v**2
    check(f"E5 dT proportional to R (R={Rr} km): {dTr:.5f} s", dTr / Rr, dT / R_s, tol=1e-15)
check("E5 volume ratio (R=50 vs R=25) is 8, delay ratio is 2", [(50 / 25)**3, 2.0], [8.0, 2.0])
check("E5 velocity perturbation is -1%", delta / v, -0.01, tol=1e-15)

print("\n=== summary ===")
print("ALL OK" if status["fail"] == 0 else f"{status['fail']} FAILURES")
