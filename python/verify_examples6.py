"""Numerical checks for examples/examples6.tex (Worked examples for Lecture 17).

Every numerical or closed-form claim made in the solutions is checked here.
Run with:  python3 verify_examples6.py
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar, brentq

np.set_printoptions(precision=6, suppress=True)
failures = 0


def check(name, got, expected, tol=1e-6, rel=False):
    global failures
    err = abs(got - expected)
    scale = max(abs(expected), 1e-300) if rel else 1.0
    ok = err <= tol * scale
    if not ok:
        failures += 1
    print(f"{'OK  ' if ok else 'FAIL'} {name}: got {got!r}, expected {expected!r}, "
          f"err {err:.3e}")
    return ok


# ---------------------------------------------------------------------------
# Generic spherical integrals (lecture eqs. for Delta and T), by quadrature.
# The integrable 1/sqrt singularity at r_t is removed with r = r_t + s^2.
# ---------------------------------------------------------------------------
def spherical_Delta_T(alpha, q, rt, b):
    # g = d/dr [r^2 - alpha^2 q^2] at r_t, so that sqrt(r^2 - alpha^2 q^2) ~ s sqrt(g)
    eps = 1e-6 * rt
    g = ((rt + eps) ** 2 - (alpha(rt + eps) * q) ** 2 - (rt - eps) ** 2 + (alpha(rt - eps) * q) ** 2) / (2 * eps)

    def den(s):
        r = rt + s * s
        arg = r * r - (alpha(r) * q) ** 2
        return np.sqrt(arg) if arg > 0 and s > 1e-6 * np.sqrt(rt) else s * np.sqrt(g)

    def fD(s):
        r = rt + s * s
        return 4 * s * alpha(r) * q / (r * den(s))

    def fT(s):
        r = rt + s * s
        return 4 * s * r / (alpha(r) * den(s))

    smax = np.sqrt(b - rt)
    D = quad(fD, 0, smax, limit=1000, epsabs=1e-12, epsrel=1e-11)[0]
    T = quad(fT, 0, smax, limit=1000, epsabs=1e-12, epsrel=1e-11)[0]
    return D, T


b = 6371.0

# ---------------------------------------------------------------------------
print("\n=== Example 1: homogeneous sphere ===")
alpha0 = 8.0
for q in [100.0, 400.0, 700.0]:
    rt = alpha0 * q
    D, T = spherical_Delta_T(lambda r: alpha0, q, rt, b)
    Dc = 2 * np.arccos(alpha0 * q / b)
    Tc = 2 * np.sqrt(b ** 2 - alpha0 ** 2 * q ** 2) / alpha0
    check(f"Delta closed form, q={q}", D, Dc, 1e-8)
    check(f"T closed form, q={q}", T, Tc, 1e-6)
    check(f"chord formula T=2b sin(Delta/2)/alpha, q={q}", Tc, 2 * b * np.sin(Dc / 2) / alpha0, 1e-9)
    # dT/dDelta = q by finite differences in q
    h = 1e-4
    dT = (2 * np.sqrt(b ** 2 - alpha0 ** 2 * (q + h) ** 2) / alpha0
          - 2 * np.sqrt(b ** 2 - alpha0 ** 2 * (q - h) ** 2) / alpha0) / (2 * h)
    dD = (2 * np.arccos(alpha0 * (q + h) / b) - 2 * np.arccos(alpha0 * (q - h) / b)) / (2 * h)
    check(f"dT/dDelta = q, q={q}", dT / dD, q, 1e-6, rel=True)
    check(f"dDelta/dq closed form, q={q}", dD, -2 * alpha0 / np.sqrt(b ** 2 - alpha0 ** 2 * q ** 2), 1e-6, rel=True)

# ---------------------------------------------------------------------------
print("\n=== Example 2: Benndorf relation and Herglotz-Wiechert for a smooth model ===")
# Smooth test model with r/alpha increasing: alpha(r) = a0 + a1 (1 - r/b)  (increases with depth)
a0, a1 = 8.0, 6.0
alpha_s = lambda r: a0 + a1 * (1 - r / b)


def rt_of_q(q, alpha, b, rmin=1.0):
    return brentq(lambda r: r / alpha(r) - q, rmin, b, xtol=1e-13)


def tau_of_q(q, alpha, b):
    rt = rt_of_q(q, alpha, b)

    def f(s):
        r = rt + s * s
        return 2 * s * 2 / (alpha(r) * r) * np.sqrt(r * r - (alpha(r) * q) ** 2)
    return quad(f, 0, np.sqrt(b - rt), limit=400, epsabs=1e-12, epsrel=1e-12)[0]


for q in [200.0, 400.0, 600.0]:
    rt = rt_of_q(q, alpha_s, b)
    D, T = spherical_Delta_T(alpha_s, q, rt, b)
    check(f"tau = T - q Delta, q={q}", tau_of_q(q, alpha_s, b), T - q * D, 1e-7, rel=True)
    h = 1e-3
    dtau = (tau_of_q(q + h, alpha_s, b) - tau_of_q(q - h, alpha_s, b)) / (2 * h)
    check(f"dtau/dq = -Delta, q={q}", dtau, -D, 1e-6, rel=True)
    Dp, Tp = spherical_Delta_T(alpha_s, q + h, rt_of_q(q + h, alpha_s, b), b)
    Dm, Tm = spherical_Delta_T(alpha_s, q - h, rt_of_q(q - h, alpha_s, b), b)
    check(f"dT/dDelta = q (smooth model), q={q}", (Tp - Tm) / (Dp - Dm), q, 1e-6, rel=True)

# Herglotz-Wiechert for the homogeneous sphere: q(Delta) = (b/alpha) cos(Delta/2)
for D1 in [np.radians(40.0), np.radians(100.0)]:
    q1 = (b / alpha0) * np.cos(D1 / 2)
    integrand = lambda D: np.arccosh((b / alpha0) * np.cos(D / 2) / q1)
    lhs = quad(integrand, 0, D1, limit=200)[0] / np.pi
    check(f"Herglotz-Wiechert, homogeneous sphere, Delta1={np.degrees(D1):.0f} deg",
          lhs, -np.log(np.cos(D1 / 2)), 1e-8)
# and for the smooth model, via numerically tabulated q(Delta)
q1 = 400.0
r1 = rt_of_q(q1, alpha_s, b)
D1 = spherical_Delta_T(alpha_s, q1, r1, b)[0]
qmax = b / alpha_s(b)


def q_of_D(D):
    # ray parameter of the ray emerging at epicentral angle D (Delta is monotonic in q)
    f = lambda q: spherical_Delta_T(alpha_s, q, rt_of_q(q, alpha_s, b), b)[0] - D
    return brentq(f, q1, qmax * (1 - 1e-9), xtol=1e-11)


# integrate over D with the substitution D = D1 - v^2 to keep the integrand smooth at D1
hw = quad(lambda v: 2 * v * np.arccosh(max(q_of_D(D1 - v * v) / q1, 1.0)), 0, np.sqrt(D1),
          limit=200, epsabs=1e-9, epsrel=1e-9)[0] / np.pi
check("Herglotz-Wiechert, smooth model, ln(b/r1)", hw, np.log(b / r1), 1e-6)

# ---------------------------------------------------------------------------
print("\n=== Example 3: flat-earth transformation and power-law model ===")


def flat_X_T(alpha_f, qf, zt):
    eps = 1e-6 * zt
    g = -((qf * alpha_f(zt + eps)) ** 2 - (qf * alpha_f(zt - eps)) ** 2) / (2 * eps)

    def den(s):
        z = zt - s * s
        arg = 1 - (qf * alpha_f(z)) ** 2
        return np.sqrt(arg) if arg > 0 and s > 1e-6 * np.sqrt(zt) else s * np.sqrt(g)

    def fX(s):
        z = zt - s * s
        return 4 * s * qf * alpha_f(z) / den(s)

    def fT(s):
        z = zt - s * s
        return 4 * s / (alpha_f(z) * den(s))
    X = quad(fX, 0, np.sqrt(zt), limit=1000, epsabs=1e-12, epsrel=1e-11)[0]
    T = quad(fT, 0, np.sqrt(zt), limit=1000, epsabs=1e-12, epsrel=1e-11)[0]
    return X, T


# (i) the transformation itself, for the smooth linear model above
for q in [200.0, 500.0]:
    rt = rt_of_q(q, alpha_s, b)
    D, T = spherical_Delta_T(alpha_s, q, rt, b)
    alpha_f = lambda z: (b / (b * np.exp(-z / b))) * alpha_s(b * np.exp(-z / b))
    zt = b * np.log(b / rt)
    X, Tf = flat_X_T(alpha_f, q / b, zt)
    check(f"flat-earth X = b Delta (linear model), q={q}", X, b * D, 1e-7, rel=True)
    check(f"flat-earth T = spherical T (linear model), q={q}", Tf, T, 1e-7, rel=True)
    check(f"turning depth consistent, q={q}", (q / b) * alpha_f(zt), 1.0, 1e-10)

# (ii) power-law model alpha = alpha0 (b/r)^nu
for nu in [0.0, 1.0, 0.6, -0.5]:
    alpha_p = lambda r, nu=nu: alpha0 * (b / r) ** nu
    for q in [300.0, 600.0]:
        rt_c = b * (alpha0 * q / b) ** (1 / (1 + nu))
        check(f"turning radius, nu={nu}, q={q}", rt_c / alpha_p(rt_c), q, 1e-9, rel=True)
        D, T = spherical_Delta_T(alpha_p, q, rt_c, b)
        Dc = 2 / (1 + nu) * np.arccos(alpha0 * q / b)
        Tc = 2 / ((1 + nu) * alpha0) * np.sqrt(b ** 2 - alpha0 ** 2 * q ** 2)
        check(f"power law Delta(q), nu={nu}, q={q}", D, Dc, 1e-8)
        check(f"power law T(q), nu={nu}, q={q}", T, Tc, 1e-6)
        check(f"power law T(Delta), nu={nu}, q={q}", 2 * b / ((1 + nu) * alpha0) * np.sin((1 + nu) * Dc / 2), Tc, 1e-9)
        # transformed model is exponential with h = b/(1+nu)
        z = 1234.5
        check(f"alpha_f exponential, nu={nu}", (b / (b * np.exp(-z / b))) * alpha_p(b * np.exp(-z / b)),
              alpha0 * np.exp((1 + nu) * z / b), 1e-9, rel=True)
    # Benndorf for the power law
    h = 1e-4
    q = 300.0
    Tq = lambda q: 2 / ((1 + nu) * alpha0) * np.sqrt(b ** 2 - alpha0 ** 2 * q ** 2)
    Dq = lambda q: 2 / (1 + nu) * np.arccos(alpha0 * q / b)
    check(f"dT/dDelta = q, power law nu={nu}", (Tq(q + h) - Tq(q - h)) / (Dq(q + h) - Dq(q - h)), q, 1e-6, rel=True)

# numbers quoted in the text for nu = 1, q = 300
nu, q = 1.0, 300.0
Dq = np.degrees(2 / (1 + nu) * np.arccos(alpha0 * q / b))
Tq = 2 / ((1 + nu) * alpha0) * np.sqrt(b ** 2 - alpha0 ** 2 * q ** 2)
rt = b * (alpha0 * q / b) ** (1 / (1 + nu))
print(f"     nu=1, q=300: Delta = {Dq:.2f} deg, T = {Tq:.1f} s, r_t = {rt:.0f} km")
check("quoted Delta = 67.87 deg", round(Dq, 2), 67.87, 1e-9)
check("quoted T = 737.7 s", round(Tq, 1), 737.7, 1e-9)
check("quoted r_t = 3910 km", round(rt), 3910, 1e-9)
check("quoted alpha at r = 3480 km, nu = 1", round(alpha0 * b / 3480.0, 1), 14.6, 1e-9)
# antiderivative used in eq. (19)
u = 0.37
F = lambda u: -np.sqrt(1 - u * u) / u
check("antiderivative of 1/(u^2 sqrt(1-u^2))", (F(u + 1e-6) - F(u - 1e-6)) / 2e-6, 1 / (u * u * np.sqrt(1 - u * u)), 1e-6, rel=True)

# ---------------------------------------------------------------------------
print("\n=== Example 4: hypocentre location ===")
alpha = 6.0
X = np.array([[30., 0., 0.], [0., 40., 0.], [-25., 15., 0.], [-10., -35., 0.], [20., -20., 0.]])
xtrue = np.array([3., -2., 12.])
sig = 0.1 * np.ones(len(X))
Tfun = lambda xh: np.linalg.norm(X - xh, axis=1) / alpha
t = 0.0 + Tfun(xtrue)
print("     synthetic arrival times:", np.round(t, 3))
for tv, te in zip(np.round(t, 3), [4.936, 7.297, 5.814, 6.241, 4.586]):
    check("quoted arrival time", tv, te, 1e-9)

# derivative check by finite differences
xh0 = np.array([1., 2., 7.])
d = X - xh0
D = np.linalg.norm(d, axis=1)
grad_analytic = -d / (alpha * D[:, None])
eps = 1e-6
for j in range(3):
    e = np.zeros(3)
    e[j] = eps
    fd = (Tfun(xh0 + e) - Tfun(xh0 - e)) / (2 * eps)
    check(f"dT/dx_h component {j+1} (max abs err over stations)",
          np.max(np.abs(fd - grad_analytic[:, j])), 0.0, 1e-7)

m = np.array([0., 0., 5., 1.])
history = []
for k in range(6):
    xh, th = m[:3], m[3]
    d = X - xh
    D = np.linalg.norm(d, axis=1)
    r = (th + D / alpha - t) / sig
    J = np.sum(r ** 2)
    history.append((m.copy(), J))
    G = np.column_stack([-d / (alpha * D[:, None]), np.ones(len(X))]) / sig[:, None]
    dm = np.linalg.solve(G.T @ G, -G.T @ r)
    m = m + dm
for k, (mk, Jk) in enumerate(history[:4]):
    print(f"     iteration {k}: x_h = {np.round(mk[:3], 3)}, t_h = {mk[3]:.4f}, J = {Jk:.4g}")
table = [([0.000, 0.000, 5.000], 1.0000), ([2.795, -1.868, 13.377], 0.0966),
         ([3.004, -2.001, 12.034], 0.0017), ([3.000, -2.000, 12.000], 0.0000)]
for k, ((xq, tq), (mk, Jk)) in enumerate(zip(table, history)):
    check(f"table row {k} hypocentre", np.max(np.abs(np.round(mk[:3], 3) - np.array(xq))), 0.0, 1e-9)
    check(f"table row {k} origin time", round(mk[3], 4), tq, 1e-9)
Js = [Jk for _, Jk in history]
check("J iteration 0 = 370", round(Js[0]), 370, 1e-9)
check("J iteration 1 = 18", round(Js[1]), 18, 1e-9)
check("J iteration 2 = 6.7e-3", Js[2], 6.7e-3, 0.05e-3)
check("J iteration 3 = 3.4e-9", Js[3], 3.4e-9, 0.05e-9)
check("converged hypocentre", np.max(np.abs(history[5][0][:3] - xtrue)), 0.0, 1e-8)
check("converged origin time", history[5][0][3], 0.0, 1e-8)

# symmetric configuration
hdep, R = 12.0, 30.0
az = np.radians([0, 72, 144, 216, 288])
Xs = np.column_stack([R * np.cos(az), R * np.sin(az), np.zeros(5)])
xh = np.array([0., 0., hdep])
d = Xs - xh
D = np.linalg.norm(d, axis=1)
G = np.column_stack([-d / (alpha * D[:, None]), np.ones(5)]) / 0.1
w, V = np.linalg.eigh(G.T @ G)
print("     eigenvalues of normal matrix:", w)
check("smallest eigenvalue is zero", w[0], 0.0, 1e-10)
check("eigenvalue 5.99 (x2)", round(w[1], 2), 5.99, 1e-9)
check("eigenvalue 5.99 (x2)", round(w[2], 2), 5.99, 1e-9)
check("eigenvalue 502", round(w[3]), 502, 1e-9)
v = np.array([0, 0, 1, -hdep / (alpha * D[0])])
check("G v = 0 for v = (0,0,1,-h/(alpha D))", np.linalg.norm(G @ v), 0.0, 1e-12)
check("h/(alpha D) = 0.0619", round(hdep / (alpha * D[0]), 4), 0.0619, 1e-9)
check("D = 32.31 km", round(D[0], 2), 32.31, 1e-9)

# ---------------------------------------------------------------------------
print("\n=== Example 5: P-wave shadow zone ===")
c, am, ac = 3480.0, 13.0, 8.0
q1 = c / am
D1 = 2 * np.arccos(c / b)
T1 = 2 * np.sqrt(b ** 2 - c ** 2) / am
print(f"     q1 = {q1:.1f} s/rad, Delta_1 = {np.degrees(D1):.2f} deg, T1 = {T1:.1f} s, "
      f"theta0 = {np.degrees(np.arcsin(c / b)):.1f} deg from vertical")
check("q1 = 267.7", round(q1, 1), 267.7, 1e-9)
check("Delta_1 = 113.8 deg", round(np.degrees(D1), 1), 113.8, 1e-9)
check("T1 = 821 s", round(T1), 821, 1e-9)
check("grazing take-off 33.1 deg", round(np.degrees(np.arcsin(c / b)), 1), 33.1, 1e-9)


def Dcore(q):
    return 2 * (np.arccos(am * q / b) - np.arccos(am * q / c)) + 2 * np.arccos(ac * q / c)


def Tcore(q):
    return (2 * (np.sqrt(b ** 2 - (am * q) ** 2) - np.sqrt(c ** 2 - (am * q) ** 2)) / am
            + 2 * np.sqrt(c ** 2 - (ac * q) ** 2) / ac)


# quadrature check of the assembled Delta and T for a core ray
for q in [100.0, 200.0, 260.0]:
    Dm = quad(lambda r: 2 * am * q / (r * np.sqrt(r * r - (am * q) ** 2)), c, b, epsabs=1e-12, epsrel=1e-12)[0]
    Tm = quad(lambda r: 2 * r / (am * np.sqrt(r * r - (am * q) ** 2)), c, b, epsabs=1e-12, epsrel=1e-12)[0]
    Dc_, Tc_ = spherical_Delta_T(lambda r: ac, q, ac * q, c)
    check(f"core-ray Delta(q) vs quadrature, q={q}", Dm + Dc_, Dcore(q), 1e-8)
    check(f"core-ray T(q) vs quadrature, q={q}", Tm + Tc_, Tcore(q), 1e-6)
    h = 1e-4
    check(f"Benndorf for core ray, q={q}", (Tcore(q + h) - Tcore(q - h)) / (Dcore(q + h) - Dcore(q - h)), q, 1e-6, rel=True)

Dgraze = Dcore(q1 * (1 - 1e-12))
check("Delta(q1-) = Delta_1 + 2 arccos(alpha_c/alpha_m)", Dgraze, D1 + 2 * np.arccos(ac / am), 1e-5)
check("2 arccos(8/13) = 104.0 deg", round(np.degrees(2 * np.arccos(ac / am)), 1), 104.0, 1e-9)
check("Delta(q1-) = 217.8 deg", round(np.degrees(Dgraze), 1), 217.8, 1e-9)
check("Delta(0) = 180 deg", np.degrees(Dcore(0.0)), 180.0, 1e-9)
check("theta_c = 38.0 deg", round(np.degrees(np.arcsin(ac / am)), 1), 38.0, 1e-9)
res = minimize_scalar(Dcore, bounds=(0, q1), method="bounded", options={"xatol": 1e-12})
print(f"     caustic: q_B = {res.x:.1f} s/rad, Delta_B = {np.degrees(res.fun):.2f} deg")
check("q_B = 157.6", round(res.x, 1), 157.6, 1e-9)
check("Delta_B = 172.1 deg", round(np.degrees(res.fun), 1), 172.1, 1e-9)
qs = np.linspace(0, q1, 400001)
dd = np.minimum(Dcore(qs), 2 * np.pi - Dcore(qs))
print(f"     min angular distance of core rays = {np.degrees(dd.min()):.2f} deg at q = {qs[dd.argmin()]:.1f}")
check("delta_min attained at q -> q1", qs[dd.argmin()], q1, 1e-9)
check("delta_min = 142.2 deg", round(np.degrees(dd.min()), 1), 142.2, 1e-9)
check("delta_min = 360 - 113.8 - 104.0", 360 - np.degrees(D1) - np.degrees(2 * np.arccos(ac / am)), np.degrees(dd.min()), 1e-6)
# arccos(c/b) = 56.9 deg
check("arccos(c/b) = 56.9 deg", round(np.degrees(np.arccos(c / b)), 1), 56.9, 1e-9)
# branch structure: Delta = 180 crossed once between q_B and q1
qa = brentq(lambda q: Dcore(q) - np.pi, res.x, q1)
print(f"     Delta = 180 deg at q_a = {qa:.1f} s/rad (three arrivals for 172.1 < delta < 180)")

print("\nALL CHECKS PASSED" if failures == 0 else f"\n{failures} CHECK(S) FAILED")
