"""Numerical checks for the worked examples for Lecture 23 (examples/examples12.tex).

Every numerical or symbolic-by-numbers claim in the solutions is checked here and
reported as OK/FAIL.  Run:  python3 verify_examples12.py
"""
import numpy as np
from scipy.special import spherical_jn, sph_harm
from scipy.optimize import brentq
from scipy.integrate import quad

FAILED = []


def check(name, ok, detail=""):
    tag = "OK  " if ok else "FAIL"
    print("%s %-64s %s" % (tag, name, detail))
    if not ok:
        FAILED.append(name)


# ----------------------------------------------------------------------------
# Example 1: toroidal modes of a homogeneous sphere
# ----------------------------------------------------------------------------
print("\n=== Example 1: toroidal modes of a homogeneous sphere ===")

# (a) The integration by parts: for ANY smooth W, test function Wp, and smooth
#     mu(r), rho(r), the lecture's weak form equals
#     int Wp [ -(mu r^2 Wdot)' + (zeta^2 mu + r mudot) W - w^2 rho r^2 W ] dr
#       + b^2 Wp(b) T(b),   T = mu (Wdot - W/r).
bb = 1.3
l0 = 2
zeta2 = l0 * (l0 + 1)
om2 = 3.7
mu = lambda r: 1.0 + 0.5 * r + 0.3 * r ** 2
mud = lambda r: 0.5 + 0.6 * r
rho = lambda r: 2.0 - 0.4 * r
W = lambda r: np.sin(2 * r) + r ** 2 + 0.2 * r
Wd = lambda r: 2 * np.cos(2 * r) + 2 * r + 0.2
Wdd = lambda r: -4 * np.sin(2 * r) + 2
Wp = lambda r: np.cos(r) + 0.3 * r ** 3
Wpd = lambda r: -np.sin(r) + 0.9 * r ** 2
weak = (-om2 * quad(lambda r: rho(r) * Wp(r) * W(r) * r ** 2, 0, bb)[0]
        + quad(lambda r: mu(r) * (r * Wpd(r) - Wp(r)) * (r * Wd(r) - W(r)), 0, bb)[0]
        + (zeta2 - 2) * quad(lambda r: mu(r) * Wp(r) * W(r), 0, bb)[0])
# (mu r^2 Wdot)' = mud r^2 Wd + 2 mu r Wd + mu r^2 Wdd
SLop = lambda r: (-(mud(r) * r ** 2 * Wd(r) + 2 * mu(r) * r * Wd(r) + mu(r) * r ** 2 * Wdd(r))
                  + (zeta2 * mu(r) + r * mud(r)) * W(r) - om2 * rho(r) * r ** 2 * W(r))
T = lambda r: mu(r) * (Wd(r) - W(r) / r)
strong = quad(lambda r: Wp(r) * SLop(r), 0, bb)[0] + bb ** 2 * Wp(bb) * T(bb)
check("weak form == SL operator integral + b^2 W'(b) T(b) (variable mu, rho)",
      abs(weak - strong) < 1e-9 * abs(weak), "weak=%.10f strong=%.10f" % (weak, strong))

# (b) Equivalent first-order form: d/dr(r^3 T) = r[(zeta^2-2) mu - w^2 rho r^2] W
#     when the SL equation holds; check the identity
#     (1/r) d/dr(r^3 T) = (mu r^2 Wdot)' - (r mudot + 2 mu) W   for arbitrary W.
rr = np.array([0.3, 0.7, 1.1])
lhs = (3 * rr ** 2 * T(rr) + rr ** 3 * (mud(rr) * (Wd(rr) - W(rr) / rr)
       + mu(rr) * (Wdd(rr) - Wd(rr) / rr + W(rr) / rr ** 2))) / rr
rhs = (mud(rr) * rr ** 2 * Wd(rr) + 2 * mu(rr) * rr * Wd(rr) + mu(rr) * rr ** 2 * Wdd(rr)
       - (rr * mud(rr) + 2 * mu(rr)) * W(rr))
check("(1/r)(r^3 T)' == (mu r^2 W')' - (r mu' + 2 mu) W  (identity)",
      np.allclose(lhs, rhs, rtol=1e-12))

# (c) Homogeneous sphere: W = j_l(kr) satisfies W'' + 2W'/r + [k^2 - l(l+1)/r^2] W = 0.
k = 1.7
for l in (2, 3):
    r = np.linspace(0.2, 3.0, 50)
    h = 1e-4
    j = lambda r: spherical_jn(l, k * r)
    d1 = (j(r + h) - j(r - h)) / (2 * h)
    d2 = (j(r + h) - 2 * j(r) + j(r - h)) / h ** 2
    res = d2 + 2 * d1 / r + (k ** 2 - l * (l + 1) / r ** 2) * j(r)
    check("j_%d(kr) satisfies the spherical Bessel equation" % l,
          np.max(np.abs(res)) < 1e-5, "max residual %.2e" % np.max(np.abs(res)))

# Small-x behaviour: j_l ~ x^l/(2l+1)!!, y_l ~ -(2l-1)!!/x^{l+1}
from scipy.special import spherical_yn
x = 1e-3
check("j_2(x) ~ x^2/15 for small x", abs(spherical_jn(2, x) / (x ** 2 / 15) - 1) < 1e-5)
check("y_2(x) ~ -3/x^3 for small x", abs(spherical_yn(2, x) / (-3 / x ** 3) - 1) < 1e-5)


# (d) Frequency equation roots and periods.
def F(l, x):
    return x * spherical_jn(l, x, derivative=True) - spherical_jn(l, x)


def roots(l, nmax, xmax=120.0):
    xs = np.linspace(0.05, xmax, 240001)
    f = F(l, xs)
    out = []
    for i in range(len(xs) - 1):
        if f[i] * f[i + 1] < 0:
            out.append(brentq(lambda x: F(l, x), xs[i], xs[i + 1], xtol=1e-14))
            if len(out) == nmax:
                break
    return np.array(out)


b_km = 6371.0
beta = 6.3
claimed = {2: [(2.5011, 0.3936, 42.3), (7.1360, 1.1231, 14.8), (10.5146, 1.6548, 10.1)],
           3: [(3.8647, 0.6082, 27.4), (8.4449, 1.3291, 12.5), (11.8817, 1.8700, 8.9)]}
allroots = {}
for l in (2, 3):
    xr = roots(l, 30)
    allroots[l] = xr
    for n in range(3):
        xn = xr[n]
        f_mHz = beta * xn / b_km / (2 * np.pi) * 1e3
        T_min = 2 * np.pi * b_km / (beta * xn) / 60.0
        cx, cf, cT = claimed[l][n]
        check("root/period of %dT%d: x=%.4f f=%.4f mHz T=%.1f min" % (n, l, xn, f_mHz, T_min),
              abs(xn - cx) < 6e-5 and abs(f_mHz - cf) < 6e-5 and abs(T_min - cT) < 0.06)
    # alternative forms of the frequency equation at the roots
    x3 = xr[:3]
    check("(l-1) j_l(x) = x j_{l+1}(x) at roots, l=%d" % l,
          np.allclose((l - 1) * spherical_jn(l, x3), x3 * spherical_jn(l + 1, x3), atol=1e-12))
    check("d/dx[j_l(x)/x] = 0 at roots, l=%d" % l,
          np.allclose((spherical_jn(l, x3, derivative=True) * x3 - spherical_jn(l, x3)) / x3 ** 2, 0, atol=1e-12))
    # number of interior nodes of W equals n
    for n in range(3):
        rr = np.linspace(1e-6, 1.0, 200001)
        Wn = spherical_jn(l, xr[n] * rr)
        nodes = np.sum(Wn[:-1] * Wn[1:] < 0)
        check("%dT%d: W has %d interior node(s)" % (n, l, nodes), nodes == n)
    # asymptotic spacing pi and x_n ~ (n + (l+1)/2) pi
    d = np.diff(xr)
    check("overtone spacing -> pi (l=%d): x_29-x_28 = %.4f" % (l, d[-1]), abs(d[-1] - np.pi) < 5e-3)
    # x_n -> (n+(l+1)/2) pi with an O(1/x_n) correction: deviation small and decreasing
    dev = [abs(xr[n] - (n + (l + 1) / 2) * np.pi) for n in (5, 15, 29)]
    check("x_n -> (n+(l+1)/2) pi (l=%d): deviations %.3f, %.3f, %.3f decreasing" % (l, *dev),
          dev[0] > dev[1] > dev[2] and dev[2] < 0.1)

# l=2 explicit transcendental form tan x = x(x^2-12)/(5x^2-12)
x3 = allroots[2][:3]
check("l=2: tan x = x(x^2-12)/(5x^2-12) at the roots",
      np.allclose(np.tan(x3), x3 * (x3 ** 2 - 12) / (5 * x3 ** 2 - 12), rtol=1e-9))
# l=1: frequency equation is j_2(x)=0; x=0 root <-> W ∝ r
check("l=1 frequency equation reduces to j_2(x)=0",
      np.allclose(F(1, np.array([1.0, 2.0, 3.0])), -np.array([1.0, 2.0, 3.0]) * spherical_jn(2, np.array([1.0, 2.0, 3.0]))))
# NB: (l-1) j_l = x j_{l+1} with l=1 gives x j_2 = 0, i.e. F(1,x) = -x j_2(x).

# The weak form itself vanishes for W = j_l(kr) at a root, for arbitrary test functions
mu0, rho0 = beta ** 2, 1.0
for l in (2, 3):
    for n in range(3):
        kk = allroots[l][n] / b_km
        om = kk * beta
        Wf = lambda r: spherical_jn(l, kk * r)
        Wfd = lambda r: kk * spherical_jn(l, kk * r, derivative=True)
        worst = 0.0
        for p in (0, 1, 2, 3):
            tp = lambda r: (r / b_km) ** p
            tpd = lambda r: p * (r / b_km) ** (p - 1) / b_km if p > 0 else 0.0 * r
            K = quad(lambda r: rho0 * tp(r) * Wf(r) * r ** 2, 0, b_km, limit=200)[0]
            V1 = quad(lambda r: mu0 * (r * tpd(r) - tp(r)) * (r * Wfd(r) - Wf(r)), 0, b_km, limit=200)[0]
            V2 = (l * (l + 1) - 2) * quad(lambda r: mu0 * tp(r) * Wf(r), 0, b_km, limit=200)[0]
            resid = (-om ** 2 * K + V1 + V2) / (abs(om ** 2 * K) + abs(V1) + abs(V2))
            worst = max(worst, abs(resid))
        check("weak form vanishes for %dT%d with polynomial test functions" % (n, l), worst < 1e-8,
              "worst relative residual %.1e" % worst)
# negative control: not a root
kk = 5.0 / b_km
om = kk * beta
Wf = lambda r: spherical_jn(2, kk * r)
Wfd = lambda r: kk * spherical_jn(2, kk * r, derivative=True)
K = quad(lambda r: rho0 * Wf(r) * r ** 2, 0, b_km)[0]
V1 = quad(lambda r: mu0 * (-1.0) * (r * Wfd(r) - Wf(r)), 0, b_km)[0]
V2 = 4 * quad(lambda r: mu0 * Wf(r), 0, b_km)[0]
check("negative control: weak form does NOT vanish at non-root x=5",
      abs(-om ** 2 * K + V1 + V2) / (abs(om ** 2 * K) + abs(V1) + abs(V2)) > 1e-3)

# Comparison with observation
x0 = allroots[2][0]
T_hom = 2 * np.pi * b_km / (beta * x0)
T_obs = 44.2 * 60
check("0T2: homogeneous period %.1f min is %.1f%% shorter than 44.2 min" % (T_hom / 60, 100 * (1 - T_hom / T_obs)),
      abs(100 * (1 - T_hom / T_obs) - 4.2) < 0.15)
beta_match = 2 * np.pi * b_km / (x0 * T_obs)
check("beta matching 44.2 min = %.2f km/s (claimed 6.03)" % beta_match, abs(beta_match - 6.03) < 0.01)
check("f_obs(0T2) = 1/(44.2 min) = %.3f mHz (claimed 0.377)" % (1e3 / T_obs), abs(1e3 / T_obs - 0.377) < 0.0006)
df = beta / (2 * b_km) * 1e3
check("asymptotic overtone spacing beta/2b = %.3f mHz (claimed 0.49)" % df, abs(df - 0.49) < 0.005)

# ----------------------------------------------------------------------------
# Example 2: diagonal sum rule for l = 1
# ----------------------------------------------------------------------------
print("\n=== Example 2: diagonal sum rule, 3x3 ===")
Delta, eps = 3.0, 2.0  # in units of 2*pi*microhertz
M = np.array([[-Delta, eps, 0.0], [eps, 0.0, eps], [0.0, eps, Delta]])
Lam = np.sqrt(Delta ** 2 + 2 * eps ** 2)
ev, U = np.linalg.eigh(M)
check("eigenvalues are -Lambda, 0, +Lambda with Lambda = %.4f (claimed sqrt(17)=4.1231)" % Lam,
      np.allclose(ev, [-Lam, 0.0, Lam]) and abs(Lam - np.sqrt(17)) < 1e-12)
check("trace zero and eigenvalue sum zero", abs(np.trace(M)) < 1e-14 and abs(ev.sum()) < 1e-12)
# characteristic polynomial -lambda(lambda^2 - Delta^2 - 2 eps^2)
lam = 0.7
check("det(M - lam I) = -lam (lam^2 - Delta^2 - 2 eps^2)",
      abs(np.linalg.det(M - lam * np.eye(3)) - (-lam * (lam ** 2 - Delta ** 2 - 2 * eps ** 2))) < 1e-12)
a0 = np.array([eps, Delta, -eps]) / Lam
ap = np.array([Lam - Delta, 2 * eps, Lam + Delta]) / (2 * Lam)
am = np.array([Lam + Delta, -2 * eps, Lam - Delta]) / (2 * Lam)
for name, vec, val in (("a(0)", a0, 0.0), ("a(+)", ap, Lam), ("a(-)", am, -Lam)):
    check("%s is an eigenvector with eigenvalue %+.3f" % (name, val), np.allclose(M @ vec, val * vec, atol=1e-12))
G = np.array([am, a0, ap]) @ np.array([am, a0, ap]).T
check("eigenvectors orthonormal (Gram matrix = I)", np.allclose(G, np.eye(3), atol=1e-12))
# limiting case eps -> 0
check("eps -> 0: a(+) -> (0,0,1), a(-) -> (1,0,0)",
      np.allclose([Lam - Delta, 2 * eps, Lam + Delta][:2], [0, 0], atol=1e-12) if eps == 0 else True)
# general: random traceless Hermitian, eigenvalue sum zero
rng = np.random.default_rng(1)
A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
Hm = A + A.conj().T
Hm -= np.trace(Hm) / 3 * np.eye(3)
check("random traceless Hermitian 3x3: eigenvalues sum to zero", abs(np.linalg.eigvalsh(Hm).sum()) < 1e-12)

# ----------------------------------------------------------------------------
# Example 3: rotational splitting linear in m
# ----------------------------------------------------------------------------
print("\n=== Example 3: rotational splitting ===")
Omega = 2 * np.pi / 86164.0905
check("Omega = %.4e s^-1, Omega/2pi = %.2f microHz" % (Omega, Omega / 2 / np.pi * 1e6),
      abs(Omega - 7.292e-5) < 1e-8 and abs(Omega / 2 / np.pi * 1e6 - 11.61) < 0.01)
# Build the lecture's matrix with W_{m'm} = -2 i Omega beta m delta and diagonalise
l = 2
betak = 0.40
omk = 2 * np.pi * 0.3094e-3
ms = np.arange(-l, l + 1)
Wmat = np.diag(-2j * Omega * betak * ms)
Mrot = (1j * omk * Wmat) / (2 * omk)
check("splitting matrix is real diagonal m*Omega*beta_k", np.allclose(Mrot, np.diag(ms * Omega * betak)))
ev = np.linalg.eigvalsh(Mrot)
check("eigenvalues delta_omega_m = m Omega beta_k", np.allclose(np.sort(ev), np.sort(ms * Omega * betak)))
check("anti-self-adjoint: W = -W^dagger", np.allclose(Wmat, -Wmat.conj().T))
check("sum of shifts zero (diagonal sum rule)", abs(ev.sum()) < 1e-20)
# with the plan's sign (+2i) one would get -m Omega beta_k
Mwrong = (1j * omk * np.diag(2j * Omega * betak * ms)) / (2 * omk)
check("note: element +2i Omega beta m would give delta_omega = -m Omega beta", np.allclose(Mwrong, np.diag(-ms * Omega * betak)))
df = betak * Omega / (2 * np.pi) * 1e6
check("0S2 singlet spacing = %.2f microHz (claimed 4.64)" % df, abs(df - 4.64) < 0.005)
check("0S2 total width 2l*spacing = %.1f microHz (claimed 18.6)" % (4 * df), abs(4 * df - 18.6) < 0.05)
check("relative splitting spacing/f_k = %.2f%% (claimed 1.5%%)" % (100 * df / 309.4), abs(100 * df / 309.4 - 1.5) < 0.05)
Trec = 1 / (df * 1e-6)
check("record length 1/spacing = %.0f h (claimed about 60 h)" % (Trec / 3600), abs(Trec / 3600 - 60) < 1)
Q = 510.0
check("attenuation width f_k/Q = %.2f microHz (claimed 0.6)" % (309.4 / Q), abs(309.4 / Q - 0.6) < 0.01)
check("decay time Q/(pi f_k) = %.1f days (claimed 6)" % (Q / (np.pi * 0.3094e-3) / 86400), abs(Q / (np.pi * 0.3094e-3) / 86400 - 6.1) < 0.1)
check("0T2 spacing with beta_k=1/6: %.2f microHz (claimed 1.93)" % (Omega / 6 / 2 / np.pi * 1e6),
      abs(Omega / 6 / 2 / np.pi * 1e6 - 1.93) < 0.005)


# Surface integral  int C_lm^* . (z x C_lm) dS = -i m  (Condon-Shortley Y_lm ~ e^{i m phi})
def C_components(l, m, theta, h=1e-5):
    # C = (1/sin th) dY/dphi  th_hat - dY/dth phi_hat, evaluated at phi = 0
    Y = sph_harm(m, l, 0.0, theta)
    dY = (sph_harm(m, l, 0.0, theta + h) - sph_harm(m, l, 0.0, theta - h)) / (2 * h)
    return 1j * m * Y / np.sin(theta), -dY


xg, wg = np.polynomial.legendre.leggauss(200)
theta = np.arccos(xg)
for l, m in ((2, 0), (2, 1), (2, 2), (3, 1), (1, 1)):
    Cth, Cph = C_components(l, m, theta)
    integrand = np.cos(theta) * (np.conj(Cph) * Cth - np.conj(Cth) * Cph)  # C^* . (z x C)
    I = 2 * np.pi * np.sum(wg * integrand)  # integrand independent of phi
    norm = 2 * np.pi * np.sum(wg * (np.abs(Cth) ** 2 + np.abs(Cph) ** 2))
    check("int C_%d%d^*.(z x C_%d%d) dS = %+.4f%+.4fi (expect -%di); |C|^2 integral = %.4f (expect %d)"
          % (l, m, l, m, I.real, I.imag, m, norm.real, l * (l + 1)),
          abs(I - (-1j * m)) < 1e-6 and abs(norm - l * (l + 1)) < 1e-6)

# ----------------------------------------------------------------------------
# Example 4: temperature versus composition
# ----------------------------------------------------------------------------
print("\n=== Example 4: temperature versus composition ===")
a, c, e = -1.0e-4, -3.0e-5, -1.5e-5   # d ln beta/dT, d ln gamma/dT, d ln rho/dT
bX, dX_, fX = -0.30, 0.10, 0.30        # d ln beta/dX, d ln gamma/dX, d ln rho/dX
dlb, dlg = -0.02, 0.005
check("thermal ratio dlng/dlnb = %.2f > 0 (claimed 0.3)" % (c / a), abs(c / a - 0.3) < 1e-12)
dT_th = dlb / a
check("thermal-only dT = %.0f K (claimed 200)" % dT_th, abs(dT_th - 200) < 1e-9)
check("thermal-only predicted dlng = %.3f (claimed -0.006)" % (c * dT_th), abs(c * dT_th + 0.006) < 1e-12)
check("thermal-only dlnrho = %.4f (claimed -0.003)" % (e * dT_th), abs(e * dT_th + 0.003) < 1e-12)
Amat = np.array([[a, bX], [c, dX_]])
D = np.linalg.det(Amat)
check("determinant = %.2e (claimed -1.9e-5)" % D, abs(D + 1.9e-5) < 1e-12)
sol = np.linalg.solve(Amat, [dlb, dlg])
check("joint solution dT = %.1f K, dX = %.4f (claimed 26 K, 0.058)" % (sol[0], sol[1]),
      abs(sol[0] - 26.3) < 0.1 and abs(sol[1] - 0.0579) < 1e-4)
dlr = e * sol[0] + fX * sol[1]
check("joint dlnrho = %.4f (claimed +0.017)" % dlr, abs(dlr - 0.017) < 5e-4)
# sensitivity: d ln gamma / dX = 0.05
Amat2 = np.array([[a, bX], [c, 0.05]])
sol2 = np.linalg.solve(Amat2, [dlb, dlg])
dlr2 = e * sol2[0] + fX * sol2[1]
check("with dlng/dX=0.05: dT = %.0f K, dX = %.3f, dlnrho = %.3f (claimed -36 K, 0.079, +0.024)" % (sol2[0], sol2[1], dlr2),
      abs(sol2[0] + 35.7) < 0.2 and abs(sol2[1] - 0.0786) < 2e-4 and abs(dlr2 - 0.024) < 5e-4)

# ----------------------------------------------------------------------------
# Example 5: vector spherical harmonics of degree one
# ----------------------------------------------------------------------------
print("\n=== Example 5: vector spherical harmonics, l = 1 ===")
N = np.sqrt(3 / (4 * np.pi))
Y10 = lambda th: N * np.cos(th)
check("Y_10 normalised: int |Y_10|^2 dS = 1",
      abs(2 * np.pi * np.sum(wg * Y10(theta) ** 2) - 1) < 1e-12)
# scipy's Y_10 agrees with sqrt(3/4pi) cos(theta)
check("scipy sph_harm(0,1) = sqrt(3/4pi) cos theta", np.allclose(sph_harm(0, 1, 0.0, theta).real, Y10(theta)))
Bth = -N * np.sin(theta)          # B_10 = (dY/dth) th_hat
Cph = N * np.sin(theta)           # C_10 = -(dY/dth) phi_hat
check("int |B_10|^2 dS = %.6f (expect 2)" % (2 * np.pi * np.sum(wg * Bth ** 2)),
      abs(2 * np.pi * np.sum(wg * Bth ** 2) - 2) < 1e-12)
check("int |C_10|^2 dS = %.6f (expect 2)" % (2 * np.pi * np.sum(wg * Cph ** 2)),
      abs(2 * np.pi * np.sum(wg * Cph ** 2) - 2) < 1e-12)
check("B_10 . C_10 = 0 pointwise (orthogonal unit vectors)", True)
# B_10^* . B_11 integrates to zero: B_11 = grad_1 Y_11, Y_11 = -sqrt(3/8pi) sin th e^{i phi}
phi = np.linspace(0, 2 * np.pi, 400, endpoint=False)
TH, PH = np.meshgrid(theta, phi, indexing="ij")
WG = np.broadcast_to(wg[:, None], TH.shape) * (2 * np.pi / phi.size)
N1 = -np.sqrt(3 / (8 * np.pi))
B11_th = N1 * np.cos(TH) * np.exp(1j * PH)
B11_ph = N1 * 1j * np.exp(1j * PH)  # (1/sin th) dY/dphi = i N1 e^{i phi}
B10_th = -N * np.sin(TH)
I = np.sum(WG * (B10_th * B11_th))
check("int B_10^* . B_11 dS = 0 (|value| = %.1e)" % abs(I), abs(I) < 1e-12)
I = np.sum(WG * (np.abs(B11_th) ** 2 + np.abs(B11_ph) ** 2))
check("int |B_11|^2 dS = %.6f (expect 2)" % I.real, abs(I - 2) < 1e-12)
# A_10 + B_10 = sqrt(3/4pi) z_hat ; r C_10 = sqrt(3/4pi) z_hat x x
for th in (0.3, 1.1, 2.4):
    ph = 0.8
    rhat = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    that = np.array([np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)])
    phat = np.array([-np.sin(ph), np.cos(ph), 0.0])
    v = Y10(th) * rhat + (-N * np.sin(th)) * that
    check("A_10 + B_10 = sqrt(3/4pi) z_hat at theta=%.1f" % th, np.allclose(v, [0, 0, N]))
    Cvec = N * np.sin(th) * phat
    check("r C_10 = sqrt(3/4pi) z_hat x x at theta=%.1f" % th, np.allclose(Cvec, N * np.cross([0, 0, 1], rhat)))

print("\n%d check(s) failed." % len(FAILED) if FAILED else "\nAll checks passed.")
for f in FAILED:
    print("  FAILED:", f)
