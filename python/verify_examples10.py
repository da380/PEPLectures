"""Numerical checks for the worked examples accompanying Lecture 21
(Equilibrium figures), examples/examples10.tex.

Every quantitative claim in the solutions is checked here and reported as an
OK/FAIL line.  Only numpy and scipy are used.
"""

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

G = 6.674e-11                    # m^3 kg^-1 s^-2
b = 6.371e6                      # mean radius of the Earth, m
M_E = 5.972e24                   # mass of the Earth, kg
rho_bar = 5515.0                 # mean density, kg m^-3
Omega = 7.2921e-5                # angular velocity, s^-1
GM = G * M_E

fails = 0


def check(name, value, expected, rtol=1e-6, atol=0.0):
    global fails
    ok = abs(value - expected) <= atol + rtol * abs(expected)
    tag = "OK  " if ok else "FAIL"
    if not ok:
        fails += 1
    print(f"{tag} {name}: {value:.6g}  (expected {expected:.6g})")


def P2(x):
    return 0.5 * (3.0 * x**2 - 1.0)


print("=" * 72)
print("Example 1: homogeneous non-rotating planet")
print("=" * 72)


def p_hom(r, rho, b):
    return (2.0 / 3.0) * np.pi * G * rho**2 * (b**2 - r**2)


# check the closed form against direct quadrature of dp/dr = -rho g_0
for r in [0.0, 0.3 * b, 0.8 * b]:
    g0 = lambda s: (4.0 / 3.0) * np.pi * G * rho_bar * s
    pq, _ = quad(lambda s: rho_bar * g0(s), r, b)
    check(f"p(r={r/b:.1f}b) closed form vs quadrature", p_hom(r, rho_bar, b), pq)

pc_hom = p_hom(0.0, rho_bar, b)
check("central pressure, homogeneous Earth (Pa)", pc_hom, 1.7256e11, rtol=1e-3)
M_hom = (4.0 / 3.0) * np.pi * rho_bar * b**3
check("mass of homogeneous model (kg)", M_hom, M_E, rtol=3e-3)
check("p_c = 3GM^2/(8 pi b^4)", 3 * G * M_hom**2 / (8 * np.pi * b**4), pc_hom)
g_surf = G * M_hom / b**2
check("p_c = 3 g^2/(8 pi G)", 3 * g_surf**2 / (8 * np.pi * G), pc_hom)
print(f"     ratio PREM/homogeneous = {3.64e11/pc_hom:.3f}")
check("PREM central pressure / homogeneous value", 3.64e11 / pc_hom, 2.11, rtol=5e-3)
# mean pressure of the homogeneous model, and the interpretation p_c = (3/2) mean
p_mean = quad(lambda s: p_hom(s, rho_bar, b) * 4 * np.pi * s**2, 0, b)[0] / ((4 / 3) * np.pi * b**3)
check("p_c / <p> for homogeneous sphere", pc_hom / p_mean, 2.5)

print()
print("=" * 72)
print("Example 2: two-layer planet")
print("=" * 72)

rho_c, rho_m, c = 11.0e3, 4.5e3, 3.48e6


def mass_2(r):
    if r <= c:
        return (4.0 / 3.0) * np.pi * rho_c * r**3
    return (4.0 / 3.0) * np.pi * (rho_c * c**3 + rho_m * (r**3 - c**3))


def g_2(r):
    return G * mass_2(r) / r**2 if r > 0 else 0.0


def rho_2(r):
    return rho_c if r <= c else rho_m


def p_2(r):
    """closed forms quoted in the solution"""
    A = (4.0 / 3.0) * np.pi * G
    if r >= c:
        return A * rho_m * (0.5 * rho_m * (b**2 - r**2)
                            + (rho_c - rho_m) * c**3 * (1.0 / r - 1.0 / b))
    p_cmb = p_2(c)
    return p_cmb + (2.0 / 3.0) * np.pi * G * rho_c**2 * (c**2 - r**2)


M2 = mass_2(b)
I2 = (8.0 / 15.0) * np.pi * (rho_c * c**5 + rho_m * (b**5 - c**5))
I2_quad = quad(lambda s: (8.0 / 3.0) * np.pi * rho_2(s) * s**4, 0, b, points=[c])[0]
check("I closed form vs quadrature", I2, I2_quad)
print(f"     M = {M2:.4e} kg   (Earth 5.972e24, ratio {M2/M_E:.4f})")
check("mass of two-layer model (kg)", M2, 6.022e24, rtol=1e-3)
print(f"     I/(M b^2) = {I2/(M2*b**2):.4f}   (Earth 0.3307)")
check("I/(M b^2)", I2 / (M2 * b**2), 0.3465, rtol=1e-3)
g_cmb_minus = G * (4 / 3) * np.pi * rho_c * c
check("g continuous at CMB", g_2(c * (1 + 1e-12)), g_cmb_minus, rtol=1e-8)
print(f"     g(CMB) = {g_2(c):.3f} m s^-2,  g(b) = {g_2(b):.3f} m s^-2")
check("g(CMB) (m s^-2)", g_2(c), 10.70, rtol=1e-3)
check("g(b)  (m s^-2)", g_2(b), 9.90, rtol=1e-3)
# closed forms against quadrature
for r in [0.0, 0.5 * c, c, 0.5 * (b + c), 0.95 * b]:
    pq = quad(lambda s: rho_2(s) * g_2(s), r, b, points=[c] if r < c else None, limit=200)[0]
    check(f"p(r={r/1e3:7.0f} km) closed form vs quadrature", p_2(r), pq, rtol=1e-7)
print(f"     p(CMB) = {p_2(c):.4e} Pa  (PREM 1.36e11)")
print(f"     p(0)   = {p_2(0):.4e} Pa  (PREM 3.64e11)")
check("p(CMB) (Pa)", p_2(c), 1.2555e11, rtol=1e-3)
check("p(0)   (Pa)", p_2(0), 3.304e11, rtol=1e-3)
# the two contributions to p(CMB) quoted in the text
term1 = (2 / 3) * np.pi * G * rho_m**2 * (b**2 - c**2)
term2 = (4 / 3) * np.pi * G * rho_m * (rho_c - rho_m) * c**2 * (1 - c / b)
check("p(CMB) split: self-weight term (Pa)", term1, 8.06e10, rtol=2e-3)
check("p(CMB) split: excess-core term (Pa)", term2, 4.49e10, rtol=2e-3)
check("core contribution p(0)-p(CMB) (Pa)", p_2(0) - p_2(c), 2.048e11, rtol=1e-3)

print()
print("=" * 72)
print("Example 3: centrifugal potential and Y_20")
print("=" * 72)

rng = np.random.default_rng(1)
Om = np.array([0.0, 0.0, Omega])
for _ in range(5):
    x = rng.normal(size=3) * b
    r = np.linalg.norm(x)
    cth = x[2] / r
    psi_def = 0.5 * (np.outer(Om, Om) - Omega**2 * np.eye(3)) @ x @ x
    psi_fmt = -(1 / 3) * Omega**2 * r**2 + (1 / 3) * Omega**2 * r**2 * P2(cth)
    check("psi tensor definition vs P_2 form", psi_def, psi_fmt, rtol=1e-12)
    check("psi = -(1/2) Omega^2 r^2 sin^2 theta", psi_def, -0.5 * Omega**2 * (x[0]**2 + x[1]**2), rtol=1e-12)
    # gradient of psi reproduces the centrifugal acceleration eps_ijk eps_klm Om_j Om_l x_m
    h = 1.0
    grad = np.array([(0.5 * (np.outer(Om, Om) - Omega**2 * np.eye(3)) @ (x + h * e) @ (x + h * e)
                      - 0.5 * (np.outer(Om, Om) - Omega**2 * np.eye(3)) @ (x - h * e) @ (x - h * e)) / (2 * h)
                     for e in np.eye(3)])
    centrif = np.cross(Om, np.cross(Om, x))
    check("grad psi = Omega x (Omega x x)", np.linalg.norm(grad - centrif), 0.0, atol=1e-9 * np.linalg.norm(centrif))
# Laplacian of psi
check("Laplacian of psi = -2 Omega^2 (from -(1/2)Omega^2(x^2+y^2))", -0.5 * Omega**2 * 4, -2 * Omega**2, rtol=1e-14)
# Y_20 normalisation with Y_20 = sqrt(5/4pi) P_2(cos theta)
norm = quad(lambda t: (5 / (4 * np.pi)) * P2(np.cos(t))**2 * 2 * np.pi * np.sin(t), 0, np.pi)[0]
check("int |Y_20|^2 dOmega", norm, 1.0)
mean_P2 = quad(lambda t: P2(np.cos(t)) * np.sin(t), 0, np.pi)[0] / 2
check("spherical average of P_2", mean_P2, 0.0, atol=1e-14)
# ratio of aspherical to spherical part at a fixed r:  -P_2 / 1
# size of the terms at the Earth's surface
print(f"     (1/3) Omega^2 b^2 = {Omega**2*b**2/3:.4e} m^2 s^-2,  GM/b = {GM/b:.4e} m^2 s^-2,"
      f"  ratio = {Omega**2*b**3/(3*GM):.3e}")
check("ratio (1/3)Omega^2 b^3/(GM) = m/3", Omega**2 * b**3 / (3 * GM), 1.15e-3, rtol=2e-3)
# surface shape r = b[1 - (2/3) eps P_2]: equatorial minus polar radius = eps b, mean radius = b
eps = 1 / 298.257
r_eq = b * (1 - (2 / 3) * eps * P2(0.0))
r_po = b * (1 - (2 / 3) * eps * P2(1.0))
check("(r_eq - r_pol)/b = eps", (r_eq - r_po) / b, eps, rtol=1e-12)
check("r_eq = b(1+eps/3)", r_eq, b * (1 + eps / 3), rtol=1e-12)
check("r_pol = b(1-2eps/3)", r_po, b * (1 - 2 * eps / 3), rtol=1e-12)
r_mean = quad(lambda t: b * (1 - (2 / 3) * eps * P2(np.cos(t))) * np.sin(t), 0, np.pi)[0] / 2
check("mean radius of the spheroid = b", r_mean, b, rtol=1e-12)
vol = quad(lambda t: (1 / 3) * (b * (1 - (2 / 3) * eps * P2(np.cos(t))))**3 * 2 * np.pi * np.sin(t), 0, np.pi)[0]
check("volume unchanged to first order (relative change O(eps^2))",
      abs(vol / ((4 / 3) * np.pi * b**3) - 1) / eps**2, 0.0, atol=1.0)

print()
print("=" * 72)
print("Example 4: Clairaut's equation for a homogeneous planet")
print("=" * 72)

m = Omega**2 * b**3 / GM
print(f"     m = Omega^2 b^3/(GM) = {m:.5e}")
check("m for the Earth", m, 3.45e-3, rtol=2e-3)
check("m = 3 Omega^2/(4 pi G rho_bar) (homogeneous model)", 3 * Omega**2 / (4 * np.pi * G * M_E / ((4 / 3) * np.pi * b**3)), m, rtol=1e-12)
eps_N = 1.25 * m
print(f"     Newton: eps = 5m/4 = {eps_N:.5e},  1/eps = {1/eps_N:.1f}   (observed 298.26)")
check("1/eps for homogeneous Earth", 1 / eps_N, 231.9, rtol=1e-3)
check("1/eps Roche limit (m/2)", 2 / m, 579.7, rtol=1e-3)


def clairaut_rhs(r, y, rho0, g0):
    """y = (eps, eps').  Clairaut: eps'' + 8 pi G rho0/g0 (eps' + eps/r) - 6 eps/r^2 = 0."""
    e, de = y
    k = 8 * np.pi * G * rho0(r) / g0(r)
    return [de, -k * (de + e / r) + 6 * e / r**2]


def shoot(rho0, g0, Mtot, breaks=()):
    """Regular solution normalised to eps(r0)=1, then scaled to meet the surface condition
    eps'(b) = (1/b)[5 Omega^2 b^3/(2 G M) - 2 eps(b)].  Returns eps(b), eps'(b) and a callable eps(r)."""
    r0 = 1.0e-3
    y = np.array([1.0, 0.0])
    edges = [r0] + sorted(breaks) + [b]
    pieces = []
    for ra, rb in zip(edges[:-1], edges[1:]):
        sol = solve_ivp(clairaut_rhs, (ra, rb), y, args=(rho0, g0), method="DOP853",
                        rtol=1e-12, atol=1e-14, dense_output=True)
        y = sol.y[:, -1]
        pieces.append((ra, rb, sol.sol))
    et_b, det_b = y
    mm = Omega**2 * b**3 / (G * Mtot)
    lam = 2.5 * mm / (b * det_b + 2 * et_b)

    def eps_of_r(r):
        for ra, rb, f in pieces:
            if ra <= r <= rb:
                return lam * f(r)[0]
        return lam * pieces[0][2](r0)[0]

    return lam * et_b, lam * det_b, eps_of_r


# (a) constant density: coefficient 8 pi G rho0/g0 = 6/r
rho0_h = lambda r: rho_bar
g0_h = lambda r: (4 / 3) * np.pi * G * rho_bar * r
for r in [0.1 * b, 0.5 * b, b]:
    check(f"8 pi G rho0/g0 = 6/r at r={r/b:.1f}b", 8 * np.pi * G * rho0_h(r) / g0_h(r), 6 / r, rtol=1e-12)
# (b) direct substitution: eps = const and eps = r^-5 both solve the reduced ODE
for r in [0.2 * b, 0.7 * b]:
    res_const = clairaut_rhs(r, [1.0, 0.0], rho0_h, g0_h)[1]
    check(f"residual of eps=const at r={r/b:.1f}b (scaled by 6/r^2)", res_const / (6 / r**2), 0.0, atol=1e-12)
    e, de, dde = r**-5, -5 * r**-6, 30 * r**-7
    res_pow = dde + (6 / r) * de
    check(f"eps = r^-5 solves eps''+6eps'/r=0 (r={r/b:.1f}b)", res_pow / dde, 0.0, atol=1e-12)
# (c) shooting solution of the full ODE with the lecture's surface condition
eps_b_h, deps_b, eps_fun_h = shoot(rho0_h, g0_h, M_hom)
check("shooting: eps(b) for constant rho0 equals 5m/4", eps_b_h, 1.25 * Omega**2 * b**3 / (G * M_hom), rtol=1e-9)
r_grid = np.linspace(1e-3, b, 400)
check("shooting: eps constant through the interior (max |eps(r)/eps(b) - 1|)",
      np.max(np.abs(np.array([eps_fun_h(r) for r in r_grid]) / eps_b_h - 1)), 0.0, atol=1e-9)
check("shooting: eps'(b) = 0", deps_b * b / eps_b_h, 0.0, atol=1e-9)
# (d) independent check via the potential of the equatorial bulge (integral form):
#     (8/9) pi G rho0 eps = (8/15) pi G rho0 eps + Omega^2/3  ->  eps = 15 Omega^2/(16 pi G rho0) = 5m/4
eps_int = 15 * Omega**2 / (16 * np.pi * G * rho_bar)
check("integral (bulge-potential) form gives 5m/4", eps_int, 1.25 * Omega**2 * b**3 / (G * M_hom), rtol=1e-12)
# The bulge's self-gravity term: without it (Roche limit) eps = m/2
check("without bulge self-gravity: eps = 3 Omega^2/(8 pi G rho0) = m/2", 3 * Omega**2 / (8 * np.pi * G * rho_bar), 0.5 * Omega**2 * b**3 / (G * M_hom), rtol=1e-12)

# (e) Radau-Darwin approximation as an illustration of the trend with central condensation
def radau_darwin(mm, Ifac):
    return 2.5 * mm / (1 + 6.25 * (1 - 1.5 * Ifac)**2)


check("Radau-Darwin reproduces 5m/4 for I/(Mb^2)=0.4", radau_darwin(m, 0.4), 1.25 * m, rtol=1e-12)
print(f"     Radau-Darwin with I/(Mb^2)=0.3307: 1/eps = {1/radau_darwin(m, 0.3307):.1f}")
check("Radau-Darwin 1/eps for the Earth (~300)", 1 / radau_darwin(m, 0.3307), 300.0, rtol=2e-3)

# (f) two-layer model of Example 2: shooting vs exact algebraic solution
rho0_2 = rho_2
eps_b_2, deps_b_2, eps_fun_2 = shoot(rho0_2, g_2, M2, breaks=(c,))
m2 = Omega**2 * b**3 / (G * M2)
# algebraic solution: eps constant (= e_c) in the core; in the mantle
# (8/9) D(r) r^2 eps(r) = -(8/15)[c^5 e_c (rho_m-rho_c) r^-3 - r^2 e_b rho_m] + Omega^2 r^2/(3 pi G)
D_b = 3 * M2 / (4 * np.pi * b**3)
A = np.array([[(8 / 9) * rho_c - (8 / 15) * (rho_c - rho_m), -(8 / 15) * rho_m],
              [(8 / 15) * (rho_m - rho_c) * (c / b)**5, (8 / 9) * D_b - (8 / 15) * rho_m]])
rhs = np.array([Omega**2 / (3 * np.pi * G), Omega**2 / (3 * np.pi * G)])
e_c, e_b = np.linalg.solve(A, rhs)
check("two-layer: shooting eps(b) vs exact", eps_b_2, e_b, rtol=1e-8)
check("two-layer: shooting eps(c) vs exact", eps_fun_2(c), e_c, rtol=1e-8)
check("two-layer: shooting eps(0) = eps(c) (constant in core)", eps_fun_2(0.3 * c), e_c, rtol=1e-8)
check("two-layer: surface condition satisfied", deps_b_2, (2.5 * m2 - 2 * eps_b_2) / b, rtol=1e-8)


def eps_mantle_exact(r):
    D = 3 * mass_2(r) / (4 * np.pi * r**3)
    return (-(8 / 15) * (c**5 * e_c * (rho_m - rho_c) / r**3 - r**2 * e_b * rho_m) + Omega**2 * r**2 / (3 * np.pi * G)) / ((8 / 9) * D * r**2)


for r in [0.6 * b, 0.8 * b]:
    check(f"two-layer: shooting eps(r={r/b:.1f}b) vs exact", eps_fun_2(r), eps_mantle_exact(r), rtol=1e-8)
print(f"     two-layer model: 1/eps(b) = {1/e_b:.1f},  1/eps(CMB) = {1/e_c:.1f},  eps(CMB)/eps(b) = {e_c/e_b:.4f}")
check("two-layer 1/eps(b)", 1 / e_b, 285.4, rtol=2e-3)
check("two-layer 1/eps(CMB)", 1 / e_c, 374.3, rtol=2e-3)
eps_prof_2 = np.array([eps_fun_2(r) for r in np.linspace(1e-3, b, 500)])
check("two-layer: eps increases outward", float(np.all(np.diff(eps_prof_2) >= -1e-12)), 1.0)
check("two-layer: eps(b) between Roche and Newton limits", float(0.5 * m2 < e_b < 1.25 * m2), 1.0)
print(f"     Newton 1/(5m2/4) = {1/(1.25*m2):.1f}, Roche 1/(m2/2) = {2/m2:.1f}")

print()
print("=" * 72)
print("Example 5: coincidence of level surfaces")
print("=" * 72)

# smooth test model rho0 = rho_c0 (1 - alpha r^2/b^2)
rho_c0, alpha = 13.0e3, 10.0 / 13.0
rho0_s = lambda r: rho_c0 * (1 - alpha * r**2 / b**2)
drho0_s = lambda r: -2 * rho_c0 * alpha * r / b**2
mass_s = lambda r: 4 * np.pi * rho_c0 * (r**3 / 3 - alpha * r**5 / (5 * b**2))
g0_s = lambda r: G * mass_s(r) / r**2
phi0_s = lambda r: quad(g0_s, 0.5 * b, r)[0]          # potential up to an irrelevant constant
p0_s = lambda r: quad(lambda s: rho0_s(s) * g0_s(s), r, b)[0]
dp0_s = lambda r: -rho0_s(r) * g0_s(r)

# first-order gravity potential perturbation of Y_20 type (amplitude chosen so that s*gamma1/(r g0) ~ eps)
gamma1 = lambda r, cth: (2 / 3) * r * g0_s(r) * P2(cth) * 1.0e-3
p1 = lambda r, cth: -rho0_s(r) * gamma1(r, cth)
rho1 = lambda r, cth: drho0_s(r) / g0_s(r) * gamma1(r, cth)

# radial component of the first-order equation: d_r p1 + rho0 d_r gamma1 + rho1 g0 = 0
for (r, cth) in [(0.4 * b, 0.3), (0.9 * b, -0.8)]:
    h = 1.0
    dp1 = (p1(r + h, cth) - p1(r - h, cth)) / (2 * h)
    dg1 = (gamma1(r + h, cth) - gamma1(r - h, cth)) / (2 * h)
    res = dp1 + rho0_s(r) * dg1 + rho1(r, cth) * g0_s(r)
    check(f"radial first-order balance at r={r/b:.1f}b (relative)", res / (rho0_s(r) * g0_s(r) * 1e-3), 0.0, atol=1e-6)
    # tangential component: d_theta (p1 + rho0 gamma1) = 0 trivially since p1 = -rho0 gamma1
    check("p1 + rho0 gamma1 = 0", p1(r, cth) + rho0_s(r) * gamma1(r, cth), 0.0, atol=1e-9)

# displaced level surfaces, found by root-finding, versus -gamma1/g0
for (rc, cth) in [(0.4 * b, 0.3), (0.9 * b, -0.8)]:
    predicted = -gamma1(rc, cth) / g0_s(rc)
    for s in [1.0, 0.5]:
        p_tot = lambda r: p0_s(r) + s * p1(r, cth)
        rho_tot = lambda r: rho0_s(r) + s * rho1(r, cth)
        gam_tot = lambda r: phi0_s(r) + s * gamma1(r, cth)
        dr_p = brentq(lambda r: p_tot(r) - p0_s(rc), 0.9 * rc, min(1.1 * rc, b), xtol=1e-9) - rc
        dr_rho = brentq(lambda r: rho_tot(r) - rho0_s(rc), 0.9 * rc, min(1.1 * rc, b), xtol=1e-9) - rc
        dr_gam = brentq(lambda r: gam_tot(r) - phi0_s(rc), 0.9 * rc, min(1.1 * rc, b), xtol=1e-9) - rc
        for name, dr in [("p", dr_p), ("rho", dr_rho), ("gamma", dr_gam)]:
            # first-order agreement: relative error should be O(s * eps) ~ 1e-3
            check(f"level surface of {name} at r={rc/b:.1f}b, s={s}: dr/(s*predicted)", dr / (s * predicted), 1.0, rtol=3e-3)
    # the three displacements agree with each other to higher accuracy than with the prediction
    check(f"dr_p = dr_gamma at r={rc/b:.1f}b (rel.)", dr_p / dr_gam, 1.0, rtol=1e-3)
    check(f"dr_rho = dr_gamma at r={rc/b:.1f}b (rel.)", dr_rho / dr_gam, 1.0, rtol=1e-3)
# geometric form of the displacement h = -(2/3) r eps P_2 with eps = gamma1/((2/3) r g0 P_2)
eps_loc = 1.0e-3
check("h = -gamma1/g0 = -(2/3) r eps P_2", -gamma1(0.5 * b, 0.3) / g0_s(0.5 * b), -(2 / 3) * 0.5 * b * eps_loc * P2(0.3), rtol=1e-12)

print()
print("=" * 72)
if fails == 0:
    print("ALL CHECKS PASSED")
else:
    print(f"{fails} CHECK(S) FAILED")

print()
print("=" * 72)
print("Additional numbers quoted in the text")
print("=" * 72)
check("r_eq - r_pol for eps = 1/298 (km)", b / 298 / 1e3, 21.4, rtol=2e-3)
check("homogeneous Earth over-flattened: 232 vs 298 -> eps ratio", 298 / 231.9, 1.285, rtol=2e-3)
check("two-layer p(CMB) short of PREM by ~8 per cent", 100 * (1 - p_2(c) / 1.36e11), 7.7, rtol=2e-2)
check("two-layer p(0)   short of PREM by ~9 per cent", 100 * (1 - p_2(0) / 3.64e11), 9.2, rtol=2e-2)
check("CMB equatorial minus polar radius in two-layer model (km)", c * e_c / 1e3, 9.3, rtol=2e-2)
check("(2/3) Omega^2 b^3/(GM) ~ 2e-3 (effective-gravity correction)", (2 / 3) * m, 2.3e-3, rtol=2e-2)
print()
print("ALL CHECKS PASSED" if fails == 0 else f"{fails} CHECK(S) FAILED")
