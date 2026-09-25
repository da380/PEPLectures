"""Numerical checks for the worked examples accompanying Lecture 16 (ray theory).

Every numerical or closed-form claim made in examples/examples5.tex is checked
here, printing OK/FAIL lines. Requires numpy and scipy only.
"""
import numpy as np
from scipy.integrate import solve_ivp

status = {"fail": 0}


def check(name, ok, detail=""):
    print(("OK   " if ok else "FAIL ") + name + ("  " + detail if detail else ""))
    if not ok:
        status["fail"] += 1


# ---------------------------------------------------------------------------
# Example 1: circular rays in alpha(z) = alpha0 + g z, z downwards
# ---------------------------------------------------------------------------
print("\nExample 1: linear velocity gradient")
alpha0, g = 5.0, 0.05            # km/s, 1/s
th0 = np.deg2rad(30.0)
px = np.sin(th0) / alpha0
R = 1.0 / (g * px)
zc = -alpha0 / g
xc = R * np.cos(th0)


def alpha(z):
    return alpha0 + g * z


def hamilton(sigma, y):
    x, z, p_x, p_z = y
    a = alpha(z)
    return [a * a * p_x, a * a * p_z, 0.0, -(p_x * p_x + p_z * p_z) * a * g]


def T_closed(theta):
    return np.log(np.tan(theta / 2) / np.tan(th0 / 2)) / g


# integrate from the origin until the ray returns to the surface (z = 0)
def hit_surface(sigma, y):
    return y[1]
hit_surface.terminal = True
hit_surface.direction = -1

T_surf = 2.0 * np.log(1.0 / np.tan(th0 / 2)) / g
sol = solve_ivp(hamilton, [0.0, 2.0 * T_surf], [0.0, 0.0, px, np.cos(th0) / alpha0],
                events=hit_surface, rtol=1e-11, atol=1e-13, dense_output=True)
sig = np.linspace(0.0, sol.t_events[0][0], 400)
x, z, p_x, p_z = sol.sol(sig)

check("p_x conserved along the ray", np.max(np.abs(p_x - px)) < 1e-12)
check("eikonal H = 1/2 conserved", np.max(np.abs(0.5 * alpha(z) ** 2 * (p_x ** 2 + p_z ** 2) - 0.5)) < 1e-9)
rad = np.sqrt((x - xc) ** 2 + (z - zc) ** 2)
check("ray is a circle of radius R = 1/(g p_x)", np.max(np.abs(rad - R)) < 1e-7,
      "R = %.3f km, centre (%.3f, %.3f) km" % (R, xc, zc))
theta = np.arctan2(p_x, p_z)     # angle of the slowness vector from the downward vertical
check("closed-form T(theta) equals generating parameter sigma",
      np.max(np.abs(T_closed(theta) - sig)) < 1e-8)
check("angle theta advances as d theta/d sigma = g sin theta",
      np.max(np.abs(np.gradient(theta, sig) - g * np.sin(theta))) < 1e-4)
X_num, T_num = sol.y_events[0][0][0], sol.t_events[0][0]
X_cf = 2.0 * alpha0 / (g * np.tan(th0))
z_t = alpha0 * (1.0 / np.sin(th0) - 1.0) / g
check("surface distance X = 2 alpha0 cot(theta0)/g", abs(X_num - X_cf) < 1e-6,
      "X = %.3f km" % X_cf)
check("surface travel time T = (2/g) ln cot(theta0/2)", abs(T_num - T_surf) < 1e-8,
      "T = %.4f s" % T_surf)
check("T = (2/g) arcsinh(g X / 2 alpha0)",
      abs(T_surf - 2.0 / g * np.arcsinh(g * X_cf / (2 * alpha0))) < 1e-12)
z_turn = sol.sol(T_closed(np.pi / 2))[1]
check("turning depth z_t = alpha0 (1/sin theta0 - 1)/g at theta = pi/2", abs(z_turn - z_t) < 1e-6,
      "z_t = %.3f km" % z_t)
check("small-offset limit T -> X/alpha0", abs(2.0 / g * np.arcsinh(g * 1e-3 / (2 * alpha0)) - 1e-3 / alpha0) < 1e-12)
# numbers quoted in the text
check("quoted numbers: R=200 km, z_t=100 km, X=346.4 km, T=52.68 s",
      abs(R - 200) < 1e-9 and abs(z_t - 100) < 1e-9 and abs(X_cf - 346.41) < 0.01 and abs(T_surf - 52.68) < 0.005)
check("straight-line time X/alpha0 = 69.28 s exceeds T", abs(X_cf / alpha0 - 69.28) < 0.005 and X_cf / alpha0 > T_surf)

# ---------------------------------------------------------------------------
# Example 2: ray equation and curvature in a 2D heterogeneous medium
# ---------------------------------------------------------------------------
print("\nExample 2: ray curvature formula")
def alpha2(x, z):
    return 6.0 + 0.03 * z - 1.5 * np.exp(-((x - 150.0) ** 2 + (z - 80.0) ** 2) / (2 * 40.0 ** 2))

def grad_alpha2(x, z, h=1e-4):
    return np.array([(alpha2(x + h, z) - alpha2(x - h, z)) / (2 * h),
                     (alpha2(x, z + h) - alpha2(x, z - h)) / (2 * h)])

def ham2(s, y):
    x, z, p_x, p_z = y
    a = alpha2(x, z)
    ga = grad_alpha2(x, z)
    n2 = p_x ** 2 + p_z ** 2
    return [a * a * p_x, a * a * p_z, -n2 * a * ga[0], -n2 * a * ga[1]]

th = np.deg2rad(55.0)
a00 = alpha2(0.0, 0.0)
sol2 = solve_ivp(ham2, [0.0, 60.0], [0.0, 0.0, np.sin(th) / a00, np.cos(th) / a00],
                 rtol=1e-11, atol=1e-13, dense_output=True)
sg = np.linspace(0.0, 60.0, 3000)
Y = sol2.sol(sg)
xx, zz = Y[0], Y[1]
a_on = alpha2(xx, zz)
dx, dz = a_on ** 2 * Y[2], a_on ** 2 * Y[3]        # dx/d sigma = alpha^2 p
speed = np.sqrt(dx ** 2 + dz ** 2)
check("ds/d sigma = alpha along the ray", np.max(np.abs(speed - a_on) / a_on) < 1e-8)
s = np.concatenate([[0.0], np.cumsum(0.5 * (speed[1:] + speed[:-1]) * np.diff(sg))])
tx, tz = dx / speed, dz / speed
kappa_num = np.sqrt(np.gradient(tx, s) ** 2 + np.gradient(tz, s) ** 2)
G = np.array([grad_alpha2(xi, zi) for xi, zi in zip(xx, zz)]) / a_on[:, None]
kappa_cf = np.abs(tx * G[:, 1] - tz * G[:, 0])       # |t x grad ln alpha| in 2D
inner = 200
check("curvature = |t x grad ln alpha|",
      np.max(np.abs(kappa_num - kappa_cf)[inner:-inner]) < 2e-4 * np.max(kappa_cf),
      "max kappa = %.4f /km" % np.max(kappa_cf))
# normal points towards decreasing alpha: n . grad alpha <= 0
nx, nz = np.gradient(tx, s) / kappa_num, np.gradient(tz, s) / kappa_num
proj = (nx * G[:, 0] + nz * G[:, 1])[inner:-inner]
check("principal normal points towards decreasing alpha", np.all(proj < 1e-6))
# consistency with Example 1: kappa = g p_x = 1/R
check("Example 1: |t x grad ln alpha| = g sin(theta)/alpha = 1/R",
      abs(g * np.sin(th0) / alpha0 - 1.0 / R) < 1e-15)

# ---------------------------------------------------------------------------
# Example 3: 1D transport equation, a0 ~ (rho c)^(-1/2)
# ---------------------------------------------------------------------------
print("\nExample 3: one-dimensional rod")
def rho3(x):
    return 2.0 + 0.8 * np.tanh((x - 5.0) / 1.5)

def E3(x):
    return 20.0 + 12.0 * np.tanh((x - 4.0) / 2.0) + 3.0 * np.sin(x / 2.0)

def c3(x):
    return np.sqrt(E3(x) / rho3(x))

xs = np.linspace(0.0, 10.0, 5)
check("E/c = rho c = sqrt(rho E)", np.max(np.abs(E3(xs) / c3(xs) - rho3(xs) * c3(xs))) < 1e-12)

# transport equation 2 E T' a' + (E T')' a = 0 satisfied by a = (rho c)^(-1/2), T' = 1/c
h = 1e-4
def a3(x):
    return (rho3(x) * c3(x)) ** -0.5

def ET(x):
    return E3(x) / c3(x)

res = 2 * ET(xs) * (a3(xs + h) - a3(xs - h)) / (2 * h) + (ET(xs + h) - ET(xs - h)) / (2 * h) * a3(xs)
scale = np.abs(2 * ET(xs) * (a3(xs + h) - a3(xs - h)) / (2 * h)).max()
check("a0 = (rho c)^(-1/2) satisfies 2 E T' a0' + (E T')' a0 = 0", np.max(np.abs(res)) < 1e-7 * scale)
# energy-flux form: E T' a0^2 = const
check("E T' (a0)^2 constant", np.ptp(ET(xs) * a3(xs) ** 2) < 1e-12)

# compare with numerical solution of -w^2 rho u - (E u')' = 0 at high frequency
from scipy.integrate import quad
def wkb(x, w):
    T = quad(lambda t: 1.0 / c3(t), 0.0, x, limit=200)[0]
    return a3(x) * np.exp(-1j * w * T)

def rod(x, y, w):
    u, q = y                       # q = E u'
    return [q / E3(x), -w * w * rho3(x) * u]

errs = []
for w in (20.0, 40.0, 80.0):
    # launch a rightward WKB wave at x = 0 (including the O(1/w) amplitude derivative)
    u0 = a3(0.0)
    du0 = (a3(h) - a3(-h)) / (2 * h) - 1j * w / c3(0.0) * a3(0.0)
    y0 = [complex(u0), complex(E3(0.0) * du0)]
    s3 = solve_ivp(rod, [0.0, 10.0], y0, args=(w,), rtol=1e-10, atol=1e-12, dense_output=True)
    xe = np.linspace(2.0, 10.0, 200)
    u_num = s3.sol(xe)[0]
    u_wkb = np.array([wkb(xi, w) for xi in xe])
    errs.append(np.max(np.abs(u_num - u_wkb)) / np.max(np.abs(u_wkb)))
check("zeroth-order ray solution agrees with exact rod solution at high frequency",
      errs[-1] < 0.02, "relative errors %s at omega = 20, 40, 80" % np.array2string(np.array(errs), precision=4))
check("error decreases roughly as 1/omega", errs[0] / errs[1] > 1.6 and errs[1] / errs[2] > 1.6)
# ray-theoretic vs exact transmitted amplitude across an impedance jump Z2 = Z1 (1 + delta)
d = 1e-3
ray_amp, exact_amp = (1.0 + d) ** -0.5, 2.0 / (2.0 + d)
check("sqrt(Z1/Z2) and 2Z1/(Z1+Z2) both equal 1 - delta/2 + O(delta^2)",
      abs(ray_amp - (1 - d / 2)) < d ** 2 and abs(exact_amp - (1 - d / 2)) < d ** 2 and abs(ray_amp - exact_amp) < d ** 2)

# ---------------------------------------------------------------------------
# Example 4: point source in a homogeneous 2D medium, transport equation
# ---------------------------------------------------------------------------
print("\nExample 4: point source and 2D amplitude")
al = 3.0
def Tp(x, y):
    return np.hypot(x, y) / al

pts = np.array([[1.0, 0.3], [-0.4, 2.0], [0.7, -0.9]])
for (x, y) in pts:
    gx = (Tp(x + h, y) - Tp(x - h, y)) / (2 * h)
    gy = (Tp(x, y + h) - Tp(x, y - h)) / (2 * h)
    check("eikonal |grad T| = 1/alpha at (%.1f, %.1f)" % (x, y), abs(np.hypot(gx, gy) - 1.0 / al) < 1e-8)

# transport identity div(A0^2 grad T) = 0 with A0 = r^(-1/2)
def flux(x, y, m):
    r = np.hypot(x, y)
    return r ** (-2 * m) * np.array([x, y]) / (r * al)

def div_flux(x, y, m):
    return (flux(x + h, y, m)[0] - flux(x - h, y, m)[0]) / (2 * h) + (flux(x, y + h, m)[1] - flux(x, y - h, m)[1]) / (2 * h)

check("div(A0^2 grad T) = 0 for A0 = r^(-1/2)", max(abs(div_flux(x, y, 0.5)) for x, y in pts) < 1e-6)
check("div(A0^2 grad T) != 0 for A0 = r^(-1)", min(abs(div_flux(x, y, 1.0)) for x, y in pts) > 1e-2)

# full tensor check of the first transport equation, eq. (19) of the lecture,
# for a homogeneous isotropic medium: p_hat . RHS = 0 with a0 = r^(-1/2) p_hat
lam, mu, rho0 = 4.0, 2.5, 1.2
alph = np.sqrt((lam + 2 * mu) / rho0)
d = np.eye(3)
A = lam * np.einsum("ij,kl->ijkl", d, d) + mu * (np.einsum("ik,jl->ijkl", d, d) + np.einsum("il,jk->ijkl", d, d))

def p_vec(X):
    r = np.hypot(X[0], X[1])
    return np.array([X[0], X[1], 0.0]) / (r * alph)

def a0_vec(X, m):
    r = np.hypot(X[0], X[1])
    return r ** (-m) * alph * p_vec(X)

def rhs19(X, m):
    F = np.zeros(3)
    for j in range(3):
        e = np.zeros(3); e[j] = h
        Xp, Xm = X + e, X - e
        # d/dx_j ( A_ijkl p_l a_k )
        F += (np.einsum("ikl,l,k->i", A[:, j], p_vec(Xp), a0_vec(Xp, m))
              - np.einsum("ikl,l,k->i", A[:, j], p_vec(Xm), a0_vec(Xm, m))) / (2 * h)
    for l in range(3):
        e = np.zeros(3); e[l] = h
        da = (a0_vec(X + e, m) - a0_vec(X - e, m)) / (2 * h)
        F += np.einsum("ijk,j,k->i", A[:, :, :, l], p_vec(X), da)
    return F

X0 = np.array([0.8, 1.1, 0.0])
ph = alph * p_vec(X0)
proj_half = ph @ rhs19(X0, 0.5)
proj_one = ph @ rhs19(X0, 1.0)
check("p_hat . RHS of first transport equation vanishes for A0 = r^(-1/2)",
      abs(proj_half) < 1e-6 * np.linalg.norm(rhs19(X0, 1.0)), "value %.2e" % proj_half)
check("... but not for A0 = r^(-1)", abs(proj_one) > 1e-2, "value %.3f" % proj_one)
check("symmetry d p_i/d x_j = d p_j/d x_i",
      abs((p_vec(X0 + [h, 0, 0])[1] - p_vec(X0 - [h, 0, 0])[1]) / (2 * h)
          - (p_vec(X0 + [0, h, 0])[0] - p_vec(X0 - [0, h, 0])[0]) / (2 * h)) < 1e-8)

# ---------------------------------------------------------------------------
# Example 5: size of the first correction term
# ---------------------------------------------------------------------------
print("\nExample 5: size of the neglected term")
L, aP, period = 1000.0, 10.0, 20.0
lam_w = aP * period
omega = 2 * np.pi / period
tau = L / aP
ratio = 1.0 / (omega * tau)
check("lambda = 200 km, tau = 100 s, omega tau = 31.4", abs(lam_w - 200) < 1e-9 and abs(tau - 100) < 1e-9 and abs(omega * tau - 31.42) < 0.005)
check("1/(omega tau) = lambda/(2 pi L) = 0.032", abs(ratio - lam_w / (2 * np.pi * L)) < 1e-15 and abs(ratio - 0.0318) < 5e-4,
      "ratio = %.4f" % ratio)
check("L = 50 km gives ratio 0.64", abs(lam_w / (2 * np.pi * 50.0) - 0.637) < 5e-4)

print("\n%d check(s) failed" % status["fail"] if status["fail"] else "\nAll checks passed")
