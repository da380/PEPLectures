"""Numerical checks for examples/examples8.tex (Worked examples for Lecture 19:
Waveform tomography).

Every quantitative claim in the solutions is checked here and an OK/FAIL line is
printed.  Only numpy and scipy are used.
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

np.set_printoptions(precision=10)
FAILS = 0


def check(name, value, target, tol, rel=True):
    global FAILS
    err = abs(value - target)
    if rel and target != 0:
        err = err / abs(target)
    ok = err <= tol
    if not ok:
        FAILS += 1
    print(f"{'OK  ' if ok else 'FAIL'} {name}: value={value!r} target={target!r} "
          f"{'rel' if rel else 'abs'}_err={err:.3e} tol={tol:.1e}")
    return ok


print("=" * 72)
print("Example 1: adjoint method for a linear system")
print("=" * 72)

# --- (a) the 2x2 example worked in the text -------------------------------
def example1_closed_form(m, d):
    A = np.array([[2.0, -1.0], [-1.0, m]])
    f = np.array([1.0, 0.0])
    P = np.array([[0.0, 1.0]])
    u = np.linalg.solve(A, f)
    r = P @ u - d
    J = 0.5 * r @ r
    # adjoint solve
    up = np.linalg.solve(A.T, P.T @ r)
    dA = np.array([[0.0, 0.0], [0.0, 1.0]])
    grad_adj = -up @ dA @ u
    return u, up, J, grad_adj


m0, d0 = 1.0, 0.5
u, up, J, g_adj = example1_closed_form(m0, d0)
check("Ex1 u = (1,1)", np.linalg.norm(u - np.array([1.0, 1.0])), 0.0, 1e-13, rel=False)
check("Ex1 u' = (0.5,1)", np.linalg.norm(up - np.array([0.5, 1.0])), 0.0, 1e-13, rel=False)
check("Ex1 J(m=1) = 1/8", J, 0.125, 1e-13)
check("Ex1 adjoint dJ/dm = -1", g_adj, -1.0, 1e-13)
# direct closed form dJ/dm = -2 (u2 - d)/(2m-1)^2
u2 = 1.0 / (2 * m0 - 1)
check("Ex1 closed form dJ/dm", -2 * (u2 - d0) / (2 * m0 - 1) ** 2, g_adj, 1e-13)
h = 1e-6
Jp = example1_closed_form(m0 + h, d0)[2]
Jm = example1_closed_form(m0 - h, d0)[2]
check("Ex1 central FD dJ/dm", (Jp - Jm) / (2 * h), g_adj, 1e-8)

# --- (b) a random larger system: adjoint gradient vs finite differences ----
rng = np.random.default_rng(1)
n, nm, nobs = 7, 4, 3
A0 = rng.standard_normal((n, n)) + 4 * np.eye(n)
B = [rng.standard_normal((n, n)) for _ in range(nm)]
f = rng.standard_normal(n)
P = rng.standard_normal((nobs, n))
d = rng.standard_normal(nobs)


def Amat(mvec):
    return A0 + sum(mk * Bk for mk, Bk in zip(mvec, B))


def Jfun(mvec):
    uu = np.linalg.solve(Amat(mvec), f)
    r = P @ uu - d
    return 0.5 * r @ r


mvec = rng.standard_normal(nm) * 0.3
A = Amat(mvec)
uu = np.linalg.solve(A, f)
r = P @ uu - d
upp = np.linalg.solve(A.T, P.T @ r)
g_adj = np.array([-upp @ Bk @ uu for Bk in B])
# direct method: m solves of A du_k = -B_k u
g_dir = np.array([r @ P @ np.linalg.solve(A, -Bk @ uu) for Bk in B])
h = 1e-6
g_fd = np.array([(Jfun(mvec + h * e) - Jfun(mvec - h * e)) / (2 * h) for e in np.eye(nm)])
check("Ex1 random system: adjoint vs direct", np.max(abs(g_adj - g_dir)), 0.0, 1e-12, rel=False)
check("Ex1 random system: adjoint vs FD", np.max(abs(g_adj - g_fd) / abs(g_fd)), 0.0, 1e-7, rel=False)

print()
print("=" * 72)
print("Example 2: adjoint method for a damped oscillator")
print("=" * 72)

T = 10.0
omega0, gamma = 2.0, 0.1              # parameters at which the gradient is evaluated
omega0_true, gamma_true = 2.2, 0.15   # parameters used to generate x_obs
RTOL, ATOL = 1e-12, 1e-14


def force(t):
    return np.exp(-(t - 2.0) ** 2)


def forward(om2, gam, dense=True):
    """Solve x'' + 2 gam x' + om2 x = f, x(0)=x'(0)=0 on [0,T]."""
    def rhs(t, y):
        return [y[1], force(t) - 2 * gam * y[1] - om2 * y[0]]
    return solve_ivp(rhs, (0.0, T), [0.0, 0.0], method="DOP853",
                     rtol=RTOL, atol=ATOL, dense_output=dense)


sol_obs = forward(omega0_true ** 2, gamma_true)
x_obs = lambda t: sol_obs.sol(t)[0]


def misfit_and_adjoint(om2, gam):
    """Return J, dJ/d(omega0^2), dJ/dgamma, dJ/df-kernel norm via the adjoint method."""
    sol = forward(om2, gam)
    x = lambda t: sol.sol(t)[0]
    xdot = lambda t: sol.sol(t)[1]
    # misfit
    J = quad(lambda t: 0.5 * (x(t) - x_obs(t)) ** 2, 0.0, T, epsabs=1e-14, epsrel=1e-13, limit=500)[0]
    # adjoint: x'' - 2 gam x' + om2 x = x - x_obs, x(T)=x'(T)=0, integrated backwards;
    # carry the two kernel integrals I1 = int x x' dt, I2 = int xdot x' dt as extra states.
    def rhs_adj(t, y):
        xp, xpd, I1, I2 = y
        return [xpd, (x(t) - x_obs(t)) + 2 * gam * xpd - om2 * xp,
                x(t) * xp, xdot(t) * xp]
    sol_adj = solve_ivp(rhs_adj, (T, 0.0), [0.0, 0.0, 0.0, 0.0], method="DOP853",
                        rtol=RTOL, atol=ATOL, dense_output=True)
    # integrating from T down to 0: I(0) = -int_0^T (...) dt
    I1 = -sol_adj.y[2, -1]
    I2 = -sol_adj.y[3, -1]
    return J, -I1, -2 * I2, sol, sol_adj


def misfit_only(om2, gam):
    sol = forward(om2, gam)
    x = lambda t: sol.sol(t)[0]
    return quad(lambda t: 0.5 * (x(t) - x_obs(t)) ** 2, 0.0, T, epsabs=1e-14, epsrel=1e-13, limit=500)[0]


J0, dJ_dom2, dJ_dgam, sol_fwd, sol_adj = misfit_and_adjoint(omega0 ** 2, gamma)
print(f"J(omega0=2, gamma=0.1)        = {J0:.10e}")
print(f"adjoint dJ/d(omega0^2)        = {dJ_dom2:.10e}")
print(f"adjoint dJ/dgamma             = {dJ_dgam:.10e}")

# central finite differences with two step sizes
for hh in (1e-3, 1e-4):
    fd_om2 = (misfit_only(omega0 ** 2 + hh, gamma) - misfit_only(omega0 ** 2 - hh, gamma)) / (2 * hh)
    fd_gam = (misfit_only(omega0 ** 2, gamma + hh) - misfit_only(omega0 ** 2, gamma - hh)) / (2 * hh)
    print(f"  h={hh:.0e}: FD dJ/d(omega0^2) = {fd_om2:.10e}   FD dJ/dgamma = {fd_gam:.10e}")
    # central differences carry an O(h^2) truncation error, so the tolerance scales with h^2
    tol = 1e-4 if hh == 1e-3 else 2e-6
    check(f"Ex2 dJ/d(omega0^2) adjoint vs central FD (h={hh:.0e})", dJ_dom2, fd_om2, tol)
    check(f"Ex2 dJ/dgamma adjoint vs central FD (h={hh:.0e})", dJ_dgam, fd_gam, tol)

# also check dJ/domega0 = 2 omega0 dJ/d(omega0^2)
hh = 1e-4
fd_om = (misfit_only((omega0 + hh) ** 2, gamma) - misfit_only((omega0 - hh) ** 2, gamma)) / (2 * hh)
check("Ex2 dJ/domega0 = 2 omega0 dJ/d(omega0^2)", 2 * omega0 * dJ_dom2, fd_om, 2e-6)

# terminal conditions and the time-reversed form: y(s) = x'(T - s) satisfies the
# *forward* damped oscillator driven by the time-reversed residual
x_f = lambda t: sol_fwd.sol(t)[0]
def rhs_rev(s, y):
    t = T - s
    return [y[1], (x_f(t) - x_obs(t)) - 2 * gamma * y[1] - omega0 ** 2 * y[0]]
sol_rev = solve_ivp(rhs_rev, (0.0, T), [0.0, 0.0], method="DOP853", rtol=RTOL, atol=ATOL, dense_output=True)
smp = np.linspace(0, T, 11)
diff = max(abs(sol_rev.sol(s)[0] - sol_adj.sol(T - s)[0]) for s in smp)
check("Ex2 time-reversed adjoint = forward oscillator with reversed source", diff, 0.0, 1e-9, rel=False)
print(f"  x'(0) = {sol_adj.sol(0.0)[0]:.6e}, x'(T) = {sol_adj.sol(T)[0]:.1e}, "
      f"max|x| = {max(abs(x_f(t)) for t in np.linspace(0, T, 2001)):.5f}, "
      f"max|x'| = {max(abs(sol_adj.sol(t)[0]) for t in np.linspace(0, T, 2001)):.5e}")

# source kernel: dJ/d(epsilon) for f -> f + eps*phi equals int x' phi dt
phi = lambda t: np.sin(3 * t)
src_adj = quad(lambda t: sol_adj.sol(t)[0] * phi(t), 0.0, T, epsabs=1e-14, epsrel=1e-13, limit=500)[0]
def misfit_force(eps):
    def rhs(t, y):
        return [y[1], force(t) + eps * phi(t) - 2 * gamma * y[1] - omega0 ** 2 * y[0]]
    s = solve_ivp(rhs, (0.0, T), [0.0, 0.0], method="DOP853", rtol=RTOL, atol=ATOL, dense_output=True)
    return quad(lambda t: 0.5 * (s.sol(t)[0] - x_obs(t)) ** 2, 0.0, T, epsabs=1e-14, epsrel=1e-13, limit=500)[0]
hh = 1e-4
src_fd = (misfit_force(hh) - misfit_force(-hh)) / (2 * hh)
print(f"  source kernel: adjoint {src_adj:.10e}   FD {src_fd:.10e}")
check("Ex2 source sensitivity int x' phi dt vs FD", src_adj, src_fd, 2e-6)

print()
print("=" * 72)
print("Example 3: steepest descent on a quadratic")
print("=" * 72)

H = np.diag([1.0, 9.0])
b = np.zeros(2)
kappa = 9.0
rate = ((kappa - 1) / (kappa + 1)) ** 2
check("Ex3 ((kappa-1)/(kappa+1))^2 for kappa=9", rate, 0.64, 1e-15)

m = np.array([9.0, 1.0])
Jq = lambda mm: 0.5 * mm @ H @ mm - b @ mm
J_hist = [Jq(m)]
m_hist = [m.copy()]
for k in range(40):
    g = H @ m - b
    lam = (g @ g) / (g @ H @ g)
    m = m - lam * g
    m_hist.append(m.copy())
    J_hist.append(Jq(m))
m_hist = np.array(m_hist)
J_hist = np.array(J_hist)
# closed form: m_k = 0.8^k (9, (-1)^k), J_k = 45 * 0.64^k, lambda = 0.2 every step
ks = np.arange(41)
m_exact = np.stack([0.8 ** ks * 9.0, 0.8 ** ks * (-1.0) ** ks], axis=1)
check("Ex3 iterates m_k = 0.8^k (9, (-1)^k)", np.max(abs(m_hist - m_exact)), 0.0, 1e-12, rel=False)
check("Ex3 J_k = 45 * 0.64^k", np.max(abs(J_hist - 45 * 0.64 ** ks) / (45 * 0.64 ** ks)), 0.0, 1e-12, rel=False)
check("Ex3 J_0 = 45", J_hist[0], 45.0, 1e-15)
g0 = H @ m_hist[0]
check("Ex3 first step length 0.2", (g0 @ g0) / (g0 @ H @ g0), 0.2, 1e-15)
# H-norm of error: ||e_k||_H^2 = 2 J_k = 90 * 0.64^k
eH2 = np.einsum("ki,ij,kj->k", m_hist, H, m_hist)
check("Ex3 ||e_k||_H^2 = 2 J_k", np.max(abs(eH2 - 2 * J_hist)), 0.0, 1e-12, rel=False)
# iterations to reduce J by 10^6
n_iter = 6 / np.log10(1 / 0.64)
print(f"  iterations for 10^-6 reduction, kappa=9:   {n_iter:.2f} -> {int(np.ceil(n_iter))}")
check("Ex3 iterations for 10^-6 (kappa=9) ~ 31", int(np.ceil(n_iter)), 31, 0, rel=False)
check("Ex3 J_31/J_0 < 1e-6 and J_30/J_0 > 1e-6", float(J_hist[31] / J_hist[0] < 1e-6 < J_hist[30] / J_hist[0]), 1.0, 0, rel=False)
rate100 = (99.0 / 101.0) ** 2
n100 = 6 / np.log10(1 / rate100)
print(f"  kappa=100: rate={rate100:.5f}, iterations {n100:.1f} -> {int(np.ceil(n100))}")
check("Ex3 iterations for 10^-6 (kappa=100) = 346", int(np.ceil(n100)), 346, 0, rel=False)
r4 = ((1e4 - 1) / (1e4 + 1)) ** 2
n4 = 6 / np.log10(1 / r4)
print(f"  kappa=1e4: iterations {n4:.0f}")
check("Ex3 iterations for 10^-6 (kappa=1e4) ~ 3.5e4", n4, 3.5e4, 2e-2)
check("Ex3 J_31/J_0 = 0.64^31 ~ 9.8e-7", 0.64 ** 31, 9.8e-7, 1e-2)
print(f"  0.64^31 = {0.64**31:.3e}, 6/log10(1/0.64) = {6/np.log10(1/0.64):.3f}")
# a start along an eigenvector converges in one step
m1 = np.array([0.0, 3.0]); g = H @ m1; m1 = m1 - (g @ g) / (g @ H @ g) * g
check("Ex3 eigenvector start converges in one step", np.linalg.norm(m1), 0.0, 1e-15, rel=False)
# Kantorovich bound for random SPD matrices and random g
rng = np.random.default_rng(3)
worst = 0.0
for _ in range(2000):
    Q, _ = np.linalg.qr(rng.standard_normal((4, 4)))
    mu = rng.uniform(0.5, 20.0, 4)
    Hr = Q @ np.diag(mu) @ Q.T
    g = rng.standard_normal(4)
    ratio = 1 - (g @ g) ** 2 / ((g @ Hr @ g) * (g @ np.linalg.solve(Hr, g)))
    bound = ((mu.max() / mu.min() - 1) / (mu.max() / mu.min() + 1)) ** 2
    worst = max(worst, ratio - bound)
check("Ex3 Kantorovich contraction bound never violated (random tests)", float(worst <= 1e-12), 1.0, 0, rel=False)

print()
print("=" * 72)
print("Example 4: cross-correlation delay time")
print("=" * 72)

t0, sig = 5.0, 0.6
Tw = 12.0  # the pulse is negligible outside [0, Tw]


def s_fun(t, eps=0.0):
    """Synthetic (windowed) waveform; eps multiplies a perturbation eta(t)."""
    return (t - t0) * np.exp(-((t - t0) / sig) ** 2) + eps * eta_fun(t)


def eta_fun(t):
    return np.cos(2.0 * t) * np.exp(-((t - t0 - 0.5) / 1.0) ** 2)


def sdot_fun(t, eps=0.0):
    base = np.exp(-((t - t0) / sig) ** 2) * (1 - 2 * (t - t0) ** 2 / sig ** 2)
    detadt = np.exp(-((t - t0 - 0.5) / 1.0) ** 2) * (-2 * np.sin(2 * t) - 2 * (t - t0 - 0.5) * np.cos(2 * t))
    return base + eps * detadt


def sddot_fun(t):
    u = (t - t0) / sig
    # d/dt [ e^{-u^2}(1-2u^2) ] with du/dt = 1/sig
    return np.exp(-u ** 2) * (-2 * u * (1 - 2 * u ** 2) - 4 * u) / sig


# check derivatives by finite difference
tt = np.linspace(0, Tw, 7)
hh = 1e-5
check("Ex4 sdot analytic vs FD", np.max(abs(sdot_fun(tt, 0.3) - (s_fun(tt + hh, 0.3) - s_fun(tt - hh, 0.3)) / (2 * hh))), 0.0, 1e-8, rel=False)
check("Ex4 sddot analytic vs FD", np.max(abs(sddot_fun(tt) - (sdot_fun(tt + hh) - sdot_fun(tt - hh)) / (2 * hh))), 0.0, 1e-7, rel=False)

Delta = 0.4       # observed arrives LATER than synthetic by Delta
a_obs, sig_obs = 1.0, sig    # set per case below


def s_obs(t):
    return a_obs * (t - Delta - t0) * np.exp(-((t - Delta - t0) / sig_obs) ** 2)


def s_obs_dot(t):
    u = (t - Delta - t0) / sig_obs
    return a_obs * np.exp(-u ** 2) * (1 - 2 * u ** 2)


def Cprime(tau, eps=0.0):
    # C'(tau) = int s_obs(t - tau) sdot(t) dt  (lecture convention C = int s_obs(t-tau) s(t) dt)
    return quad(lambda t: s_obs(t - tau) * sdot_fun(t, eps), 0.0, Tw, epsabs=1e-14, epsrel=1e-13, limit=400)[0]


def C(tau, eps=0.0):
    return quad(lambda t: s_obs(t - tau) * s_fun(t, eps), 0.0, Tw, epsabs=1e-14, epsrel=1e-13, limit=400)[0]


def taubar(eps=0.0):
    # locate the global maximum of C on a grid, then refine the root of C' by bisection
    grid = np.linspace(-1.5, 0.5, 81)
    Cg = np.array([C(tau, eps) for tau in grid])
    tau_star = grid[np.argmax(Cg)]
    return brentq(lambda tau: Cprime(tau, eps), tau_star - 0.05, tau_star + 0.05, xtol=1e-14, rtol=1e-14)


def run_case(a, sg, label):
    """Compare d taubar/d eps from FD, the exact implicit-differentiation formula and
    the approximate formula (s_obs(t - taubar) ~ s(t))."""
    global a_obs, sig_obs
    a_obs, sig_obs = a, sg
    tb = taubar()
    q = dict(epsabs=1e-14, epsrel=1e-13, limit=400)
    num_exact = quad(lambda t: s_obs_dot(t - tb) * eta_fun(t), 0.0, Tw, **q)[0]
    den_exact = quad(lambda t: s_obs(t - tb) * sddot_fun(t), 0.0, Tw, **q)[0]
    den_exact2 = -quad(lambda t: s_obs_dot(t - tb) * sdot_fun(t), 0.0, Tw, **q)[0]
    check(f"Ex4 [{label}] C''(taubar) two forms agree", den_exact, den_exact2, 1e-9)
    check(f"Ex4 [{label}] C''(taubar) < 0 (maximum)", float(den_exact < 0), 1.0, 0, rel=False)
    dtau_exact = num_exact / den_exact
    num_apx = quad(lambda t: sdot_fun(t) * eta_fun(t), 0.0, Tw, **q)[0]
    den_apx = quad(lambda t: sdot_fun(t) ** 2, 0.0, Tw, **q)[0]
    dtau_apx = -num_apx / den_apx
    hh = 1e-4
    dtau_fd = (taubar(hh) - taubar(-hh)) / (2 * hh)
    print(f"  [{label}] taubar = {tb:.10f};  d taubar/d eps: FD = {dtau_fd:.8f}, exact = {dtau_exact:.8f}, "
          f"approx = {dtau_apx:.8f} (ratio approx/FD = {dtau_apx/dtau_fd:.4f})")
    print(f"  [{label}] C''(taubar) = {den_exact:.8f},  -int sdot^2 dt = {-den_apx:.8f}")
    check(f"Ex4 [{label}] exact formula vs FD", dtau_exact, dtau_fd, 1e-6)
    return tb, dtau_fd, dtau_exact, dtau_apx


# (i) identical shapes: taubar = -Delta exactly and both formulae agree
tb, fd_, ex_, ap_ = run_case(1.0, sig, "identical")
check("Ex4 identical waveforms: taubar = -Delta (lecture sign convention)", tb, -Delta, 1e-10)
check("Ex4 identical waveforms: approximate formula exact", ap_, fd_, 1e-6)
# (ii) amplitude mismatch only: the amplitude cancels, approximate formula still exact
tb, fd_, ex_, ap_ = run_case(1.2, sig, "amplitude x1.2")
check("Ex4 amplitude mismatch: taubar = -Delta", tb, -Delta, 1e-10)
check("Ex4 amplitude mismatch: approximate formula exact", ap_, fd_, 1e-6)
# (iii) width mismatch (sigma_obs = 1.1 sigma): approximate formula is only approximate
tb, fd_, ex_, ap_ = run_case(1.2, 1.1 * sig, "width x1.1")
check("Ex4 width mismatch: taubar still = -Delta (symmetric pulses)", tb, -Delta, 1e-9)
check("Ex4 width mismatch: approximate formula within 10% of FD", abs(ap_ / fd_ - 1) < 0.10, True, 0, rel=False)
check("Ex4 width mismatch: approximate formula off by 8.6%", ap_ / fd_, 0.914, 2e-3)
tb, fd_, ex_, ap_ = run_case(1.0, 1.02 * sig, "width x1.02")
check("Ex4 mild width mismatch: approximate formula within 2% of FD", abs(ap_ / fd_ - 1) < 0.02, True, 0, rel=False)

# (iv) pure delay perturbation eta = -sdot gives d taubar/d eps = +1 in the approximate formula
a_obs, sig_obs = 1.0, sig
den_apx = quad(lambda t: sdot_fun(t) ** 2, 0.0, Tw, epsabs=1e-14, epsrel=1e-13)[0]
check("Ex4 pure-delay perturbation gives +1 (approximate formula)",
      -quad(lambda t: sdot_fun(t) * (-sdot_fun(t)), 0, Tw)[0] / den_apx, 1.0, 1e-12)
tb0 = taubar()
val = quad(lambda t: s_obs_dot(t - tb0) * (-sdot_fun(t)), 0, Tw, epsabs=1e-14, epsrel=1e-13)[0] / \
      quad(lambda t: s_obs(t - tb0) * sddot_fun(t), 0, Tw, epsabs=1e-14, epsrel=1e-13)[0]
check("Ex4 pure-delay perturbation gives +1 (exact formula, identical shapes)", val, 1.0, 1e-9)

print()
print("=" * 72)
print("Example 5: cost accounting")
print("=" * 72)

c_min = 10.0          # minutes per simulation
n_ev = 1000
n_it = 15
minutes_per_year = 365.25 * 24 * 60
adj2 = 2 * n_ev * n_it * c_min
adj3 = 3 * n_ev * n_it * c_min
fd4 = (10 ** 4 + 1) * n_ev * n_it * c_min
fd6 = (10 ** 6 + 1) * n_ev * n_it * c_min
print(f"  adjoint (2 sims): {adj2:.3e} min = {adj2/60:.0f} h = {adj2/1440:.1f} days = {adj2/minutes_per_year:.3f} yr")
print(f"  adjoint (3 sims): {adj3:.3e} min = {adj3/60:.0f} h = {adj3/1440:.1f} days = {adj3/minutes_per_year:.3f} yr")
print(f"  FD, m=1e4: {fd4:.4e} min = {fd4/minutes_per_year:.0f} yr")
print(f"  FD, m=1e6: {fd6:.4e} min = {fd6/minutes_per_year:.0f} yr")
check("Ex5 adjoint 2 sims = 3e5 min", adj2, 3e5, 1e-15)
check("Ex5 adjoint 2 sims = 5000 h", adj2 / 60, 5000.0, 1e-15)
check("Ex5 adjoint 2 sims ~ 208 days", adj2 / 1440, 208.3, 1e-3)
check("Ex5 adjoint 3 sims ~ 312 days", adj3 / 1440, 312.5, 1e-3)
check("Ex5 FD m=1e4 ~ 2850 yr", fd4 / minutes_per_year, 2854.0, 2e-3)
check("Ex5 FD m=1e6 ~ 2.85e5 yr", fd6 / minutes_per_year, 2.854e5, 2e-3)
check("Ex5 ratio 3/(m+1), m=1e4", adj3 / fd4, 3 / (10 ** 4 + 1), 1e-15)
# per-iteration wall clock for one event
check("Ex5 per event per iteration (3 sims + 1 trial step) = 40 min", 4 * c_min, 40.0, 1e-15)
check("Ex5 4 sims total ~ 420 days", 4 * n_ev * n_it * c_min / 1440, 420.0, 1e-2)
check("Ex5 100 concurrent simulations -> a few days", 3 * n_ev * n_it * c_min / 1440 / 100, 3.1, 1e-2)
# storage estimate: 1e9 points x 3 comps x 1e4 steps x 4 bytes
check("Ex5 storage 1e9*3*1e4*4 B = 1.2e14 B", 1e9 * 3 * 1e4 * 4, 1.2e14, 1e-15)

print()
print("=" * 72)
print("ALL CHECKS PASSED" if FAILS == 0 else f"{FAILS} CHECK(S) FAILED")
print("=" * 72)
