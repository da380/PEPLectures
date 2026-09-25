"""Numerical checks for examples/examples11.tex (Lecture 22: Free oscillations).

Every numerical or closed-form claim in the solutions is checked here.
Prints OK/FAIL lines; exits with status 1 if anything fails.
"""
import sys
import numpy as np
from scipy import linalg, integrate, optimize

FAILS = 0


def check(name, ok, detail=""):
    global FAILS
    tag = "OK  " if ok else "FAIL"
    if not ok:
        FAILS += 1
    print(f"{tag} {name}  {detail}")


# ----------------------------------------------------------------------
# Example 1: two-degree-of-freedom analogue
# ----------------------------------------------------------------------
print("\n--- Example 1: two-degree-of-freedom system ---")
m, k, F0 = 1.3, 2.7, 0.8            # arbitrary sample parameters
M = np.diag([m, 2 * m])
K = k * np.array([[3.0, -2.0], [-2.0, 4.0]])

w2, S = linalg.eigh(K, M)           # generalised problem, M-orthonormal columns
check("eigenvalues k/m, 4k/m", np.allclose(np.sort(w2), [k / m, 4 * k / m]),
      f"omega^2 = {w2}, expected {[k/m, 4*k/m]}")

s1 = np.array([1.0, 1.0]) / np.sqrt(3 * m)
s2 = np.array([2.0, -1.0]) / np.sqrt(6 * m)
check("s1 = (1,1)/sqrt(3m) satisfies K s = (k/m) M s",
      np.allclose(K @ s1, (k / m) * M @ s1))
check("s2 = (2,-1)/sqrt(6m) satisfies K s = (4k/m) M s",
      np.allclose(K @ s2, (4 * k / m) * M @ s2))
check("M-orthogonality s1^T M s2 = 0", abs(s1 @ M @ s2) < 1e-14, f"{s1 @ M @ s2:.2e}")
check("normalisation s1^T M s1 = 1, s2^T M s2 = 1",
      np.allclose([s1 @ M @ s1, s2 @ M @ s2], 1.0))

# closed-form step response vs direct integration of the ODE
w1 = np.sqrt(k / m)


def u_closed(t):
    u1 = F0 / (3 * k) * (1 - np.cos(w1 * t)) - F0 / (12 * k) * (1 - np.cos(2 * w1 * t))
    u2 = F0 / (3 * k) * (1 - np.cos(w1 * t)) + F0 / (24 * k) * (1 - np.cos(2 * w1 * t))
    return np.array([u1, u2])


Minv = np.linalg.inv(M)


def rhs(t, y):
    u, v = y[:2], y[2:]
    return np.concatenate([v, Minv @ (np.array([0.0, F0]) - K @ u)])


tmax = 6 * 2 * np.pi / w1
sol = integrate.solve_ivp(rhs, (0, tmax), np.zeros(4), rtol=1e-11, atol=1e-13,
                          dense_output=True)
tt = np.linspace(0, tmax, 400)
err = np.max(np.abs(sol.sol(tt)[:2] - u_closed(tt)))
check("closed-form step response matches ODE integration", err < 1e-8, f"max err {err:.2e}")

ustat = np.linalg.solve(K, np.array([0.0, F0]))
check("static limit F0/4k, 3F0/8k", np.allclose(ustat, [F0 / (4 * k), 3 * F0 / (8 * k)]),
      f"{ustat} vs {[F0/(4*k), 3*F0/(8*k)]}")
tper = np.linspace(0, 2 * np.pi / w1, 20001)
tavg = integrate.trapezoid(u_closed(tper), x=tper, axis=1) / (2 * np.pi / w1)
check("time average of response over one period equals static solution",
      np.allclose(tavg, ustat, atol=1e-9), f"{tavg}")

# Fourier transform of (1-cos wk t)/wk^2 with the lecture's convention,
# evaluated at complex omega with negative imaginary part (convergence).
wk = 1.7
om = 0.9 - 0.25j
tgrid = np.linspace(0, 200.0, 400001)
integrand = (1 - np.cos(wk * tgrid)) / wk ** 2 * np.exp(-1j * om * tgrid)
val = integrate.simpson(integrand, x=tgrid)
pred = 1.0 / (1j * om * (wk ** 2 - om ** 2))
check("FT[(1-cos wk t)/wk^2] = 1/(i w (wk^2-w^2))", abs(val - pred) < 1e-6,
      f"{val:.8f} vs {pred:.8f}")
val_H = integrate.simpson(np.exp(-1j * om * tgrid), x=tgrid)
check("FT[H(t)] = 1/(i w)", abs(val_H - 1 / (1j * om)) < 1e-6, f"{val_H:.8f}")

# ----------------------------------------------------------------------
# Example 2: free-free rod with a point stress glut
# ----------------------------------------------------------------------
print("\n--- Example 2: free-free rod ---")
rho, E, L, M0, xs = 2.0, 5.0, 3.0, 0.7, 1.1
c = np.sqrt(E / rho)


def ket(kk, x):
    if kk == 0:
        return np.full_like(np.asarray(x, dtype=float), 1 / np.sqrt(rho * L))
    return np.sqrt(2 / (rho * L)) * np.cos(kk * np.pi * x / L)


xq = np.linspace(0, L, 200001)
G = np.array([[integrate.simpson(rho * ket(a, xq) * ket(b, xq), x=xq) for b in range(5)]
              for a in range(5)])
check("<k|P|k'> = delta_kk' for k,k'<5", np.allclose(G, np.eye(5), atol=1e-10),
      f"max dev {np.max(np.abs(G-np.eye(5))):.1e}")

# cosine coefficients of the static step g(x) = (M0/E)[H(x-xs) - (L-xs)/L]
g = (M0 / E) * ((xq > xs).astype(float) - (L - xs) / L)
gfun = lambda x: (M0 / E) * (float(x > xs) - (L - xs) / L)
a0 = (1 / L) * integrate.quad(gfun, 0, L, points=[xs])[0]
ak = np.array([(2 / L) * integrate.simpson(g * np.cos(kk * np.pi * xq / L), x=xq)
               for kk in range(1, 9)])
ak_pred = np.array([-(2 * M0 / (np.pi * E * kk)) * np.sin(kk * np.pi * xs / L)
                    for kk in range(1, 9)])
check("static step has zero mean (k=0 coefficient)", abs(a0) < 1e-8, f"{a0:.1e}")
check("cosine coefficients of step = -(2 M0/(pi E k)) sin(k pi xs/L)",
      np.allclose(ak, ak_pred, atol=1e-6), f"max dev {np.max(np.abs(ak-ak_pred)):.1e}")


def u_modesum(x, t, N=4000, force=None):
    """Mode sum for the point source (force=None) or a smoothed source."""
    tot = np.zeros_like(np.asarray(x, dtype=float))
    for kk in range(1, N + 1):
        wk_ = kk * np.pi * c / L
        if force is None:
            fk = -M0 * np.sqrt(2 / (rho * L)) * (kk * np.pi / L) * np.sin(kk * np.pi * xs / L)
        else:
            fk = force(kk)
        tot += fk * (1 - np.cos(wk_ * t)) / wk_ ** 2 * ket(kk, x)
    return tot


# (i) time-averaged point-source response -> static step
xpts = np.array([0.3, 0.8, 1.6, 2.5])
step = (M0 / E) * ((xpts > xs).astype(float) - (L - xs) / L)
avg = np.zeros_like(xpts)
for kk in range(1, 20001):
    wk_ = kk * np.pi * c / L
    fk = -M0 * np.sqrt(2 / (rho * L)) * (kk * np.pi / L) * np.sin(kk * np.pi * xs / L)
    avg += fk / wk_ ** 2 * ket(kk, xpts)
check("time-averaged mode sum -> step of size M0/E across xs (20000 terms)",
      np.allclose(avg, step, atol=2e-4), f"{avg} vs {step}")

# (ii) mode sum vs finite-difference time stepping for a Gaussian-smoothed glut
sig = 0.06


def Sbar(x):
    return M0 * np.exp(-(x - xs) ** 2 / (2 * sig ** 2)) / (np.sqrt(2 * np.pi) * sig)


def force_smooth(kk):
    dket = -np.sqrt(2 / (rho * L)) * (kk * np.pi / L) * np.sin(kk * np.pi * xq / L)
    return integrate.simpson(Sbar(xq) * dket, x=xq)


def fd_rod(Nn):
    """Second-order finite-difference solution of the rod with the smoothed glut."""
    h = L / Nn
    xn = np.linspace(0, L, Nn + 1)
    xh = 0.5 * (xn[1:] + xn[:-1])
    Sh = Sbar(xh)
    mass = np.full(Nn + 1, rho * h)
    mass[0] = mass[-1] = rho * h / 2

    def rhs_rod(t, y):
        u, v = y[:Nn + 1], y[Nn + 1:]
        sigma = E * np.diff(u) / h - Sh                 # E u' - Sbar at half nodes
        acc = np.zeros(Nn + 1)
        acc[:-1] += sigma                               # m_j u_j'' = sigma_{j+1/2} - sigma_{j-1/2}
        acc[1:] -= sigma                                # (sigma = 0 outside the rod)
        return np.concatenate([v, acc / mass])

    T_end = 2.3 * L / c
    solr = integrate.solve_ivp(rhs_rod, (0, T_end), np.zeros(2 * (Nn + 1)), rtol=1e-9,
                               atol=1e-12, t_eval=[0.7 * L / c, T_end], method="DOP853")
    return xn, solr.t, solr.y[:Nn + 1, :]


errs = []
for Nn in (600, 1200):
    xn, tvals, ufd = fd_rod(Nn)
    ums = np.array([u_modesum(xn, tv, N=300, force=force_smooth) for tv in tvals]).T
    scale = np.max(np.abs(ufd))
    errs.append(np.max(np.abs(ums - ufd)) / scale)
check("mode sum (smoothed source) is the limit of the finite-difference rod solution",
      errs[1] < 2.5e-3 and errs[0] / errs[1] > 3.0,
      f"rel err {errs[0]:.2e} (N=600) -> {errs[1]:.2e} (N=1200), ratio {errs[0]/errs[1]:.2f}")
mom = [integrate.simpson(rho * ufd[:, i], x=xn) for i in range(ufd.shape[1])]
check("centre of mass does not move: int rho u dx = 0", max(abs(v) for v in mom) < 1e-6 * scale,
      f"{mom}")

# (iii) cosine-series (completeness) representation of delta'(x-xs), tested weakly
phi = np.exp(-(xq - 1.9) ** 2 / (2 * 0.2 ** 2)) * (xq - 1.2) ** 2      # smooth test fn
dphi = np.gradient(phi, xq)
lhs = -dphi[np.argmin(np.abs(xq - xs))]                                  # int delta' phi = -phi'(xs)
rhs_series = sum((2 * np.pi / L ** 2) * kk * np.sin(kk * np.pi * xs / L)
                 * integrate.simpson(np.cos(kk * np.pi * xq / L) * phi, x=xq) for kk in range(1, 3000))
check("delta'(x-xs) = (2 pi/L^2) sum k sin(k pi xs/L) cos(k pi x/L) (weakly)",
      abs(lhs - rhs_series) < 1e-4 * max(1, abs(lhs)), f"{lhs:.6f} vs {rhs_series:.6f}")

# ----------------------------------------------------------------------
# Example 3: Duhamel mode sum against direct integration (2-dof system)
# ----------------------------------------------------------------------
print("\n--- Example 3: Duhamel form of the mode sum ---")


def fvec(t):
    return np.array([np.sin(0.7 * t), t * np.exp(-0.4 * t)])


def rhs_f(t, y):
    u, v = y[:2], y[2:]
    return np.concatenate([v, Minv @ (fvec(t) - K @ u)])


tmax = 25.0
solf = integrate.solve_ivp(rhs_f, (0, tmax), np.zeros(4), rtol=1e-11, atol=1e-13,
                           dense_output=True)
tchk = np.array([3.1, 9.4, 17.7, 24.0])
u_duh = np.zeros((2, len(tchk)))
for j, tv in enumerate(tchk):
    for s_, w2_ in ((s1, k / m), (s2, 4 * k / m)):
        wk_ = np.sqrt(w2_)
        q = integrate.quad(lambda tp: np.sin(wk_ * (tv - tp)) / wk_ * (s_ @ fvec(tp)), 0, tv,
                           limit=400, epsabs=1e-12, epsrel=1e-12)[0]
        u_duh[:, j] += q * s_
errd = np.max(np.abs(u_duh - solf.sol(tchk)[:2]))
check("Duhamel mode sum = ODE solution for a general force", errd < 1e-8, f"max err {errd:.2e}")

# transform of sin(wk t)/wk
integrand = np.sin(wk * tgrid) / wk * np.exp(-1j * om * tgrid)
val = integrate.simpson(integrand, x=tgrid)
check("FT[sin(wk t)/wk] = 1/(wk^2 - w^2)", abs(val - 1 / (wk ** 2 - om ** 2)) < 1e-6,
      f"{val:.8f} vs {1/(wk**2-om**2):.8f}")

# ----------------------------------------------------------------------
# Example 4: spectral resolution
# ----------------------------------------------------------------------
print("\n--- Example 4: spectral resolution ---")
T = 1.0
ff = np.linspace(0, 3, 3000001)
x = np.pi * ff * T
amp = np.abs(np.sinc(x / np.pi))          # |h~|/T as a function of f
pw = amp ** 2
f_half_amp = ff[np.argmax(amp < 0.5)]
f_half_pw = ff[np.argmax(pw < 0.5)]
check("FWHM of |h~|^2 = 0.886/T", abs(2 * f_half_pw - 0.886) < 1e-3, f"{2*f_half_pw:.4f}/T")
check("FWHM of |h~| = 1.207/T", abs(2 * f_half_amp - 1.207) < 1e-3, f"{2*f_half_amp:.4f}/T")
xr = optimize.brentq(lambda z: np.sin(z) / z - 1 / np.sqrt(2), 0.5, 2.5)
check("sin x/x = 1/sqrt2 at x = 1.3916", abs(xr - 1.3916) < 1e-4, f"{xr:.5f}, 2x/pi = {2*xr/np.pi:.4f}")
xa = optimize.brentq(lambda z: np.sin(z) / z - 0.5, 0.5, 2.5)
check("sin x/x = 1/2 at x = 1.8955", abs(xa - 1.8955) < 1e-4, f"{xa:.5f}, 2x/pi = {2*xa/np.pi:.4f}")
first_zero = ff[1:][np.argmax(amp[1:] < 1e-6)]
check("first zero of |h~| at f = 1/T", abs(first_zero - 1.0) < 1e-5, f"{first_zero:.6f}")

# windowed cosine: analytic form vs direct quadrature
A_, w0, Tw = 1.3, 2.0 * np.pi * 0.37, 41.0


def htil(w, T_):
    return (1 - np.exp(-1j * w * T_)) / (1j * w)


for wtest in (0.9 * w0, w0 + 0.05, 1.2 * w0):
    tg = np.linspace(0, Tw, 200001)
    direct = integrate.simpson(A_ * np.cos(w0 * tg) * np.exp(-1j * wtest * tg), x=tg)
    pred = 0.5 * A_ * (htil(wtest - w0, Tw) + htil(wtest + w0, Tw))
    check(f"windowed cosine transform at w={wtest:.3f}", abs(direct - pred) < 1e-7,
          f"{direct:.6f} vs {pred:.6f}")

# 0S2 numbers
f0, dfs, Q = 0.3094e-3, 4.6e-6, 510.0
Tneed = 1 / dfs
tau = Q / (np.pi * f0)
check("record length 1/(4.6 microHz) = 60 h", abs(Tneed / 3600 - 60.4) < 0.5, f"{Tneed/3600:.1f} h")
check("1/(24 h) = 11.6 microHz", abs(1e6 / (24 * 3600) - 11.6) < 0.05, f"{1e6/86400:.2f} microHz")
check("1/(120 h) = 2.3 microHz", abs(1e6 / (120 * 3600) - 2.31) < 0.01, f"{1e6/432000:.2f} microHz")
check("tau = Q/(pi f0) = 6.1 days", abs(tau / 86400 - 6.07) < 0.03, f"{tau/86400:.2f} days")
check("Lorentzian FWHM f0/Q = 0.61 microHz", abs(f0 / Q * 1e6 - 0.607) < 0.005, f"{f0/Q*1e6:.3f}")
check("decay over 120 h: exp(-T/tau) = 0.44", abs(np.exp(-120 * 3600 / tau) - 0.44) < 0.01,
      f"{np.exp(-120*3600/tau):.3f}")
check("omega0 T for T = 24 h is about 170", abs(2 * np.pi * f0 * 86400 - 168) < 1)

# resolution in the figure: count maxima of the two-cosine spectrum near f0
f1_, f2_ = f0 - dfs / 2, f0 + dfs / 2
fg = np.linspace(0.300e-3, 0.319e-3, 40001)


def spec(f, T_):
    w = 2 * np.pi * f
    tot = 0
    for fk in (f1_, f2_):
        wk_ = 2 * np.pi * fk
        tot = tot + 0.5 * (htil(w - wk_, T_) + htil(w + wk_, T_))
    return np.abs(tot)


for Th, expect in ((24.0, 1), (120.0, 2)):
    s = spec(fg, Th * 3600)
    s = s / s.max()
    peaks = np.where((s[1:-1] > s[:-2]) & (s[1:-1] > s[2:]) & (s[1:-1] > 0.5))[0]
    check(f"T = {Th:.0f} h: {expect} peak(s) above half maximum in [0.300,0.319] mHz",
          len(peaks) == expect, f"found {len(peaks)} at {fg[peaks+1]*1e3}")

# ----------------------------------------------------------------------
# Example 5: pointwise stability of an isotropic solid
# ----------------------------------------------------------------------
print("\n--- Example 5: pointwise stability ---")
rng = np.random.default_rng(1)
d = np.eye(3)


def A_iso(lam, mu):
    return (lam * np.einsum("ij,kl->ijkl", d, d)
            + mu * (np.einsum("ik,jl->ijkl", d, d) + np.einsum("il,jk->ijkl", d, d)))


lam, mu = 1.7, 0.9
A = A_iso(lam, mu)
kap = lam + 2 * mu / 3
e = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
e = 0.5 * (e + e.T)
tr = np.trace(e)
dev = e - tr / 3 * d
lhs = np.einsum("ijkl,ij,kl->", A, e.conj(), e)
rhs = kap * abs(tr) ** 2 + 2 * mu * np.sum(dev.conj() * dev)
check("A e* e = kappa |tr e|^2 + 2 mu d*:d (random complex symmetric e)",
      abs(lhs - rhs) < 1e-12, f"{lhs:.6f} vs {rhs:.6f}")
check("e*:e = |tr e|^2/3 + d*:d", abs(np.sum(e.conj() * e) - (abs(tr) ** 2 / 3 + np.sum(dev.conj() * dev))) < 1e-12)

# eigenvalues of A on symmetric tensors: 3 kappa (x1), 2 mu (x5)
basis = []
for i in range(3):
    b = np.zeros((3, 3)); b[i, i] = 1; basis.append(b)
for (i, j) in ((0, 1), (0, 2), (1, 2)):
    b = np.zeros((3, 3)); b[i, j] = b[j, i] = 1 / np.sqrt(2); basis.append(b)
Asym = np.array([[np.einsum("ijkl,ij,kl->", A, p, q) for q in basis] for p in basis])
ev = np.sort(np.linalg.eigvalsh(Asym))
check("eigenvalues of A on symmetric tensors = {3kappa, 2mu x5}",
      np.allclose(ev, np.sort([3 * kap] + [2 * mu] * 5)), f"{ev}")
check("largest c0 = min(3 kappa, 2 mu)", abs(ev[0] - min(3 * kap, 2 * mu)) < 1e-12)

# counter-example lambda = -0.9 mu
mu = 1.0; lam = -0.9 * mu; rho_ = 1.0
kap = lam + 2 * mu / 3
alpha = np.sqrt((lam + 2 * mu) / rho_); beta = np.sqrt(mu / rho_)
check("lambda=-0.9mu: kappa = -7mu/30 < 0", abs(kap + 7 / 30) < 1e-12 and kap < 0, f"kappa/mu = {kap:.4f}")
check("lambda=-0.9mu: alpha = sqrt(1.1) beta = 1.049 beta > beta > 0",
      abs(alpha / beta - np.sqrt(1.1)) < 1e-12 and alpha > beta > 0, f"alpha/beta = {alpha/beta:.4f}")
A = A_iso(lam, mu)
check("uniform dilatation e=delta gives A e e = 9 kappa = -2.1 mu < 0",
      abs(np.einsum("ijkl,ij,kl->", A, d, d) - 9 * kap) < 1e-12 and 9 * kap < 0, f"{9*kap:.3f}")
check("kappa + 4mu/3 = rho alpha^2 > 0", abs(kap + 4 * mu / 3 - rho_ * alpha ** 2) < 1e-12)
check("kappa + mu/3 = rho(alpha^2 - beta^2) > 0", abs(kap + mu / 3 - rho_ * (alpha ** 2 - beta ** 2)) < 1e-12)

# rank-one strains: A e e = mu|a|^2 + (kappa + mu/3)(a.p)^2 = rho a.Gamma a
for _ in range(3):
    a = rng.normal(size=3); p = rng.normal(size=3); p /= np.linalg.norm(p)
    e = 0.5 * (np.outer(a, p) + np.outer(p, a))
    lhs = np.einsum("ijkl,ij,kl->", A, e, e)
    rhs = mu * a @ a + (kap + mu / 3) * (a @ p) ** 2
    Gam = np.einsum("ijkl,j,l->ik", A, p, p) / rho_
    check("rank-one strain: A e e = mu|a|^2+(kappa+mu/3)(a.p)^2 = rho a.Gamma(p) a",
          abs(lhs - rhs) < 1e-12 and abs(lhs - rho_ * a @ Gam @ a) < 1e-12,
          f"{lhs:.6f} {rhs:.6f} {rho_*a@Gam@a:.6f}")
a = p.copy()
e = np.outer(p, p)
check("longitudinal rank-one strain a=p: tr e = 1, d:d = 2/3",
      abs(np.trace(e) - 1) < 1e-12 and abs(np.sum((e - d / 3) ** 2) - 2 / 3) < 1e-12)
check("Christoffel eigenvalues for lambda=-0.9mu are positive: alpha^2, beta^2",
      np.allclose(np.sort(np.linalg.eigvalsh(Gam)), [beta ** 2, beta ** 2, alpha ** 2]))

print("\nALL PASSED" if FAILS == 0 else f"\n{FAILS} CHECK(S) FAILED")
sys.exit(1 if FAILS else 0)
