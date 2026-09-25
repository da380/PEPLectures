"""
Numerical checks for the worked examples accompanying Lecture 15
(Plane wave propagation).  Every numerical or closed-form claim in
examples/examples4.tex is checked here; the script prints OK/FAIL lines.
"""
import numpy as np

np.set_printoptions(precision=6, suppress=True)
status = {"fail": 0}


def check(name, ok, detail=""):
    tag = "OK  " if ok else "FAIL"
    if not ok:
        status["fail"] += 1
    print(f"{tag} {name}  {detail}")


def close(a, b, rtol=1e-9, atol=0.0):
    return np.allclose(a, b, rtol=rtol, atol=atol)


d3 = np.eye(3)


def christoffel(A, rho, phat):
    """Gamma_ik = A_ijkl p_j p_l / rho."""
    return np.einsum("ijkl,j,l->ik", A, phat, phat) / rho


def isotropic(lam, mu):
    return (lam * np.einsum("ij,kl->ijkl", d3, d3)
            + mu * (np.einsum("ik,jl->ijkl", d3, d3) + np.einsum("il,jk->ijkl", d3, d3)))


def cubic(c11, c12, c44):
    A = isotropic(c12, c44)
    dd = c11 - c12 - 2 * c44
    for n in range(3):
        A[n, n, n, n] += dd
    return A


def cubic_perturbation(c11, c12, c44):
    dA = np.zeros((3, 3, 3, 3))
    dd = c11 - c12 - 2 * c44
    for n in range(3):
        dA[n, n, n, n] = dd
    return dA


def transversely_isotropic(lam, mu, gam, xi, zeta, nu):
    A = isotropic(lam, mu)
    A += 8 * gam * np.einsum("i,j,k,l->ijkl", nu, nu, nu, nu)
    A += 4 * xi * (np.einsum("i,j,kl->ijkl", nu, nu, d3) + np.einsum("ij,k,l->ijkl", d3, nu, nu))
    A -= zeta * (np.einsum("i,k,jl->ijkl", nu, nu, d3) + np.einsum("j,k,il->ijkl", nu, nu, d3)
                 + np.einsum("j,l,ik->ijkl", nu, nu, d3) + np.einsum("i,l,jk->ijkl", nu, nu, d3))
    return A


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


# ---------------------------------------------------------------- Example 1
print("\n=== Example 1: PREM-like isotropic values at 100 km depth ===")
rho1, kappa, mu1 = 3.38e3, 1.30e11, 6.8e10
lam1 = kappa - 2 * mu1 / 3
alpha1 = np.sqrt((lam1 + 2 * mu1) / rho1)
beta1 = np.sqrt(mu1 / rho1)
nu_p = lam1 / (2 * (lam1 + mu1))
print(f"     lambda = {lam1:.4e} Pa, alpha = {alpha1:.1f} m/s, beta = {beta1:.1f} m/s,"
      f" alpha/beta = {alpha1/beta1:.4f}, Poisson ratio = {nu_p:.4f}")
check("lambda = 8.47e10 Pa", close(lam1, 8.4667e10, rtol=1e-4), f"{lam1:.4e}")
check("alpha = 8.08 km/s", abs(alpha1 - 8080) < 5, f"{alpha1:.1f}")
check("beta = 4.49 km/s", abs(beta1 - 4485) < 5, f"{beta1:.1f}")
check("alpha/beta = 1.80", abs(alpha1 / beta1 - 1.801) < 2e-3, f"{alpha1/beta1:.4f}")
check("Poisson ratio = 0.277", abs(nu_p - 0.2773) < 5e-4, f"{nu_p:.4f}")
# alpha/beta from Poisson ratio
ratio_from_nu = np.sqrt(2 * (1 - nu_p) / (1 - 2 * nu_p))
check("alpha/beta = sqrt(2(1-nu)/(1-2nu))", close(ratio_from_nu, alpha1 / beta1), f"{ratio_from_nu:.4f}")
# Poisson solid
check("Poisson solid: alpha/beta = sqrt(3), nu = 1/4",
      close(np.sqrt(3.0), np.sqrt(2 * (1 - 0.25) / (1 - 0.5))))
# eigen-decomposition of the isotropic Christoffel matrix in a random direction
A1 = isotropic(lam1, mu1)
rng = np.random.default_rng(1)
ph = unit(rng.normal(size=3))
ev = np.sort(np.linalg.eigvalsh(christoffel(A1, rho1, ph)))
check("isotropic Christoffel eigenvalues = beta^2, beta^2, alpha^2",
      close(ev, [beta1**2, beta1**2, alpha1**2]), f"{np.sqrt(ev)}")
# kappa = lambda + 2mu/3 consistency with the constants quoted
check("kappa/rho and mu/rho used consistently", close(alpha1**2, (kappa + 4 * mu1 / 3) / rho1))

# ---------------------------------------------------------------- Example 2
print("\n=== Example 2: cubic crystal ===")
c11, c12, c44, rho = 320e9, 70e9, 80e9, 3300.0
dd = c11 - c12 - 2 * c44
A2 = cubic(c11, c12, c44)
# symmetries
check("cubic A has minor and major symmetries",
      close(A2, np.transpose(A2, (1, 0, 2, 3))) and close(A2, np.transpose(A2, (0, 1, 3, 2)))
      and close(A2, np.transpose(A2, (2, 3, 0, 1))))
# rho Gamma = (c12 + c44) p p^T + c44 I + d diag(p_n^2)
for _ in range(3):
    ph = unit(rng.normal(size=3))
    G_formula = ((c12 + c44) * np.outer(ph, ph) + c44 * d3 + dd * np.diag(ph**2)) / rho
    check("closed form of cubic Christoffel matrix", close(christoffel(A2, rho, ph), G_formula))


def speeds(A, rho, ph):
    w, v = np.linalg.eigh(christoffel(A, rho, unit(ph)))
    return np.sqrt(w), v


for label, ph, rc2, pols in [
    ("[100]", [1, 0, 0], [c44, c44, c11], None),
    ("[110]", [1, 1, 0], [c44, 0.5 * (c11 - c12), 0.5 * (c11 + c12 + 2 * c44)],
     [[0, 0, 1], [1, -1, 0], [1, 1, 0]]),
    ("[111]", [1, 1, 1], [(c11 - c12 + c44) / 3] * 2 + [(c11 + 2 * c12 + 4 * c44) / 3], None),
]:
    c, v = speeds(A2, rho, ph)
    pred = np.sqrt(np.sort(np.array(rc2)) / rho)
    check(f"{label} phase speeds", close(c, pred), f"{c/1e3} km/s (predicted {pred/1e3})")
    if pols is not None:
        for k, pol in enumerate(pols):
            pol = unit(pol)
            check(f"{label} polarisation {k}: {pol}", abs(abs(pol @ v[:, k]) - 1) < 1e-9)
    else:
        # fastest wave polarised along phat; others orthogonal to it
        check(f"{label} fastest polarised along p", abs(abs(unit(ph) @ v[:, 2]) - 1) < 1e-9)
print(f"     [100]: {np.sqrt(np.array([c11, c44]) / rho) / 1e3} km/s")
print(f"     [110]: {np.sqrt(np.array([0.5*(c11+c12+2*c44), 0.5*(c11-c12), c44]) / rho) / 1e3} km/s")
print(f"     [111]: {np.sqrt(np.array([(c11+2*c12+4*c44)/3, (c11-c12+c44)/3]) / rho) / 1e3} km/s")
check("[111] rho c^2 values are 260 and 110 GPa",
      close([(c11 + 2 * c12 + 4 * c44) / 3, (c11 - c12 + c44) / 3], [260e9, 110e9]))
check("Zener ratio 2c44/(c11-c12) = 0.64", close(2 * c44 / (c11 - c12), 0.64))
# Explicit numbers quoted in the text (three significant figures, km/s)
quoted = {"c11": 9.85, "c44": 4.92, "(c11+c12+2c44)/2": 9.13, "(c11-c12)/2": 6.15,
          "(c11+2c12+4c44)/3": 8.88, "(c11-c12+c44)/3": 5.77}
values = {"c11": c11, "c44": c44, "(c11+c12+2c44)/2": 0.5 * (c11 + c12 + 2 * c44),
          "(c11-c12)/2": 0.5 * (c11 - c12), "(c11+2c12+4c44)/3": (c11 + 2 * c12 + 4 * c44) / 3,
          "(c11-c12+c44)/3": (c11 - c12 + c44) / 3}
for k in quoted:
    v = np.sqrt(values[k] / rho) / 1e3
    check(f"quoted speed sqrt({k}/rho) = {quoted[k]} km/s", abs(v - quoted[k]) < 0.005, f"{v:.4f}")

# ---------------------------------------------------------------- Example 3
print("\n=== Example 3: SH wave in a transversely isotropic medium ===")
lam3, mu3, gam3, xi3, zeta3, rho3 = 80e9, 70e9, 2e9, 1e9, 8e9, 3300.0
nu3 = np.array([0.0, 0.0, 1.0])
A3 = transversely_isotropic(lam3, mu3, gam3, xi3, zeta3, nu3)
check("TI A has minor and major symmetries",
      close(A3, np.transpose(A3, (1, 0, 2, 3))) and close(A3, np.transpose(A3, (0, 1, 3, 2)))
      and close(A3, np.transpose(A3, (2, 3, 0, 1))))
e2 = np.array([0.0, 1.0, 0.0])
ok_vec, ok_val, ok_pd, ok_pert = True, True, True, True
max_err_pert = 0.0
for th in np.linspace(0, np.pi, 181):
    ph = np.array([np.sin(th), 0.0, np.cos(th)])
    G = christoffel(A3, rho3, ph)
    Ge2 = G @ e2
    lam_sh = (mu3 - zeta3 * np.cos(th)**2) / rho3
    ok_vec &= close(Ge2, lam_sh * e2, atol=1e-6 * abs(lam_sh))
    w = np.linalg.eigvalsh(G)
    ok_val &= np.min(np.abs(w - lam_sh)) < 1e-6 * lam_sh
    ok_pd &= np.all(w > 0)
    # first-order degenerate perturbation theory for the SH branch vs exact
    beta3 = np.sqrt(mu3 / rho3)
    dG = christoffel(A3 - isotropic(lam3, mu3), rho3, ph)
    t1 = np.cross(e2, ph)
    M = np.array([[t1 @ dG @ t1, t1 @ dG @ e2], [e2 @ dG @ t1, e2 @ dG @ e2]])
    ok_pert &= abs(M[0, 1]) < 1e-6 * mu3 / rho3   # e2 decouples exactly
    db = M[1, 1] / (2 * beta3)
    max_err_pert = max(max_err_pert, abs(beta3 + db - np.sqrt(lam_sh)))
check("e_2 is an exact eigenvector of Gamma for all theta", ok_vec)
check("its eigenvalue is (mu - zeta cos^2 theta)/rho", ok_val)
check("sample TI medium has positive-definite Gamma for all theta in the (1,3)-plane", ok_pd)
check("SH branch decouples in the 2x2 degenerate problem (off-diagonal = 0)", ok_pert)
print(f"     max |first-order SH speed - exact| over theta = {max_err_pert:.2f} m/s"
      f" (beta = {np.sqrt(mu3/rho3):.1f} m/s; zeta/mu = {zeta3/mu3:.3f})")
check("first-order SH speed within 0.2% of exact for zeta/mu = 0.114",
      max_err_pert / np.sqrt(mu3 / rho3) < 2e-3)
# The SH sheet is the ellipse mu p1^2 + (mu - zeta) p3^2 = rho
ok_ell = True
for th in np.linspace(0, 2 * np.pi, 97):
    ph = np.array([np.sin(th), 0.0, np.cos(th)])
    c = np.sqrt((mu3 - zeta3 * np.cos(th)**2) / rho3)
    p = ph / c
    ok_ell &= abs(mu3 * p[0]**2 + (mu3 - zeta3) * p[2]**2 - rho3) < 1e-6 * rho3
check("SH slowness curve is the ellipse mu p1^2 + (mu - zeta) p3^2 = rho", ok_ell)
# Voigt constants of the sample medium
Avo = lam3 + 2 * mu3
Cvo = lam3 + 2 * mu3 + 8 * gam3 + 8 * xi3 - 4 * zeta3
Fvo = lam3 + 4 * xi3
Lvo = mu3 - zeta3
Nvo = mu3
check("Voigt A = A_1111", close(A3[0, 0, 0, 0], Avo), f"{Avo/1e9} GPa")
check("Voigt C = A_3333", close(A3[2, 2, 2, 2], Cvo), f"{Cvo/1e9} GPa")
check("Voigt F = A_1133", close(A3[0, 0, 2, 2], Fvo), f"{Fvo/1e9} GPa")
check("Voigt L = A_1313", close(A3[0, 2, 0, 2], Lvo), f"{Lvo/1e9} GPa")
check("Voigt N = A_1212", close(A3[0, 1, 0, 1], Nvo), f"{Nvo/1e9} GPa")
print(f"     SH speeds: horizontal sqrt(N/rho) = {np.sqrt(Nvo/rho3):.1f} m/s,"
      f" vertical sqrt(L/rho) = {np.sqrt(Lvo/rho3):.1f} m/s")
check("quoted SH speeds 4.61 and 4.33 km/s",
      abs(np.sqrt(Nvo / rho3) / 1e3 - 4.61) < 0.005 and abs(np.sqrt(Lvo / rho3) / 1e3 - 4.33) < 0.005)
# along the axis both S-waves are degenerate with rho c^2 = mu - zeta
w_axis = np.sort(np.linalg.eigvalsh(christoffel(A3, rho3, nu3)))
check("along the symmetry axis both S-waves have rho c^2 = mu - zeta, P has C",
      close(w_axis * rho3, [Lvo, Lvo, Cvo]))
w_h = np.sort(np.linalg.eigvalsh(christoffel(A3, rho3, np.array([1.0, 0, 0]))))
check("horizontally: rho c^2 = L (SV), N (SH), A (P)", close(w_h * rho3, [Lvo, Nvo, Avo]))

# ---------------------------------------------------------------- Example 4
print("\n=== Example 4: perturbation theory for the cubic medium ===")
A_iso = isotropic(c12, c44)
dA = cubic_perturbation(c11, c12, c44)
check("cubic A = isotropic(lambda=c12, mu=c44) + dA", close(A2, A_iso + dA))
alpha = np.sqrt((c12 + 2 * c44) / rho)
beta = np.sqrt(c44 / rho)
print(f"     isotropic reference: alpha = {alpha:.1f} m/s, beta = {beta:.1f} m/s, d = {dd/1e9} GPa")
check("quoted alpha_0 = 8.35 km/s, beta_0 = 4.92 km/s", abs(alpha / 1e3 - 8.35) < 0.005 and abs(beta / 1e3 - 4.92) < 0.005)

# generic first-order formulae
# rho dGamma = d diag(p_n^2)
for _ in range(3):
    ph = unit(rng.normal(size=3))
    check("rho dGamma = d diag(p_n^2)", close(christoffel(dA, rho, ph), dd * np.diag(ph**2) / rho))
    # first-order qP: rho c^2 = c12 + 2 c44 + d sum p_n^4
    da = (ph @ christoffel(dA, rho, ph) @ ph) / (2 * alpha)
    check("delta alpha = d sum p_n^4 /(2 rho alpha)", close(da, dd * np.sum(ph**4) / (2 * rho * alpha)))


def first_order(s, ph):
    """Return (first-order qP speed, first-order qS speeds (2,), exact speeds sorted (3,))."""
    ph = unit(ph)
    dG = s * christoffel(dA, rho, ph)
    da = (ph @ dG @ ph) / (2 * alpha)
    # orthonormal basis of the plane orthogonal to phat
    t1 = np.cross(ph, [0, 0, 1.0])
    if np.linalg.norm(t1) < 1e-8:
        t1 = np.cross(ph, [1.0, 0, 0])
    t1 = unit(t1)
    t2 = np.cross(ph, t1)
    M = np.array([[t1 @ dG @ t1, t1 @ dG @ t2], [t2 @ dG @ t1, t2 @ dG @ t2]])
    m = np.linalg.eigvalsh(M)
    db = m / (2 * beta)
    exact = np.sqrt(np.linalg.eigvalsh(christoffel(A_iso + s * dA, rho, ph)))
    return alpha + da, beta + db, exact, (ph @ dG @ ph), m


# (a) [110] quasi-P
ph110 = unit([1, 1, 0])
cP1, cS1, ex, pdGp, m = first_order(1.0, ph110)
check("[110]: p.dGamma.p = d/(2 rho)", close(pdGp, dd / (2 * rho)))
check("[110]: delta alpha = d/(4 rho alpha)", close(cP1 - alpha, dd / (4 * rho * alpha)))
check("[110]: exact rho c_qP^2 = rho alpha^2 + d/2 (first order exact for c^2)",
      close(rho * ex[2]**2, rho * alpha**2 + dd / 2))
print(f"     [110] qP: first order {cP1:.1f} m/s, exact {ex[2]:.1f} m/s, error {cP1-ex[2]:.1f} m/s"
      f" ({100*(cP1-ex[2])/ex[2]:.2f} %)")
check("[110] qP quoted: 9.17 (first order) vs 9.13 km/s (exact), error 37 m/s (0.4%)",
      abs(cP1 / 1e3 - 9.17) < 0.005 and abs(ex[2] / 1e3 - 9.13) < 0.005 and abs((cP1 - ex[2]) - 37) < 1)
second = dd**2 / (32 * rho**2 * alpha**3)
print(f"     second-order estimate of the square-root error d^2/(32 rho^2 alpha^3) = {second:.1f} m/s")
check("square-root error agrees with -d^2/(32 rho^2 alpha^3) to ~10%", abs((cP1 - ex[2]) - second) / second < 0.12)
check("expansion parameter d/(2 rho alpha^2) = 0.196", close(dd / (2 * rho * alpha**2), 90 / 460))

# (b) [110] quasi-S
check("[110]: 2x2 matrix eigenvalues are 0 and d/(2 rho)", close(np.sort(m), [0, dd / (2 * rho)]))
t1, t2 = np.array([0, 0, 1.0]), unit([1, -1, 0])
dG110 = christoffel(dA, rho, ph110)
check("[110]: t1.dG.t1 = 0, t2.dG.t2 = d/(2rho), t1.dG.t2 = 0",
      close(t1 @ dG110 @ t1, 0, atol=1e-12) and close(t2 @ dG110 @ t2, dd / (2 * rho))
      and close(t1 @ dG110 @ t2, 0, atol=1e-12))
check("[110]: exact rho c^2 for S-waves = c44 and (c11-c12)/2, first order exact for c^2",
      close(rho * ex[:2]**2, [c44, 0.5 * (c11 - c12)]))
print(f"     [110] qS2: first order {cS1[1]:.1f} m/s, exact {ex[1]:.1f} m/s, error {cS1[1]-ex[1]:.1f} m/s"
      f" ({100*(cS1[1]-ex[1])/ex[1]:.2f} %); qS1: {cS1[0]:.1f} vs {ex[0]:.1f}")
check("[110] qS2 quoted: delta beta = 1.38 km/s, 6.31 (first order) vs 6.15 km/s (exact)",
      abs((cS1[1] - beta) / 1e3 - 1.38) < 0.005 and abs(cS1[1] / 1e3 - 6.31) < 0.005 and abs(ex[1] / 1e3 - 6.15) < 0.005)
check("expansion parameter d/(2 rho beta^2) = 0.5625", close(dd / (2 * rho * beta**2), 0.5625))

# (c) generic direction (2,3,6)/7
phg = np.array([2, 3, 6.0]) / 7
check("(2,3,6)/7 is a unit vector", close(phg @ phg, 1))
check("sum p_n^4 = 1393/2401", close(np.sum(phg**4), 1393 / 2401))
dG1 = christoffel(dA, rho, phg)
t1 = unit(np.cross(phg, [0, 0, 1.0])); t2 = np.cross(phg, t1)
M = np.array([[t1 @ dG1 @ t1, t1 @ dG1 @ t2], [t2 @ dG1 @ t1, t2 @ dG1 @ t2]])
m = np.linalg.eigvalsh(M)
print(f"     generic direction: rho p.dG.p = {rho*(phg@dG1@phg)/1e9:.3f} GPa,"
      f" rho M = {rho*M/1e9} GPa, eigenvalues rho m = {rho*m/1e9} GPa")
check("generic: first-order rho c_qP^2 = 230 + 90*1393/2401 = 282.22 GPa",
      close(rho * alpha**2 + rho * (phg @ dG1 @ phg), 230e9 + 90e9 * 1393 / 2401))
check("generic: rho M has trace d(1 - sum p_n^4) = 37.78 GPa", close(np.trace(rho * M), dd * (1 - np.sum(phg**4))))
check("generic: quoted 2x2 eigenvalues 9.45 and 28.34 GPa", close(rho * m / 1e9, [9.446, 28.338], rtol=2e-4))
print("     first-order vs exact eigenvalues rho c^2 (GPa) in direction (2,3,6)/7:")
print("       s     qS1(1st)  qS1(exact)  err     qS2(1st)  qS2(exact)  err     qP(1st)  qP(exact)  err")
errs = {}
for s in [1.0, 0.5, 0.25]:
    fo = np.array([c44 + s * rho * m[0], c44 + s * rho * m[1], rho * alpha**2 + s * rho * (phg @ dG1 @ phg)]) / 1e9
    exv = rho * np.linalg.eigvalsh(christoffel(A_iso + s * dA, rho, phg)) / 1e9
    errs[s] = fo - exv
    print(f"      {s:4.2f}  {fo[0]:8.2f}  {exv[0]:8.2f}  {fo[0]-exv[0]:6.2f}   {fo[1]:8.2f}  {exv[1]:8.2f}  {fo[1]-exv[1]:6.2f}"
          f"   {fo[2]:8.2f}  {exv[2]:8.2f}  {fo[2]-exv[2]:6.2f}")
    check(f"s={s}: first-order eigenvalues sum to the exact trace", close(np.sum(fo), np.sum(exv)))
check("generic s=1 quoted eigenvalues: 1st order 89.45, 108.34, 282.22; exact 89.18, 105.58, 285.25 GPa",
      close(errs[1.0] + rho * np.linalg.eigvalsh(christoffel(A_iso + dA, rho, phg)) / 1e9, [89.45, 108.34, 282.22], rtol=1e-4)
      and close(rho * np.linalg.eigvalsh(christoffel(A_iso + dA, rho, phg)) / 1e9, [89.18, 105.58, 285.25], rtol=1e-4))
check("generic s=1 quoted errors +0.27, +2.76, -3.03 GPa", close(errs[1.0], [0.27, 2.76, -3.03], atol=0.006))
r = errs[0.5] / errs[0.25]
check("eigenvalue error scales as s^2 (ratio for s=0.5 vs 0.25 within 10% of 4)", np.all(np.abs(r / 4 - 1) < 0.1), f"ratios {r}")
check("qP pushed up, qS pushed down at second order (signs)", errs[1.0][2] < 0 and np.all(errs[1.0][:2] > 0))
# second-order coefficient for c^2 of qP: sum |t.dG.p|^2/(alpha^2-beta^2)
sec = sum((t @ dG1 @ phg)**2 for t in (t1, t2)) / (alpha**2 - beta**2)
c2_exact = lambda s: np.linalg.eigvalsh(christoffel(A_iso + s * dA, rho, phg))[2]
s = 1e-3
num_sec = (c2_exact(s) - alpha**2 - s * (phg @ dG1 @ phg)) / s**2
check("second-order coefficient of c_qP^2 matches Rayleigh-Schroedinger formula", close(num_sec, sec, rtol=1e-3),
      f"{num_sec:.4e} vs {sec:.4e}  (rho * coefficient = {rho*sec/1e9:.2f} GPa)")
check("quoted second-order coefficient 3.60 GPa; at s=1/4 gives 0.225 GPa vs observed 0.22 GPa",
      abs(rho * sec / 1e9 - 3.60) < 0.005 and abs(rho * sec / 1e9 / 16 - 0.225) < 0.001 and abs(errs[0.25][2] + 0.22) < 0.005)
check("quoted error ratios lie between 3.8 and 3.9", np.all((r > 3.8) & (r < 3.9)))
check("third-order square-root term for [110] qP is ~4 m/s (explains 40 vs 37)",
      abs((dd / (4 * rho * alpha))**3 / (2 * alpha**2) - 3.9) < 0.3)
# speeds at s = 1: lecture formula c = alpha + delta alpha etc. against exact
cP, cS, ex, _, _ = first_order(1.0, phg)
print(f"     s=1 speeds (m/s): qP first order {cP:.1f}, exact {ex[2]:.1f}; qS1 {cS[0]:.1f} vs {ex[0]:.1f}; qS2 {cS[1]:.1f} vs {ex[1]:.1f}")
check("generic s=1 quoted speeds: qP 9.30 vs 9.30, qS1 5.21 vs 5.20, qS2 5.80 vs 5.66 km/s",
      abs(cP / 1e3 - 9.30) < 0.005 and abs(ex[2] / 1e3 - 9.30) < 0.005 and abs(cS[0] / 1e3 - 5.21) < 0.005
      and abs(ex[0] / 1e3 - 5.20) < 0.005 and abs(cS[1] / 1e3 - 5.80) < 0.005 and abs(ex[1] / 1e3 - 5.66) < 0.005)
# decomposition of the qP speed error at s=1: second-order eigenvalue shift vs square-root linearisation
sq_err = -(phg @ dG1 @ phg)**2 / (8 * alpha**3)
print(f"     qP speed error budget: +second-order eigenvalue shift/(2 alpha) = {sec/(2*alpha):.1f} m/s,"
      f" square-root term = {sq_err:.1f} m/s, actual (exact - first order) = {ex[2]-cP:.1f} m/s")

# ---------------------------------------------------------------- Example 5
print("\n=== Example 5: energy flux of a plane wave ===")


def flux_check(A, rho, ph, label):
    w, v = np.linalg.eigh(christoffel(A, rho, unit(ph)))
    ok = True
    for k in range(3):
        c = np.sqrt(w[k]); a = v[:, k]; p = unit(ph) / c
        # energy velocity  v_j = A_ijkl a_i a_k p_l / rho  (per unit f'^2, E = rho f'^2)
        vE = np.einsum("ijkl,i,k,l->j", A, a, a, p) / rho
        E_over_fp2 = 0.5 * rho + 0.5 * np.einsum("ijkl,i,j,k,l", A, a, p, a, p)
        ok &= close(E_over_fp2, rho)
        ok &= close(vE @ p, 1.0)
    check(f"{label}: E = rho f'^2 and v_E . p = 1 for all three waves", ok)
    return w, v


flux_check(A2, rho, [2, 3, 6], "cubic, generic direction")
flux_check(A3, rho3, [np.sin(0.7), 0, np.cos(0.7)], "TI, theta = 0.7")
# isotropic P: flux = rho alpha <udot^2> phat ; S: rho beta phat
ph = unit(rng.normal(size=3))
a = ph; p = ph / alpha1
vP = np.einsum("ijkl,i,k,l->j", A1, a, a, p)
check("isotropic P: A_ijkl a_i a_k p_l = rho alpha phat", close(vP, rho1 * alpha1 * ph))
t = unit(np.cross(ph, [0, 0, 1.0])); p = ph / beta1
vS = np.einsum("ijkl,i,k,l->j", A1, t, t, p)
check("isotropic S: A_ijkl a_i a_k p_l = rho beta phat", close(vS, rho1 * beta1 * ph))
# TI SH: energy velocity (mu sin th, 0, (mu - zeta) cos th)/(rho c), not parallel to phat
th = 0.7
ph = np.array([np.sin(th), 0, np.cos(th)])
c = np.sqrt((mu3 - zeta3 * np.cos(th)**2) / rho3)
vE = np.einsum("ijkl,i,k,l->j", A3, e2, e2, ph / c) / rho3
vE_formula = np.array([mu3 * np.sin(th), 0, (mu3 - zeta3) * np.cos(th)]) / (rho3 * c)
check("TI SH energy velocity = (mu sin, 0, (mu-zeta) cos)/(rho c)", close(vE, vE_formula), f"{vE}")
check("TI SH energy velocity . phat = c", close(vE @ ph, c))
ang = np.degrees(np.arccos(vE @ ph / np.linalg.norm(vE)))
print(f"     TI SH at theta = 40.1 deg: angle between energy velocity and phat = {ang:.2f} deg")
# normal to the ellipse mu p1^2 + (mu - zeta) p3^2 = rho is parallel to vE
p = ph / c
normal = np.array([2 * mu3 * p[0], 0, 2 * (mu3 - zeta3) * p[2]])
check("TI SH energy velocity is normal to the slowness curve", close(unit(normal), unit(vE)))
# maximum deviation angle over theta:  tan(psi) with psi = theta - theta_E
ths = np.linspace(0, np.pi / 2, 9001)
angs = []
for th in ths:
    ph = np.array([np.sin(th), 0, np.cos(th)])
    vE = np.array([mu3 * np.sin(th), 0, (mu3 - zeta3) * np.cos(th)])
    angs.append(np.degrees(np.arccos(unit(vE) @ ph)))
print(f"     maximum angle between SH energy velocity and phat = {max(angs):.2f} deg at theta = {np.degrees(ths[int(np.argmax(angs))]):.1f} deg")
check("quoted maximum deviation about 3.5 deg", abs(max(angs) - 3.5) < 0.1)

print()
print("ALL CHECKS PASSED" if status["fail"] == 0 else f"{status['fail']} CHECK(S) FAILED")
