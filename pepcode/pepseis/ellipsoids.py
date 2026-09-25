"""Classical equilibrium figures of a homogeneous, rotating, self-gravitating fluid: the
Maclaurin spheroids and Jacobi ellipsoids (Chandrasekhar, Ellipsoidal Figures of Equilibrium,
1969).  Index symbols A_i = a1 a2 a3 int_0^inf du / ((a_i^2 + u) Delta),
B_ij = a1 a2 a3 int_0^inf u du / ((a_i^2 + u)(a_j^2 + u) Delta), Delta = sqrt(prod (a_k^2 + u)).
together with A_ij = a1 a2 a3 int_0^inf du / ((a_i^2 + u)(a_j^2 + u) Delta).  A Jacobi
ellipsoid with semi-axes a1 >= a2 >= a3 rotating about the a3 axis satisfies
    a1^2 a2^2 A_12 = a3^2 A_3,     Omega^2 / (pi G rho) = 2 B_12,
while a Maclaurin spheroid (a1 = a2) has Omega^2 / (pi G rho) = 2 (A_1 - (a3/a1)^2 A_3).  Lengths are scaled so that
a1 a2 a3 = 1 (i.e. the volume, and hence the mass, is fixed).
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def _delta(u, a):
    return np.sqrt((a[0] ** 2 + u) * (a[1] ** 2 + u) * (a[2] ** 2 + u))


def A(i, a):
    f = lambda u: 1.0 / ((a[i] ** 2 + u) * _delta(u, a))
    return a[0] * a[1] * a[2] * (quad(f, 0, 1.0)[0] + quad(f, 1.0, np.inf)[0])


def A2(i, j, a):
    f = lambda u: 1.0 / ((a[i] ** 2 + u) * (a[j] ** 2 + u) * _delta(u, a))
    return a[0] * a[1] * a[2] * (quad(f, 0, 1.0)[0] + quad(f, 1.0, np.inf)[0])


def B(i, j, a):
    f = lambda u: u / ((a[i] ** 2 + u) * (a[j] ** 2 + u) * _delta(u, a))
    return a[0] * a[1] * a[2] * (quad(f, 0, 1.0)[0] + quad(f, 1.0, np.inf)[0])


def omega2(a):
    """Omega^2 / (pi G rho) for an equilibrium ellipsoid with semi-axes a (rotation about a[2])."""
    if abs(a[0] - a[1]) < 1e-9 * a[0]:
        return 2.0 * (A(0, a) - (a[2] / a[0]) ** 2 * A(2, a))
    return 2.0 * B(0, 1, a)


def angular_momentum(a):
    """L / sqrt(G M^3 abar) with abar = (a1 a2 a3)^{1/3} and rho = M / (4 pi abar^3 / 3)."""
    abar = (a[0] * a[1] * a[2]) ** (1.0 / 3.0)
    om = np.sqrt(omega2(a) * np.pi * 3.0 / (4.0 * np.pi))          # Omega in units sqrt(G M / abar^3)
    return 0.2 * (a[0] ** 2 + a[1] ** 2) / abar ** 2 * om            # L = (2/5) M (a1^2 + a2^2) Omega


def maclaurin(e):
    """Exact Maclaurin relation for meridional eccentricity e."""
    return 2.0 * np.sqrt(1 - e ** 2) / e ** 3 * (3 - 2 * e ** 2) * np.arcsin(e) - 6.0 * (1 - e ** 2) / e ** 2


def maclaurin_axes(e):
    a = (1 - e ** 2) ** (-1.0 / 6.0)          # a1 = a2 = a, a3 = a sqrt(1-e^2), a^2 a3 = 1
    return np.array([a, a, a * np.sqrt(1 - e ** 2)])


def jacobi_axes(ratio):
    """Semi-axes of the Jacobi ellipsoid with a2/a1 = ratio (< 1), unit volume."""
    def resid(c_over_a):
        a = np.array([1.0, ratio, c_over_a]); a = a / (a.prod()) ** (1.0 / 3.0)
        return a[0] ** 2 * a[1] ** 2 * A2(0, 1, a) - a[2] ** 2 * A(2, a)
    c = brentq(resid, 0.05, 0.99)
    a = np.array([1.0, ratio, c]); return a / (a.prod()) ** (1.0 / 3.0)


if __name__ == "__main__":
    for e in [0.5, 0.81267, 0.92995]:
        a = maclaurin_axes(e); print("Maclaurin e=%.5f: closed form %.5f, from index symbols %.5f, L=%.4f" % (e, maclaurin(e), omega2(a), angular_momentum(a)))
    a = jacobi_axes(0.999999); print("Jacobi limit a2/a1->1: c/a=%.4f, Omega^2=%.4f (expect 0.5827, 0.3742)" % (a[2] / a[0], omega2(a)))
    for r in [0.8, 0.6, 0.43218, 0.3]:
        a = jacobi_axes(r); print("Jacobi a2/a1=%.4f: c/a=%.4f, Omega^2=%.4f, L=%.4f" % (r, a[2] / a[0], omega2(a), angular_momentum(a)))
