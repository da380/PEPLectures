"""Helpers for anisotropic elastic media: Voigt matrices, the fourth-order tensor, and
slowness sections from the Christoffel equation.  Used by the figure scripts."""
import numpy as np


def voigt_to_tensor(Cv):
    idx = {(0, 0): 0, (1, 1): 1, (2, 2): 2, (1, 2): 3, (2, 1): 3, (0, 2): 4, (2, 0): 4, (0, 1): 5, (1, 0): 5}
    A = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    A[i, j, k, l] = Cv[idx[(i, j)], idx[(k, l)]]
    return A

def ti_voigt(A, C, F, L, N):
    """Transversely isotropic medium with symmetry axis x3 (Love notation), in Voigt form."""
    Cv = np.zeros((6, 6))
    Cv[0, 0] = Cv[1, 1] = A; Cv[2, 2] = C
    Cv[0, 1] = Cv[1, 0] = A - 2 * N
    Cv[0, 2] = Cv[2, 0] = Cv[1, 2] = Cv[2, 1] = F
    Cv[3, 3] = Cv[4, 4] = L; Cv[5, 5] = N
    return Cv

def olivine_voigt():
    """San Carlos olivine at room conditions (GPa), Abramson et al. (1997)."""
    C11, C22, C33, C44, C55, C66 = 320.5, 196.5, 233.5, 64.0, 77.0, 78.7
    C12, C13, C23 = 68.1, 71.6, 76.8
    Cv = np.zeros((6, 6))
    Cv[0, 0], Cv[1, 1], Cv[2, 2] = C11, C22, C33
    Cv[3, 3], Cv[4, 4], Cv[5, 5] = C44, C55, C66
    Cv[0, 1] = Cv[1, 0] = C12; Cv[0, 2] = Cv[2, 0] = C13; Cv[1, 2] = Cv[2, 1] = C23
    return Cv

def slowness_section(A, rho, e1, e2, n=721):
    """Slowness (s/km) of the three sheets for directions in the plane spanned by e1, e2."""
    phi = np.linspace(0, 2 * np.pi, n)
    out = np.zeros((3, n))
    for m, ph in enumerate(phi):
        p = np.cos(ph) * e1 + np.sin(ph) * e2
        G = np.einsum("ijkl,j,l->ik", A, p, p) / rho       # GPa / (kg m^-3) = 1e9 m^2 s^-2 ... handled below
        c2 = np.sort(np.linalg.eigvalsh(G))                # ascending: slow to fast
        out[:, m] = 1.0 / np.sqrt(c2 * 1e9) * 1e3          # s/km
    return phi, out

