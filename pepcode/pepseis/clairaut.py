"""Clairaut's equation for the hydrostatic ellipticity of a slowly rotating, spherically
stratified planet (Lecture 22, eqs. 42-43):
    eps'' + (8 pi G rho0 / g0) (eps' + eps / r) - 6 eps / r^2 = 0,
    eps'(0) = 0,   eps'(b) = [ 5 Omega^2 b^3 / (2 G M) - 2 eps(b) ] / b.
The equation is linear and homogeneous, so the regular solution is integrated outwards from
the centre (eps = const there) layer by layer, with eps and eps' continuous across density
discontinuities, and then scaled to satisfy the surface condition.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from pepseis import prem


def ellipticity(omega=7.292115e-5, n_g=600):
    b = prem.R_EARTH
    rg = np.linspace(1.0, b, n_g)
    g_tab = interp1d(rg, [prem.gravity(r) for r in rg], kind="cubic", fill_value="extrapolate")
    M = prem.mass(b)

    # work in km throughout: eps' has units 1/km, and 8 pi G rho / g has units 1/m -> convert to 1/km
    def rhs_km(r, y, lo, hi):
        eps, deps = y
        rho = prem.density(min(max(r, lo + 1e-6), hi - 1e-6))
        g = g_tab(r)
        coeff = 8 * np.pi * prem.G * rho / g * 1e3          # per km
        return [deps, -coeff * (deps + eps / r) + 6 * eps / r ** 2]

    r_all, e_all = [], []
    y = [1.0, 0.0]; r0 = 1.0
    for lo, hi in zip(prem.BOUNDARIES[:-1], prem.BOUNDARIES[1:]):
        start = max(lo, r0)
        sol = solve_ivp(rhs_km, (start, hi), y, args=(lo, hi), rtol=1e-9, atol=1e-12, dense_output=True, max_step=5.0)
        rr = np.linspace(start, hi, 60); r_all.append(rr); e_all.append(sol.sol(rr)[0])
        y = sol.y[:, -1]
    r = np.concatenate(r_all); eps_h = np.concatenate(e_all)
    eps_b, deps_b = y                                        # unscaled surface values (per km)
    q = omega ** 2 * (b * 1e3) ** 3 / (prem.G * M)
    # surface condition: A (deps_b + 2 eps_b / b) = 5 q / (2 b)
    A = (5 * q / (2 * b)) / (deps_b + 2 * eps_b / b)
    return r, A * eps_h


if __name__ == "__main__":
    r, eps = ellipticity()
    print("1/eps at surface: %.1f   at centre: %.1f" % (1 / eps[-1], 1 / eps[0]))
