"""The Preliminary Reference Earth Model (Dziewonski & Anderson 1981), isotropic version.

Density and wave speeds are polynomials in x = r / 6371 km on thirteen layers.  The
anisotropic layer between 24.4 km and 220 km depth is represented by the isotropic
velocities given in the paper.  Units: km, km/s, kg/m^3 (density), m/s^2 (gravity), Pa.
"""
import numpy as np
from scipy.integrate import quad

R_EARTH = 6371.0
R_CMB = 3480.0
R_ICB = 1221.5
G = 6.674e-11

# layer boundaries (radius, km) and polynomial coefficients (lowest order first) for
# density [g/cm^3], Vp [km/s], Vs [km/s]
_LAYERS = [
    (0.0,     1221.5, [13.0885, 0, -8.8381],            [11.2622, 0, -6.3640],                [3.6678, 0, -4.4475]),
    (1221.5,  3480.0, [12.5815, -1.2638, -3.6426, -5.5281], [11.0487, -4.0362, 4.8023, -13.5732], [0.0]),
    (3480.0,  3630.0, [7.9565, -6.4761, 5.5283, -3.0807], [15.3891, -5.3181, 5.5242, -2.5514],  [6.9254, 1.4672, -2.0834, 0.9783]),
    (3630.0,  5600.0, [7.9565, -6.4761, 5.5283, -3.0807], [24.9520, -40.4673, 51.4832, -26.6419], [11.1671, -13.7818, 17.4575, -9.2777]),
    (5600.0,  5701.0, [7.9565, -6.4761, 5.5283, -3.0807], [29.2766, -23.6027, 5.5242, -2.5514],  [22.3459, -17.2473, -2.0834, 0.9783]),
    (5701.0,  5771.0, [5.3197, -1.4836],                 [19.0957, -9.8672],                   [9.9839, -4.9324]),
    (5771.0,  5971.0, [11.2494, -8.0298],                [39.7027, -32.6166],                  [22.3512, -18.5856]),
    (5971.0,  6151.0, [7.1089, -3.8045],                 [20.3926, -12.2569],                  [8.9496, -4.4597]),
    (6151.0,  6346.6, [2.6910, 0.6924],                  [4.1875, 3.9382],                     [2.1519, 2.3481]),
    (6346.6,  6356.0, [2.9],                             [6.8],                                [3.9]),
    (6356.0,  6368.0, [2.6],                             [5.8],                                [3.2]),
    (6368.0,  6371.0, [1.02],                            [1.45],                               [0.0]),
]
BOUNDARIES = np.array([l[0] for l in _LAYERS] + [R_EARTH])


def _layer_index(r):
    r = np.asarray(r, float)
    idx = np.searchsorted(BOUNDARIES, r, side="right") - 1
    return np.clip(idx, 0, len(_LAYERS) - 1)


def _poly(coeffs, x):
    return sum(c * x ** n for n, c in enumerate(coeffs))


def _evaluate(r, which, side="upper"):
    """Evaluate rho (g/cm^3), Vp or Vs at radius r; at a boundary 'upper' takes the layer
    above (larger r), 'lower' the layer below."""
    r = np.asarray(r, float)
    scalar = r.ndim == 0
    r = np.atleast_1d(r)
    idx = _layer_index(r)
    if side == "lower":
        on_boundary = np.isin(r, BOUNDARIES[1:-1])
        idx = np.where(on_boundary, idx - 1, idx)
    out = np.zeros_like(r)
    x = r / R_EARTH
    for k, lay in enumerate(_LAYERS):
        m = idx == k
        if m.any():
            out[m] = _poly(lay[which], x[m])
    return out[0] if scalar else out


def density(r, side="upper"):
    """Density in kg/m^3."""
    return 1000.0 * _evaluate(r, 2, side)

def vp(r, side="upper"):
    return _evaluate(r, 3, side)

def vs(r, side="upper"):
    return _evaluate(r, 4, side)

def moduli(r, side="upper"):
    """Bulk and shear moduli in Pa."""
    rho, a, b = density(r, side), 1e3 * vp(r, side), 1e3 * vs(r, side)
    mu = rho * b ** 2
    kappa = rho * a ** 2 - 4.0 * mu / 3.0
    return kappa, mu


def mass(r):
    """Mass (kg) within radius r (km)."""
    r = float(r)
    total = 0.0
    for lay in _LAYERS:
        lo, hi = lay[0], min(lay[1], r)
        if hi <= lo:
            break
        total += quad(lambda s: density(s) * 4 * np.pi * (s * 1e3) ** 2 * 1e3, lo, hi)[0]
    return total


def gravity(r):
    """Gravitational acceleration (m/s^2) at radius r (km)."""
    r = float(r)
    return G * mass(r) / (r * 1e3) ** 2 if r > 0 else 0.0


def pressure(r):
    """Hydrostatic pressure (Pa) at radius r (km)."""
    return quad(lambda s: density(s) * gravity(s) * 1e3, r, R_EARTH, limit=200,
                points=[b for b in BOUNDARIES if r < b < R_EARTH])[0]


if __name__ == "__main__":
    print("mass  %.4e kg  (PREM: 5.974e24)" % mass(R_EARTH))
    print("g(surface) %.3f  g(CMB) %.3f m/s^2 (PREM: 9.82, 10.68)" % (gravity(R_EARTH), gravity(R_CMB)))
    print("p(CMB) %.1f GPa  p(centre) %.1f GPa (PREM: 135.8, 364)" % (pressure(R_CMB) / 1e9, pressure(1.0) / 1e9))
    for b in BOUNDARIES[1:-1]:
        print("r=%7.1f  rho %6.0f|%6.0f  vp %6.3f|%6.3f  vs %6.3f|%6.3f" % (
            b, density(b, "lower"), density(b, "upper"), vp(b, "lower"), vp(b, "upper"), vs(b, "lower"), vs(b, "upper")))
