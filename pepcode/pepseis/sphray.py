"""Ray tracing in a spherically symmetric earth model (PREM), using the radial integrals
of Lecture 18:  dDelta/dr = p v / (r sqrt(r^2 - p^2 v^2)),  dT/dr = r / (v sqrt(r^2 - p^2 v^2)),
with p = r sin(theta) / v the ray parameter (s/rad), conserved along the ray and across
interfaces (Snell's law).  A seismic phase is a sequence of legs through the mantle, outer
core and inner core; each leg either turns within its region ('turn'), or crosses it one way
('down' or 'up').  The integrable singularity at a turning point is removed by the
substitution r = r_t + s^2.
"""
import warnings
import numpy as np
from scipy.integrate import quad, IntegrationWarning
from scipy.optimize import brentq
from pepseis import prem

warnings.simplefilter("ignore", IntegrationWarning)   # the turning-point substitution converges; quad merely reports subdivisions

R_TOP = 6368.0          # base of the ocean layer: rays start here ("surface" focus)
REGIONS = {"mantle": (prem.R_CMB, R_TOP), "outer": (prem.R_ICB, prem.R_CMB), "inner": (0.0, prem.R_ICB)}


def _speed(wave):
    return prem.vp if wave == "P" else prem.vs


def _sub_layers(r_lo, r_hi):
    """Split [r_lo, r_hi] at the PREM boundaries."""
    pts = [r_lo] + [b for b in prem.BOUNDARIES if r_lo < b < r_hi] + [r_hi]
    return list(zip(pts[:-1], pts[1:]))


def turning_radius(p, wave, r_lo, r_hi):
    """Radius at which r / v(r) = p, scanning downwards from r_hi.  Returns (r_t, status) with
    status 'turn' (turning point found), 'through' (ray reaches r_lo) or 'reflect' (the ray
    cannot enter a layer: total reflection at its top)."""
    v = _speed(wave)
    for a, b in reversed(_sub_layers(r_lo, r_hi)):
        f_top = b / v(b, "lower") - p          # just inside the layer at its top
        f_bot = a / v(a, "upper") - p          # just inside at its bottom
        if f_top <= 0:
            return b, "reflect"
        if f_bot > 0:
            continue
        rt = brentq(lambda r: r / v(r) - p, a + 1e-9, b - 1e-9)
        return rt, "turn"
    return r_lo, "through"


def _integrands(p, wave):
    v = _speed(wave)
    def dDelta(r):
        vv = v(r); q = r * r - p * p * vv * vv
        return p * vv / (r * np.sqrt(q)) if q > 0 else 0.0
    def dT(r):
        vv = v(r); q = r * r - p * p * vv * vv
        return r / (vv * np.sqrt(q)) if q > 0 else 0.0
    return dDelta, dT


def _one_way(p, wave, r_lo, r_hi, r_turn=None):
    """Delta and T for one traverse between r_lo and r_hi (r_lo may be a turning point)."""
    dDelta, dT = _integrands(p, wave)
    D = T = 0.0
    for a, b in _sub_layers(r_lo, r_hi):
        if r_turn is not None and abs(a - r_turn) < 1e-9:
            # substitution r = r_t + s^2 removes the 1/sqrt singularity at the turning point
            smax = np.sqrt(b - a)
            D += quad(lambda s: dDelta(a + s * s) * 2 * s, 0, smax, limit=100)[0]
            T += quad(lambda s: dT(a + s * s) * 2 * s, 0, smax, limit=100)[0]
        else:
            D += quad(dDelta, a, b, limit=100)[0]
            T += quad(dT, a, b, limit=100)[0]
    return D, T


def leg(p, wave, region, kind):
    """(Delta, T) for one leg, or None if the leg is not possible for this ray parameter."""
    r_lo, r_hi = REGIONS[region]
    rt, status = turning_radius(p, wave, r_lo, r_hi)
    if kind == "turn":
        if status != "turn":
            return None
        D, T = _one_way(p, wave, rt, r_hi, r_turn=rt)
        return 2 * D, 2 * T
    if status != "through":
        return None
    return _one_way(p, wave, r_lo, r_hi)


PHASES = {
    "P":     [("P", "mantle", "turn")],
    "PP":    [("P", "mantle", "turn")] * 2,
    "PcP":   [("P", "mantle", "down"), ("P", "mantle", "up")],
    "PKP":   [("P", "mantle", "down"), ("P", "outer", "turn"), ("P", "mantle", "up")],
    "PKiKP": [("P", "mantle", "down"), ("P", "outer", "down"), ("P", "outer", "up"), ("P", "mantle", "up")],
    "PKIKP": [("P", "mantle", "down"), ("P", "outer", "down"), ("P", "inner", "turn"), ("P", "outer", "up"), ("P", "mantle", "up")],
    "S":     [("S", "mantle", "turn")],
    "SS":    [("S", "mantle", "turn")] * 2,
    "ScS":   [("S", "mantle", "down"), ("S", "mantle", "up")],
    "SKS":   [("S", "mantle", "down"), ("P", "outer", "turn"), ("S", "mantle", "up")],
    "SKKS":  [("S", "mantle", "down"), ("P", "outer", "turn"), ("P", "outer", "turn"), ("S", "mantle", "up")],
    "SKIKS": [("S", "mantle", "down"), ("P", "outer", "down"), ("P", "inner", "turn"), ("P", "outer", "up"), ("S", "mantle", "up")],
}


def trace(p, phase):
    """Total (Delta [rad], T [s]) for a phase at ray parameter p, or None."""
    D = T = 0.0
    for wave, region, kind in PHASES[phase]:
        res = leg(p, wave, region, kind)
        if res is None:
            return None
        D += res[0]; T += res[1]
    return D, T


def travel_time_curve(phase, n=400, p_max=None):
    """Arrays of (Delta [deg], T [s], p) sampled over the ray parameters for which the phase exists."""
    p_max = p_max or R_TOP / prem.vp(R_TOP, "lower") * 1.001
    out = []
    for p in np.linspace(1e-3, p_max, n):
        res = trace(p, phase)
        if res is not None:
            out.append((np.degrees(res[0]), res[1], p))
    return np.array(out).T if out else np.zeros((3, 0))


def ray_path(p, phase, n_per_layer=25):
    """Points (Delta [rad], r [km]) along the ray path of a phase, for plotting."""
    if trace(p, phase) is None:
        return None
    pts = [(0.0, R_TOP)]
    D = 0.0
    for wave, region, kind in PHASES[phase]:
        r_lo, r_hi = REGIONS[region]
        rt, status = turning_radius(p, wave, r_lo, r_hi)
        r_bot = rt if kind == "turn" else r_lo
        segs = _sub_layers(r_bot, r_hi)
        def down_points():
            nonlocal D
            for a, b in reversed(segs):
                rs = np.linspace(b, a, n_per_layer)
                if kind == "turn" and abs(a - rt) < 1e-9:
                    rs = rt + (b - rt) * np.linspace(1, 0, n_per_layer) ** 2
                for r0, r1 in zip(rs[:-1], rs[1:]):
                    D += _one_way(p, wave, r1, r0, r_turn=rt if (kind == "turn" and abs(r1 - rt) < 1e-9) else None)[0]
                    pts.append((D, r1))
        def up_points():
            nonlocal D
            for a, b in segs:
                rs = np.linspace(a, b, n_per_layer)
                if kind == "turn" and abs(a - rt) < 1e-9:
                    rs = rt + (b - rt) * np.linspace(0, 1, n_per_layer) ** 2
                for r0, r1 in zip(rs[:-1], rs[1:]):
                    D += _one_way(p, wave, r0, r1, r_turn=rt if (kind == "turn" and abs(r0 - rt) < 1e-9) else None)[0]
                    pts.append((D, r1))
        if kind == "turn":
            down_points(); up_points()
        elif kind == "down":
            down_points()
        else:
            up_points()
    return np.array(pts)


if __name__ == "__main__":
    # checks against standard tables (ak135 / PREM, surface focus)
    for phase, target in [("PcP", 0.0), ("ScS", 0.0), ("PKIKP", 180.0), ("P", 60.0), ("S", 60.0), ("PKiKP", 0.0)]:
        d, t, pp = travel_time_curve(phase, n=2000)
        k = np.argmin(abs(d - target))
        print("%-6s Delta=%6.1f  T=%7.1f s = %2d min %4.1f s" % (phase, d[k], t[k], t[k] // 60, t[k] % 60))
    d, t, pp = travel_time_curve("P", n=2000); print("P reaches Delta_max = %.1f deg" % d.max())
    d, t, pp = travel_time_curve("PKP", n=2000); print("PKP Delta range %.1f - %.1f deg, min %.1f" % (d.min(), d.max(), d.min()))
    d, t, pp = travel_time_curve("PKIKP", n=2000); print("PKIKP Delta range %.1f - %.1f deg" % (d.min(), d.max()))
