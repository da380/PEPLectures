"""Normal modes of PREM using Alex Myhill's mode_lab_2 (spectral-element solver for a layered,
self-gravitating, spherically symmetric earth model), plus the toroidal-mode sensitivity
kernels of Lecture 24.

The package returns eigenfunctions in its own nondimensional, mass-normalised form.  For the
kernels we renormalise dimensionally using the PREM density of pepseis.prem, so that with
    N = int rho W^2 r^2 dr = 1 ,
the toroidal eigenvalue problem of Lecture 24 (eq. 18 there) gives, by Rayleigh's principle,
    delta(omega^2) = int [ (r W' - W)^2 + (zeta^2 - 2) W^2 ] delta mu dr - omega^2 int W^2 r^2 delta rho dr,
i.e. delta omega = int (K_mu delta mu + K_rho delta rho) dr with
    K_mu = [ (r W' - W)^2 + (zeta^2 - 2) W^2 ] / (2 omega),   K_rho = - omega W^2 r^2 / 2 .
Checks: int K_mu mu dr = omega/2 and int K_rho rho dr = -omega/2 (scaling of mu or rho).
"""
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
import mode_lab_2 as ml
from pepseis import prem

DECK = Path.home() / "dev/mode_lab_2/data/models/prem_1s.deck"
CATALOGUES = Path.home() / "dev/mode_lab_2/data/results"
REFERENCES = Path.home() / "dev/mode_lab_2/data/references"


class PremModes:
    def __init__(self, order=5, element_size_km=100.0):
        self.model = ml.read_mineos_deck(DECK)
        self.mesh = ml.build_mesh(self.model, ml.MeshSpec(polynomial_order=order, maximum_element_size_m=element_size_km * 1e3))
        self._top = None

    # --- toroidal -------------------------------------------------------------------------
    @property
    def toroidal_operator(self):
        if self._top is None:
            self._top = ml.assemble_toroidal(self.model, self.mesh)
        return self._top

    def toroidal(self, l, n_max):
        """Frequencies (Hz) and radial eigenfunctions W(r) (r in km) of nT_l, n = 0..n_max."""
        sp = ml.solve_toroidal(self.toroidal_operator, l, n_max + 1)
        r = np.asarray(self.toroidal_operator.radius, float)          # km, one entry per dof
        order = np.argsort(r)
        return sp.frequencies_hz[: n_max + 1], r[order], sp.vectors[order, : n_max + 1]

    def toroidal_kernels(self, l, n):
        """K_mu(r) and K_rho(r) for nT_l, in units of s^-1 per Pa per m and s^-1 per kg m^-3 per m."""
        f, r, W = self.toroidal(l, n)
        omega = 2 * np.pi * f[n]; w = W[:, n]
        zeta2 = l * (l + 1)
        rm = r * 1e3
        rho = prem.density(np.clip(r, 3480.0 + 1e-6, 6368.0 - 1e-6))
        # normalise: int rho W^2 r^2 dr = 1  (SI)
        norm = np.trapezoid(rho * w ** 2 * rm ** 2, rm)
        w = w / np.sqrt(norm)
        # derivative by a spline through the nodes of each element is awkward at the interfaces;
        # a global cubic spline is adequate for plotting
        dw = CubicSpline(rm, w).derivative()(rm)
        K_mu = ((rm * dw - w) ** 2 + (zeta2 - 2) * w ** 2) / (2 * omega)
        K_rho = -omega * w ** 2 * rm ** 2 / 2
        return f[n], r, w, K_mu, K_rho

    # --- spheroidal -----------------------------------------------------------------------
    def spheroidal(self, l, f_max_hz=5e-3):
        """All finite spheroidal states of degree l below f_max (frequencies in Hz), with the
        assembled pencil for sampling; states are not classified."""
        if l == 0:
            sp = ml.solve_radial(self.model, self.mesh, maximum_frequency_hz=f_max_hz)
            pen = ml.assemble_radial(self.model, self.mesh)
        else:
            sp = ml.solve_spheroidal(self.model, self.mesh, l, maximum_frequency_hz=f_max_hz)
            pen = ml.assemble_spheroidal(self.model, self.mesh, l)
        good = np.isfinite(sp.frequencies_hz) & ~sp.numerical_zero
        return sp.frequencies_hz[good], sp.vectors[:, good], pen

    def spheroidal_mode(self, l, n, f_max_hz=5e-3, radii_km=None):
        """U(r), V(r), P(r) of nS_l, identified by the MINEOS reference frequency."""
        f_ref = reference_frequency("radial" if l == 0 else "spheroidal", l, n)
        f, vec, pen = self.spheroidal(l, f_max_hz=max(f_max_hz, 1.3 * f_ref))
        k = int(np.argmin(np.abs(f - f_ref)))
        if abs(f[k] - f_ref) > 0.02 * f_ref:
            raise RuntimeError("no state near the reference frequency of %dS%d" % (n, l))
        r = np.linspace(1.0, 6367.999, 600) if radii_km is None else np.asarray(radii_km, float)
        U, V, P = ml.sample_spheroidal_state(pen, vec[:, k], r * 1e3)
        return f[k], r, U, V, P


def reference_frequency(family, l, n):
    """Frequency (Hz) of a mode from the MINEOS reference tables shipped with mode_lab_2."""
    import csv
    name = {"spheroidal": "spheroidal_mineos_elastic.csv", "toroidal": "toroidal_mineos_elastic.csv", "radial": "radial_mineos_elastic.csv"}[family]
    with open(REFERENCES / name, newline="") as fh:
        for row in csv.DictReader(fh):
            if int(row["l"]) == l and int(row["n"]) == n and row["frequency_hz"]:
                return float(row["frequency_hz"])
    raise KeyError((family, l, n))


def catalogue(family):
    """(l, n, f_hz) arrays from the reviewed mode_lab_2 catalogues."""
    import csv
    if family == "toroidal":
        path = CATALOGUES / "03_toroidal_elastic_catalogue/catalogue.csv"; key = "sem_frequency_hz"; nkey = "n"
        rows = [r for r in csv.DictReader(open(path, newline=""))]
    else:
        path = CATALOGUES / "04_spheroidal_elastic_catalogue/s2_catalogue.csv"; key = "frequency_hz"; nkey = "comparison_n"
        rows = [r for r in csv.DictReader(open(path, newline="")) if r["accepted"] == "True"]
    l = np.array([int(r["l"]) for r in rows]); n = np.array([int(r[nkey]) if r[nkey] not in ("", "None") else -1 for r in rows])
    f = np.array([float(r[key]) for r in rows])
    return l, n, f


def _drop_zero_states(sp):
    """Copy of a SpheroidalSpectrum without its numerical-zero / non-finite states."""
    from mode_lab_2.spheroidal import SpheroidalSpectrum
    keep = np.isfinite(sp.frequencies_hz) & (sp.frequencies_hz > 0) & ~np.asarray(sp.numerical_zero, bool)
    return SpheroidalSpectrum(
        degree=sp.degree, raw_indices=sp.raw_indices[keep].copy(), eigenvalues=sp.eigenvalues[keep].copy(),
        frequencies_hz=sp.frequencies_hz[keep].copy(), vectors=sp.vectors[:, keep].copy(),
        action_residuals=sp.action_residuals[keep].copy(), scale_residuals=sp.scale_residuals[keep].copy(),
        zero_thresholds=sp.zero_thresholds[keep].copy(), numerical_zero=sp.numerical_zero[keep].copy())


def spheroidal_extension(l_min=201, l_max=500, f_max_hz=50e-3, order=5, element_size_km=50.0,
                         cache=None, verbose=False):
    """Extend the reviewed spheroidal catalogue (which stops at degree 200) to higher degrees, applying
    the same acceptance tests as mode_lab_2's catalogue (p- and h-refinement Rayleigh errors below 1%
    and less than half the energy in the ocean).  Returns (l, f_hz) arrays of accepted states, cached
    in an .npz file (default pepcode/data/modes/spheroidal_l%d-%d_%dmHz.npz)."""
    import os
    from mode_lab_2.mesh import same_boundaries_higher_order, one_extra_element_per_layer
    if cache is None:
        d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "modes")
        os.makedirs(d, exist_ok=True)
        cache = os.path.join(d, "spheroidal_l%d-%d_%dmHz.npz" % (l_min, l_max, round(1e3 * f_max_hz)))
    if os.path.exists(cache):
        z = np.load(cache)
        return z["l"], z["f"]
    pm = PremModes(order=order, element_size_km=element_size_km)
    mesh_p = same_boundaries_higher_order(pm.mesh, order + 1)
    mesh_h = one_extra_element_per_layer(pm.mesh)
    L, F = [], []
    for l in range(l_min, l_max + 1):
        sp = ml.solve_spheroidal(pm.model, pm.mesh, l, maximum_frequency_hz=f_max_hz)
        sp = _drop_zero_states(sp)
        base = ml.assemble_spheroidal(pm.model, pm.mesh, l)
        diag = ml.diagnose_spectrum(base, sp, p_target=ml.assemble_spheroidal(pm.model, mesh_p, l),
                                    h_target=ml.assemble_spheroidal(pm.model, mesh_h, l))
        for d in diag:
            if d.accepted:
                L.append(l); F.append(d.frequency_hz)
        if verbose:
            print(l, sum(d.accepted for d in diag), "of", len(diag), flush=True)
    L, F = np.array(L), np.array(F)
    np.savez(cache, l=L, f=F)
    return L, F


# --- spheroidal sensitivity kernels by re-assembly ------------------------------------------------
# mode_lab_2 assembles its matrices from a sampled Material; the matrices are linear in the material
# parameters, so Rayleigh's principle applied to the unperturbed eigenvector with a pencil assembled
# for a perturbed material gives the exact first-order frequency change.  Perturbing one element at a
# time gives the kernels averaged over each element.

def _gravity_nd(mesh, rho):
    """Nondimensional gravity consistent with a nodal density field (matches material.gravity)."""
    from mode_lab_2.spheroidal import G_ND
    g = np.zeros_like(rho); cum = 0.0
    for e in range(len(mesh.element_layers)):
        r = np.asarray(mesh.radius[e], float); w = float(mesh.jacobian[e]) * np.asarray(mesh.gll.weights, float)
        f = rho[e] * r ** 2
        # cumulative integral inside the element via the GLL nodes (integral to each node by trapezoid on
        # top of the exact element total keeps the element totals consistent)
        inner = cum + np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(r))))
        total = cum + float(np.sum(w * f))
        inner = cum + (inner - cum) * (total - cum) / max(inner[-1] - cum, 1e-300)
        with np.errstate(divide="ignore", invalid="ignore"):
            g[e] = np.where(r > 0, 4.0 * np.pi * G_ND * inner / r ** 2, 0.0)
        cum = total
    return g


def _perturbed_pencil(pen, material, mesh, l):
    from mode_lab_2.spheroidal import SpheroidalPencil, _assemble_nonradial, _assemble_radial, _dof_partitions
    if l == 0:
        K, M = _assemble_radial(mesh, material)
    else:
        K, M = _assemble_nonradial(mesh, material, l)
    dyn, sta, disp, pot = _dof_partitions(mesh, M, radial=(l == 0))
    return SpheroidalPencil(l, K, M, mesh, material, dyn, sta, disp, pot)


def spheroidal_kernels(pm, l, n, f_max_hz=None, elements=None):
    """Element-averaged kernels K_kappa, K_mu, K_rho (SI: s^-1 per Pa per m, s^-1 per kg m^-3 per m)
    for nS_l in PREM such that delta omega = int (K_kappa d kappa + K_mu d mu + K_rho d rho) dr.
    Returns radius (km) at element centres and the three kernels (NaN for mu in fluid elements)."""
    import dataclasses
    from mode_lab_2.model import RHO_NORM, STRESS_NORM_PA
    from mode_lab_2.catalogue import rayleigh_frequency_hz
    f_ref = reference_frequency("radial" if l == 0 else "spheroidal", l, n)
    if f_max_hz is None:
        f_max_hz = 1.3 * f_ref
    f, vec, pen = pm.spheroidal(l, f_max_hz)
    k = int(np.argmin(np.abs(f - f_ref)))
    if abs(f[k] / f_ref - 1) > 2e-3:
        raise RuntimeError("no state near the reference frequency of %dS%d" % (n, l))
    s = vec[:, k]; mesh = pm.mesh; mat = pen.material
    f0 = rayleigh_frequency_hz(pen, s)
    ne = len(mesh.element_layers)
    elements = range(ne) if elements is None else elements
    rc = np.array([0.5 * (mesh.element_bounds_m[e, 0] + mesh.element_bounds_m[e, 1]) for e in elements]) / 1e3
    h = np.array([mesh.element_bounds_m[e, 1] - mesh.element_bounds_m[e, 0] for e in elements])
    out = {"kappa": [], "mu": [], "rho": []}
    eps = 1e-3
    for e in elements:
        solid = bool(mat.solid_element[e])
        for name in out:
            if name == "mu" and not solid:
                out[name].append(np.nan); continue
            rho, A, C, F, L, N = (np.array(getattr(mat, x), float) for x in ("rho", "A", "C", "F", "L", "N"))
            if name == "rho":
                d = eps * rho[e].mean(); rho[e] += d; g = _gravity_nd(mesh, rho)
                m2 = dataclasses.replace(mat, rho=rho, gravity=g); dsi = d * RHO_NORM
            elif name == "kappa":
                d = eps * A[e].mean(); A[e] += d; C[e] += d; F[e] += d
                m2 = dataclasses.replace(mat, A=A, C=C, F=F); dsi = d * STRESS_NORM_PA
            else:
                d = eps * L[e].mean(); A[e] += 4 * d / 3; C[e] += 4 * d / 3; F[e] -= 2 * d / 3; L[e] += d; N[e] += d
                m2 = dataclasses.replace(mat, A=A, C=C, F=F, L=L, N=N); dsi = d * STRESS_NORM_PA
            f1 = rayleigh_frequency_hz(_perturbed_pencil(pen, m2, mesh, l), s)
            out[name].append(2 * np.pi * (f1 - f0) / dsi)
    K = {name: np.array(v) / h for name, v in out.items()}
    return f0, rc, K["kappa"], K["mu"], K["rho"], (s, pen, k)


def modal_q_and_dispersion(pm, l, n):
    """Quality factor of nS_l from the kernels and PREM's Q_kappa, Q_mu, and the frequency corrected
    for physical dispersion from the 1 s reference period of the elastic model (constant-Q absorption
    band: f = f_e [1 + ln(f_e T_ref) / (pi Q)]).  Returns (f_elastic, Q, f_anelastic) in Hz."""
    f0, rc, Kk, Km, Kr, _ = spheroidal_kernels(pm, l, n)
    rows = pm.model.rows
    rr = np.array([row.radius_m for row in rows]) / 1e3
    qk = np.array([row.q_kappa for row in rows]); qm = np.array([row.q_mu for row in rows])
    inv_qk = np.interp(rc, rr, np.where(qk > 0, 1 / np.where(qk > 0, qk, 1), 0.0))
    inv_qm = np.interp(rc, rr, np.where(qm > 0, 1 / np.where(qm > 0, qm, 1), 0.0))
    kap, mu = prem.moduli(np.clip(rc, 1e-3, 6368 - 1e-6))
    om = 2 * np.pi * f0; dr = np.gradient(rc) * 1e3
    inv_q = 2 * np.sum((Kk * kap * inv_qk + np.nan_to_num(Km) * mu * inv_qm) * dr) / om
    Q = 1 / inv_q
    f1 = f0 * (1 + np.log(f0 * pm.model.reference_period_s) / (np.pi * Q))
    return f0, Q, f1
