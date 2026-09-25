"""Two-dimensional Hamiltonian ray tracing for isotropic P-waves.

The ray Hamiltonian is H(x, p) = (1/2) c(x)^2 |p|^2, and rays are the projections onto x of
solutions of Hamilton's equations dx/ds = dH/dp, dp/ds = -dH/dx with H = 1/2, so that the
generating parameter s is the travel time.  The wave speed c(x, z) is given on a grid and
interpolated with a bicubic spline, which also supplies the gradient needed for dp/ds.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline


class Medium:
    def __init__(self, x, z, c):
        """x, z: 1-D grids; c: wave speed array of shape (len(x), len(z))."""
        self.x, self.z = x, z
        self.c = c
        self._spl = RectBivariateSpline(x, z, c, kx=3, ky=3)

    def speed(self, x, z):
        return self._spl(x, z, grid=False)

    def grad(self, x, z):
        return (self._spl(x, z, dx=1, grid=False), self._spl(x, z, dy=1, grid=False))


def trace_ray(medium, x0, p0, t_max, bounds, max_step=None):
    """Integrate Hamilton's equations from (x0, p0) for travel time t_max, stopping when
    the ray leaves the rectangle bounds = (xmin, xmax, zmin, zmax).  Returns (t, x, z)."""
    xmin, xmax, zmin, zmax = bounds

    def rhs(t, y):
        x, z, px, pz = y
        c = medium.speed(x, z)
        cx, cz = medium.grad(x, z)
        # dx/dt = c^2 p ;  dp/dt = -c |p|^2 grad c = -(1/c) grad c   (using |p| = 1/c)
        return [c * c * px, c * c * pz, -cx / c, -cz / c]

    def leave(t, y):
        x, z = y[0], y[1]
        return min(x - xmin, xmax - x, z - zmin, zmax - z)
    leave.terminal = True
    leave.direction = -1

    sol = solve_ivp(rhs, (0.0, t_max), [x0[0], x0[1], p0[0], p0[1]], events=leave,
                    max_step=max_step or t_max / 400, rtol=1e-8, atol=1e-10, dense_output=False)
    return sol.t, sol.y[0], sol.y[1]


def plane_wave_fan(medium, z_start, x_start, n_rays, t_max, bounds):
    """Rays leaving the line x = x_start (a plane wavefront) travelling in +x."""
    rays = []
    for z0 in np.linspace(z_start[0], z_start[1], n_rays):
        c0 = float(medium.speed(x_start, z0))
        rays.append(trace_ray(medium, (x_start, z0), (1.0 / c0, 0.0), t_max, bounds))
    return rays


def gaussian_lens(x, z, c0=1.0, dc=-0.25, centre=(2.0, 1.5), width=0.4):
    """Smooth low-velocity anomaly: c = c0 [1 + dc exp(-r^2 / 2 width^2)]."""
    X, Z = np.meshgrid(x, z, indexing="ij")
    r2 = (X - centre[0]) ** 2 + (Z - centre[1]) ** 2
    return c0 * (1.0 + dc * np.exp(-r2 / (2.0 * width ** 2)))


def random_medium(x, z, c0=1.0, rms=0.05, corr_len=0.4, seed=3, x_onset=None):
    """Gaussian random field with a Gaussian correlation function, made by filtering white
    noise in the Fourier domain; rms is the relative velocity fluctuation.  If x_onset is
    given, the fluctuations are switched on smoothly over x_onset[0] < x < x_onset[1] so that
    the medium is homogeneous to the left."""
    rng = np.random.default_rng(seed)
    nx, nz = len(x), len(z)
    dx, dz = x[1] - x[0], z[1] - z[0]
    kx = 2 * np.pi * np.fft.fftfreq(nx, dx)
    kz = 2 * np.pi * np.fft.fftfreq(nz, dz)
    KX, KZ = np.meshgrid(kx, kz, indexing="ij")
    filt = np.exp(-(KX ** 2 + KZ ** 2) * corr_len ** 2 / 4.0)
    noise = rng.standard_normal((nx, nz))
    field = np.real(np.fft.ifft2(np.fft.fft2(noise) * filt))
    field *= rms / field.std()
    if x_onset is not None:
        ramp = np.clip((x - x_onset[0]) / (x_onset[1] - x_onset[0]), 0.0, 1.0)
        ramp = 0.5 - 0.5 * np.cos(np.pi * ramp)
        field *= ramp[:, None]
    return c0 * (1.0 + field)
