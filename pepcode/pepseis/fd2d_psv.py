"""A small 2D elastic (P-SV) finite-difference solver on a staggered grid.

Velocity-stress formulation (Virieux 1986), fourth-order differences in space (Levander
1988) and second-order leap-frog in time.  The grid is periodic in both directions, with
absorbing sponge layers at the left and right edges.  Written for the lecture figures;
plane-wave initial conditions are provided, and the displacement is integrated alongside
the velocity so that snapshots of the displacement field can be drawn.

Fields (all arrays of shape (nx, nz), index [i, j] for x and z):
    vx  at (i+1/2, j)      vz  at (i, j+1/2)
    txx, tzz at (i, j)     txz at (i+1/2, j+1/2)
"""
import numpy as np

C1, C2 = 9.0 / 8.0, 1.0 / 24.0


def dxf(f, h):   # forward difference in x, result at i+1/2
    return (C1 * (np.roll(f, -1, 0) - f) - C2 * (np.roll(f, -2, 0) - np.roll(f, 1, 0))) / h

def dxb(f, h):   # backward difference in x, result at i-1/2
    return (C1 * (f - np.roll(f, 1, 0)) - C2 * (np.roll(f, -1, 0) - np.roll(f, 2, 0))) / h

def dzf(f, h):
    return (C1 * (np.roll(f, -1, 1) - f) - C2 * (np.roll(f, -2, 1) - np.roll(f, 1, 1))) / h

def dzb(f, h):
    return (C1 * (f - np.roll(f, 1, 1)) - C2 * (np.roll(f, -1, 1) - np.roll(f, 2, 1))) / h


class PSVSolver:
    def __init__(self, x, z, alpha, beta, rho, sponge_width=0.15, sponge_strength=0.015, dt=None, sponge_z=False):
        """sponge_z: also absorb at the z edges (by default the grid is periodic in z, which
        suits plane waves travelling in x)."""
        self.x, self.z = x, z
        self.h = x[1] - x[0]
        self.nx, self.nz = len(x), len(z)
        mu = rho * beta ** 2
        lam = rho * alpha ** 2 - 2.0 * mu
        self.lam, self.mu, self.rho = lam, mu, rho
        # material parameters at the staggered positions (arithmetic averages)
        self.mu_xz = 0.25 * (mu + np.roll(mu, -1, 0) + np.roll(mu, -1, 1) + np.roll(np.roll(mu, -1, 0), -1, 1))
        self.rho_x = 0.5 * (rho + np.roll(rho, -1, 0))
        self.rho_z = 0.5 * (rho + np.roll(rho, -1, 1))
        vmax = np.max(alpha)
        self.dt = dt if dt is not None else 0.5 * self.h / vmax
        # Cerjan-type sponge at the x edges
        n = int(sponge_width / self.h)
        d = np.ones(self.nx)
        k = np.arange(n)
        d[:n] = np.exp(-(sponge_strength * (n - k)) ** 2)
        d[-n:] = np.exp(-(sponge_strength * (k + 1)) ** 2)
        self.damp = d[:, None]
        if sponge_z:
            nzs = int(sponge_width / self.h); dz_ = np.ones(self.nz); kz = np.arange(nzs)
            dz_[:nzs] = np.exp(-(sponge_strength * (nzs - kz)) ** 2); dz_[-nzs:] = np.exp(-(sponge_strength * (kz + 1)) ** 2)
            self.damp = d[:, None] * dz_[None, :]
        self.t = 0.0
        self.vx = np.zeros((self.nx, self.nz)); self.vz = np.zeros_like(self.vx)
        self.txx = np.zeros_like(self.vx); self.tzz = np.zeros_like(self.vx); self.txz = np.zeros_like(self.vx)
        self.ux = np.zeros_like(self.vx); self.uz = np.zeros_like(self.vx)

    def plane_p_wave(self, x0, width, alpha0):
        """Right-going plane P-wave with Gaussian displacement pulse u_x = F(x - alpha t)."""
        X = self.x[:, None] * np.ones((1, self.nz))
        F = np.exp(-(X - x0) ** 2 / (2 * width ** 2))
        Fp = -(X - x0) / width ** 2 * F
        # the leap-frog scheme holds the velocity half a time step behind the stresses, so
        # the velocity is initialised at t = -dt/2 (argument x + alpha dt/2) to avoid a
        # spurious static residual at the starting position
        # (vx also sits half a cell to the right of the stresses on the staggered grid)
        Xh = X + 0.5 * self.h + 0.5 * alpha0 * self.dt
        Fph = -(Xh - x0) / width ** 2 * np.exp(-(Xh - x0) ** 2 / (2 * width ** 2))
        # the displacement is integrated from vx, so it is sampled at the same staggered position
        self.ux = np.exp(-(X + 0.5 * self.h - x0) ** 2 / (2 * width ** 2))
        self.vx = -alpha0 * Fph
        self.txx = (self.lam + 2 * self.mu) * Fp
        self.tzz = self.lam * Fp

    def step(self, force=None, explosion=None):
        """Advance one time step.  force = (ix, iz, fx, fz) adds a point body force (per unit
        mass, i.e. an acceleration) at grid point (ix, iz); explosion = (ix, iz, q) adds q to
        both normal stresses there (an isotropic moment-rate source)."""
        h, dt = self.h, self.dt
        self.vx += dt / self.rho_x * (dxf(self.txx, h) + dzb(self.txz, h))
        self.vz += dt / self.rho_z * (dxb(self.txz, h) + dzf(self.tzz, h))
        if force is not None:
            ix, iz, fx, fz = force
            self.vx[ix, iz] += dt * fx; self.vz[ix, iz] += dt * fz
        self.vx *= self.damp; self.vz *= self.damp
        exx = dxb(self.vx, h); ezz = dzb(self.vz, h)
        self.txx += dt * ((self.lam + 2 * self.mu) * exx + self.lam * ezz)
        self.tzz += dt * (self.lam * exx + (self.lam + 2 * self.mu) * ezz)
        self.txz += dt * self.mu_xz * (dzf(self.vx, h) + dxf(self.vz, h))
        if explosion is not None:
            ix, iz, q = explosion
            self.txx[ix, iz] += dt * q; self.tzz[ix, iz] += dt * q
        self.txx *= self.damp; self.tzz *= self.damp; self.txz *= self.damp
        self.ux += dt * self.vx; self.uz += dt * self.vz
        self.t += dt

    def divergence(self):
        """div u, obtained from the normal stresses: txx + tzz = 2 (lambda + mu) div u."""
        return (self.txx + self.tzz) / (2.0 * (self.lam + self.mu))

    def run(self, t_end):
        while self.t < t_end - 1e-12:
            self.step()
        return self.ux, self.uz


def gaussian_anomaly(x, z, centre, width, frac):
    """Fractional velocity perturbation frac * exp(-r^2 / 2 width^2)."""
    X, Z = np.meshgrid(x, z, indexing="ij")
    return frac * np.exp(-((X - centre[0]) ** 2 + (Z - centre[1]) ** 2) / (2 * width ** 2))


def _draw(ax, solver, anomaly, quiver_step, vmax):
    x, z = solver.x, solver.z
    mag = np.hypot(solver.ux, solver.uz)
    im = ax.pcolormesh(x, z, mag.T, cmap="Blues", vmin=0, vmax=vmax or mag.max(), shading="auto", rasterized=True)
    q = quiver_step
    U = solver.ux[q // 2::q, q // 2::q].T.copy(); V = solver.uz[q // 2::q, q // 2::q].T.copy()
    small = np.hypot(U, V) < 0.02 * (vmax or mag.max())
    U[small] = np.nan; V[small] = np.nan
    ax.quiver(x[q // 2::q], z[q // 2::q], U, V, color="k", scale=16, width=0.002, headwidth=4, alpha=0.9)
    ax.contour(x, z, anomaly.T, levels=[0.5 * anomaly.min()], colors="k", linestyles="--", linewidths=0.8)
    ax.set_aspect("equal"); ax.set_xlim(x[0], x[-1]); ax.set_ylim(z[0], z[-1])
    ax.set_xticks([]); ax.set_yticks([])
    return im


def animate(solver, anomaly, fname, t_end, n_frames=60, quiver_step=18, vmax=1.0):
    """Save an animated GIF of the displacement field from the current time to t_end."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})
    times = np.linspace(solver.t, t_end, n_frames)
    fig, ax = plt.subplots(figsize=(7.4, 3.9))

    def frame(k):
        solver.run(times[k])
        ax.clear()
        _draw(ax, solver, anomaly, quiver_step, vmax)
        ax.set_title(f"t = {solver.t:.2f}", fontsize=10)
    anim = FuncAnimation(fig, frame, frames=n_frames)
    anim.save(fname + ".gif", writer=PillowWriter(fps=12), dpi=80)
    plt.close(fig)


def snapshot_figure(solver, anomaly, fname, quiver_step=18, vmax=None):
    """Displacement magnitude with arrows for the particle motion; the dashed line is the
    half-maximum contour of the velocity anomaly."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})
    x, z = solver.x, solver.z
    mag = np.hypot(solver.ux, solver.uz)
    fig, ax = plt.subplots(figsize=(7.4, 3.9))
    im = _draw(ax, solver, anomaly, quiver_step, vmax)
    cb = fig.colorbar(im, ax=ax, fraction=0.024, pad=0.02)
    cb.set_label("displacement magnitude")
    fig.tight_layout()
    fig.savefig(fname + ".pdf", dpi=200)
    fig.savefig(fname + ".png", dpi=200)
