"""Construction of a cross-correlation delay-time sensitivity kernel in two dimensions with
the adjoint method, following Lecture 20.

An explosive point source radiates a P-wave that is recorded at a receiver.  The delay time
tau_bar of the P arrival is measured by cross-correlation (eq. 36 of the lecture), and its
functional derivative with respect to the P-wave speed alpha is
    K_alpha(x) = -2 rho alpha int_0^T (div u)(x, t) (div u')(x, t) dt,
which follows from eq. (32) with delta A_ijkl = delta lambda delta_ij delta_kl and
delta lambda = 2 rho alpha delta alpha (beta and rho held fixed).  The adjoint field u'
satisfies the elastic wave equation with terminal conditions at t = T, driven by the adjoint
source h'_i = -w(t) s_dot(t) nu_i delta(x - x_r) / N with N = int w s_dot^2 dt; it is computed
by running the solver forward in the reversed time tau = T - t with the time-reversed source.
"""
import numpy as np
from pepseis.fd2d_psv import PSVSolver


def ricker(t, t0, f0):
    a = (np.pi * f0 * (t - t0)) ** 2
    return (1.0 - 2.0 * a) * np.exp(-a)


class KernelExperiment:
    def __init__(self, h=1.0 / 200, alpha0=1.0, beta0=0.58, src=(0.4, 0.5), rec=(1.6, 0.5),
                 f0=4.0, t_end=2.0, snap_every=4, alpha_pert=None, pad=0.5):
        # the computational box is padded by `pad` on every side so that the absorbing
        # boundaries lie well outside the region that is plotted (x in [0, 2], z in [0, 1])
        self.pad = pad
        self.x = np.arange(-pad, 2.0 + pad, h); self.z = np.arange(-pad, 1.0 + pad, h)
        self.plot_x = (self.x >= 0) & (self.x <= 2.0); self.plot_z = (self.z >= 0) & (self.z <= 1.0)
        self.h, self.alpha0, self.beta0, self.f0, self.t_end = h, alpha0, beta0, f0, t_end
        self.src, self.rec = src, rec
        self.isrc = (int(round((src[0] + pad) / h)), int(round((src[1] + pad) / h)))
        self.irec = (int(round((rec[0] + pad) / h)), int(round((rec[1] + pad) / h)))
        self.snap_every = snap_every
        X, Z = np.meshgrid(self.x, self.z, indexing="ij")
        self.alpha = alpha0 * np.ones_like(X) if alpha_pert is None else alpha0 * (1 + alpha_pert(X, Z))
        self.beta = beta0 * np.ones_like(X); self.rho = np.ones_like(X)

    def _solver(self):
        return PSVSolver(self.x, self.z, self.alpha, self.beta, self.rho, dt=0.45 * self.h / self.alpha0, sponge_z=True)

    def forward(self):
        """Run the forward problem; store the receiver displacement and snapshots of div u."""
        s = self._solver(); dt = s.dt
        n = int(np.ceil(self.t_end / dt / self.snap_every)) * self.snap_every   # a whole number of snapshots
        self.t_end = n * dt
        self.dt, self.nt = dt, n
        self.t = np.arange(n + 1) * dt
        self.s_rec = np.zeros(n + 1)
        self.snap_t, self.snap_div, self.snap_ux = [], [], []
        for k in range(n):
            q = 0.1 * ricker(self.t[k], 0.3, self.f0)
            s.step(explosion=(self.isrc[0], self.isrc[1], q))
            self.s_rec[k + 1] = s.ux[self.irec]
            if (k + 1) % self.snap_every == 0:
                self.snap_t.append(s.t); self.snap_div.append(s.divergence().astype(np.float32))
                self.snap_ux.append(s.ux.astype(np.float32))
        self.snap_t = np.array(self.snap_t)
        return self.t, self.s_rec

    def window(self, half_width=0.2):
        """Cosine-tapered window around the P arrival at the receiver."""
        dist = np.hypot(self.rec[0] - self.src[0], self.rec[1] - self.src[1])
        tc = 0.3 + dist / self.alpha0
        w = 0.5 * (1 + np.cos(np.pi * (self.t - tc) / half_width))
        w[np.abs(self.t - tc) > half_width] = 0.0
        return w

    def adjoint(self):
        """Run the adjoint problem in reversed time; store snapshots of div u' on the same
        time levels as the forward snapshots (matched by t = T - tau)."""
        w = self.window(); sdot = np.gradient(self.s_rec, self.dt)
        N = np.sum(w * sdot ** 2) * self.dt
        hsrc = -w * sdot / N                      # adjoint source time function h'(t) (x-component)
        self.hsrc = hsrc
        s = self._solver(); n = self.nt; m = len(self.snap_t)
        self.snap_div_adj = [None] * m; self.snap_ux_adj = [None] * m
        for k in range(n):
            # adjoint step k advances tau from k dt to (k+1) dt, i.e. forward time from T - k dt to T - (k+1) dt
            f = hsrc[n - k]                       # h' at forward time T - k dt
            s.step(force=(self.irec[0], self.irec[1], f / (self.h ** 2), 0.0))   # point force: delta function -> 1/h^2
            j = n - (k + 1)                       # forward time index reached
            if j % self.snap_every == 0 and j > 0:
                i = j // self.snap_every - 1      # forward snapshot index with t = j dt
                self.snap_div_adj[i] = s.divergence().astype(np.float32)
                self.snap_ux_adj[i] = s.ux.astype(np.float32)
        return hsrc

    def trim(self, f):
        """Restrict a field on the padded grid to the plotting region."""
        return f[np.ix_(self.plot_x, self.plot_z)]

    def integrand(self, i):
        """Instantaneous contribution to K_alpha at forward snapshot i:
        -2 rho alpha div u(x, t_i) div u^dagger(x, T - t_i)."""
        return -2.0 * self.rho * self.alpha * self.snap_div[i] * self.snap_div_adj[i]

    def kernel(self, upto=None):
        """K_alpha accumulated over forward snapshots up to index `upto` (all if None)."""
        K = np.zeros((len(self.x), len(self.z)))
        idx = range(len(self.snap_t)) if upto is None else range(upto)
        for i in idx:
            if self.snap_div_adj[i] is not None:
                K += self.snap_div[i] * self.snap_div_adj[i]
        return -2.0 * self.rho * self.alpha * K * self.dt * self.snap_every

    def delay_time(self, s_obs):
        """Cross-correlation delay time tau_bar of s_obs relative to the synthetic s_rec,
        using the lecture's definition C(tau) = int s_obs(t - tau) s(t) dt."""
        w = self.window()
        s = w * self.s_rec; so = w * s_obs
        lags = np.arange(-200, 201)
        C = np.array([np.sum(np.roll(so, lag) * s) for lag in lags])   # roll by lag: s_obs(t - lag dt)
        k = np.argmax(C)
        # parabolic refinement
        if 0 < k < len(C) - 1:
            d = 0.5 * (C[k - 1] - C[k + 1]) / (C[k - 1] - 2 * C[k] + C[k + 1])
        else:
            d = 0.0
        return (lags[k] + d) * self.dt
