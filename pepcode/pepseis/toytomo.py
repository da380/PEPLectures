"""A toy two-dimensional delay-time tomography problem with straight rays, used for the
Lecture 19 figures.  The unit square is divided into n x n cells; the model is the fractional
slowness perturbation in each cell; the data are delay times along straight rays between
randomly placed sources and receivers, with Gaussian errors.  Regularised least squares with
norm damping (B = I) is used throughout, so that the posterior covariance of the Bayesian
section is (A^T C^-1 A + lambda I)^-1.
"""
import numpy as np


class ToyTomography:
    def __init__(self, n=24, n_src=20, n_rec=14, seed=6, noise=0.02):
        self.n = n
        rng = np.random.default_rng(seed)
        self.rng = rng
        xc = (np.arange(n) + 0.5) / n
        # input model: smooth random field, +/- 0.4 per cent
        kx = 2 * np.pi * np.fft.fftfreq(n, 1 / n); KX, KY = np.meshgrid(kx, kx, indexing="ij")
        field = np.real(np.fft.ifft2(np.fft.fft2(rng.standard_normal((n, n))) * np.exp(-(KX ** 2 + KY ** 2) * 0.12 ** 2 / 4)))
        self.m_true = 0.004 * field / np.abs(field).max()
        # geometry and path-length matrix
        self.src = rng.uniform(0.05, 0.95, (n_src, 2)); self.rec = rng.uniform(0.05, 0.95, (n_rec, 2))
        self.rays = [(s, r) for s in self.src for r in self.rec]
        A = np.zeros((len(self.rays), n * n))
        for i, (s, r) in enumerate(self.rays):
            L = np.linalg.norm(r - s); t = np.linspace(0, 1, 2000)
            pts = s[None, :] + t[:, None] * (r - s)[None, :]
            cells = np.clip((pts * n).astype(int), 0, n - 1)
            np.add.at(A[i], cells[:, 0] * n + cells[:, 1], L / len(t))
        self.A = A
        d_true = A @ self.m_true.ravel()
        self.sigma = noise * np.abs(d_true).max()
        self.d = d_true + self.sigma * rng.standard_normal(d_true.shape)
        self.AtA = A.T @ A / self.sigma ** 2
        self.Atd = A.T @ self.d / self.sigma ** 2
        self.coverage = (A != 0).any(axis=0).reshape(n, n)

    # --- inversion ---------------------------------------------------------------------------
    def solve(self, lam, d=None):
        rhs = self.Atd if d is None else self.A.T @ d / self.sigma ** 2
        return np.linalg.solve(self.AtA + lam * np.eye(self.n ** 2), rhs).reshape(self.n, self.n)

    def chi2_per_datum(self, m):
        return np.sum((self.A @ m.ravel() - self.d) ** 2) / self.sigma ** 2 / len(self.d)

    def lambda_for_chi2(self, target=1.0, lo=1e-6, hi=1e12):
        for _ in range(50):
            lam = np.sqrt(lo * hi)
            lo, hi = (lam, hi) if self.chi2_per_datum(self.solve(lam)) < target else (lo, lam)
        return lam

    def posterior_std(self, lam):
        C = np.linalg.inv(self.AtA + lam * np.eye(self.n ** 2))
        return np.sqrt(np.diag(C)).reshape(self.n, self.n)

    def prior_covariance(self, std=0.003, corr_len=0.15):
        """Gaussian prior covariance between cells: std^2 exp(-r^2 / 2 corr_len^2)."""
        xc = (np.arange(self.n) + 0.5) / self.n
        X, Y = np.meshgrid(xc, xc, indexing="ij")
        P = np.column_stack([X.ravel(), Y.ravel()])
        r2 = ((P[:, None, :] - P[None, :, :]) ** 2).sum(-1)
        return std ** 2 * np.exp(-r2 / (2 * corr_len ** 2))

    def posterior(self, Cm):
        """Posterior mean and covariance for a Gaussian prior N(0, Cm) and data covariance
        sigma^2 I, computed in data space: m_p = Cm A^T (A Cm A^T + C)^-1 d."""
        A = self.A
        S = A @ Cm @ A.T + self.sigma ** 2 * np.eye(A.shape[0])
        K = Cm @ A.T @ np.linalg.inv(S)
        m_p = (K @ self.d).reshape(self.n, self.n)
        C_p = Cm - K @ A @ Cm
        return m_p, np.sqrt(np.clip(np.diag(C_p), 0, None)).reshape(self.n, self.n)

    def null_space_component(self, m):
        """Projection of m onto the null space of A (m minus its minimum-norm reconstruction)."""
        m_range = np.linalg.lstsq(self.A, self.A @ m.ravel(), rcond=None)[0]
        return (m.ravel() - m_range).reshape(self.n, self.n)

    def checkerboard(self, period=8, amp=0.004):
        i, j = np.meshgrid(np.arange(self.n), np.arange(self.n), indexing="ij")
        return amp * np.where(((i // period) + (j // period)) % 2 == 0, 1.0, -1.0)

    # --- drawing -----------------------------------------------------------------------------
    def draw(self, ax, m, vmax=None, rays=True, cmap="RdBu_r", title=None, vmin=None):
        n = self.n
        e = np.linspace(0, 1, n + 1)
        v = vmax if vmax is not None else np.abs(m).max()
        lo = vmin if vmin is not None else (-v if cmap == "RdBu_r" else 0)
        im = ax.pcolormesh(e, e, m.T, cmap=cmap, vmin=lo, vmax=v, shading="auto", rasterized=True)
        if rays:
            for s, r in self.rays:
                ax.plot([s[0], r[0]], [s[1], r[1]], "k", lw=0.2, alpha=0.3)
        ax.plot(self.src[:, 0], self.src[:, 1], "k*", ms=6); ax.plot(self.rec[:, 0], self.rec[:, 1], "kv", ms=5)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        if title:
            ax.set_title(title, fontsize=11)
        return im


def figure_style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})
    return plt
