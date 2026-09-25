"""Check the 2D delay-time kernel against a finite-difference perturbation."""
import numpy as np, time
from pepseis.kernel2d import KernelExperiment
t0 = time.time()
E = KernelExperiment(); E.forward(); E.adjoint(); K = E.kernel()
print("forward+adjoint in %.0f s; K range %.3g .. %.3g" % (time.time() - t0, K.min(), K.max()))
for centre in [(1.0, 0.5), (1.0, 0.62), (1.0, 0.75)]:
    pert = lambda X, Z, c=centre: 0.02 * np.exp(-((X - c[0]) ** 2 + (Z - c[1]) ** 2) / (2 * 0.05 ** 2))
    Ep = KernelExperiment(alpha_pert=pert); Ep.forward()
    tau_meas = Ep.delay_time(E.s_rec)      # synthetic computed in the perturbed model, reference field as "observed"
    X, Z = np.meshgrid(E.x, E.z, indexing="ij")
    dalpha = E.alpha0 * pert(X, Z)
    tau_pred = np.sum(K * dalpha) * E.h ** 2   # (integral over the whole padded box)
    print("blob at %s: measured change in tau_bar %+.4f  predicted from kernel %+.4f" % (centre, tau_meas, tau_pred))
