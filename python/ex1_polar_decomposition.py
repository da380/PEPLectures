"""Figure for Example 2 of the worked examples for Lecture 12.

Polar decomposition F = R U = V R of the plane deformation gradient
F = [[2, 1], [0, 1]], drawn in the spirit of Fig. 1 of Lecture 13: the unit
circle (centre of the upper row) is mapped by U to an ellipse whose principal
axes lie along the eigenvectors of U (right), by R to a rotated circle (left),
and by F to the final ellipse (bottom), whose axes lie along the eigenvectors
of V = R U R^T.  The dash-dotted lines are the eigenvectors of U and their
images under each map.  The matrices are verified numerically before plotting.

Output: examples/figures_ex/ex1_polar.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "examples", "figures_ex", "ex1_polar.png")

# ---------------------------------------------------------------- matrices
F = np.array([[2.0, 1.0], [0.0, 1.0]])
C = F.T @ F
lam, Q = np.linalg.eigh(C)                 # ascending eigenvalues
U = Q @ np.diag(np.sqrt(lam)) @ Q.T        # U = sqrt(C) by eigendecomposition
R = F @ np.linalg.inv(U)
V = R @ U @ R.T

# closed forms from the worked solution
U_exact = np.array([[6.0, 2.0], [2.0, 4.0]]) / np.sqrt(10.0)
R_exact = np.array([[3.0, 1.0], [-1.0, 3.0]]) / np.sqrt(10.0)
V_exact = np.array([[7.0, 1.0], [1.0, 3.0]]) / np.sqrt(10.0)

def check(name, ok, val=""):
    print(f"{'OK  ' if ok else 'FAIL'} {name} {val}")
    assert ok, name

tol = 1e-12
check("eigenvalues of C = 3 -/+ sqrt5", np.allclose(lam, [3 - np.sqrt(5), 3 + np.sqrt(5)]), lam)
check("U = sqrt(C) closed form", np.allclose(U, U_exact, atol=tol))
check("U^2 = C", np.allclose(U @ U, C, atol=tol))
check("R orthogonal", np.allclose(R.T @ R, np.eye(2), atol=tol))
check("det R = 1", abs(np.linalg.det(R) - 1) < tol, np.linalg.det(R))
check("R closed form", np.allclose(R, R_exact, atol=tol))
check("V = R U R^T closed form", np.allclose(V, V_exact, atol=tol))
check("F = R U", np.allclose(R @ U, F, atol=tol))
check("F = V R", np.allclose(V @ R, F, atol=tol))
theta = np.arctan2(R[1, 0], R[0, 0])
check("tan(theta) = -1/3", abs(np.tan(theta) + 1.0 / 3.0) < tol, f"theta = {np.degrees(theta):.3f} deg")

# eigenvectors of U (columns of Q, ascending eigenvalue order)
q_minus, q_plus = Q[:, 0], Q[:, 1]

# ---------------------------------------------------------------- drawing
phi = np.linspace(0.0, 2.0 * np.pi, 400)
circle = np.vstack([np.cos(phi), np.sin(phi)])            # 2 x N
axes_ref = [np.outer(q_plus, [-1.0, 1.0]), np.outer(q_minus, [-1.0, 1.0])]  # segments

centres = {
    "ref": np.array([0.0, 0.0]),
    "R":   np.array([-6.5, 0.0]),
    "U":   np.array([6.5, 0.0]),
    "F":   np.array([0.0, -6.0]),
}
maps = {"ref": np.eye(2), "R": R, "U": U, "F": F}

fig, ax = plt.subplots(figsize=(9.0, 4.8), facecolor="white")
ax.set_facecolor("white")

for key, A in maps.items():
    c = centres[key]
    curve = A @ circle
    ax.plot(curve[0] + c[0], curve[1] + c[1], color="black", lw=1.6)
    for seg in axes_ref:
        s = A @ seg
        ax.plot(s[0] + c[0], s[1] + c[1], color="0.55", lw=0.9, ls="-.")

# arrows between the shapes; (start, end, label, label offset)
def arrow(p, q, label, off, shrink_p, shrink_q):
    p, q = np.asarray(p, float), np.asarray(q, float)
    d = q - p
    d /= np.linalg.norm(d)
    p2 = p + shrink_p * d
    q2 = q - shrink_q * d
    ax.annotate("", xy=q2, xytext=p2,
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.6,
                                mutation_scale=16))
    m = 0.5 * (p2 + q2) + np.asarray(off, float)
    ax.text(m[0], m[1], label, fontsize=17, fontweight="bold",
            ha="center", va="center")

arrow(centres["ref"], centres["R"], r"$\mathbf{R}$", (0.0, 0.7), 1.4, 1.5)
arrow(centres["ref"], centres["U"], r"$\mathbf{U}$", (0.0, 0.7), 1.4, 2.6)
arrow(centres["ref"], centres["F"], r"$\mathbf{F}$", (-0.7, 0.0), 1.4, 1.5)
arrow(centres["R"], centres["F"], r"$\mathbf{V}$", (-0.8, -0.3), 1.6, 2.6)
arrow(centres["U"], centres["F"], r"$\mathbf{R}$", (0.9, -0.2), 2.0, 2.6)

ax.set_aspect("equal")
ax.set_xlim(-9.6, 9.6)
ax.set_ylim(-7.3, 2.4)
ax.axis("off")
fig.tight_layout(pad=0.2)
fig.savefig(OUT, dpi=200, facecolor="white")
print("wrote", os.path.abspath(OUT))
