"""
fem_convergence.py

P1 finite element solver for  -div(a grad u) + c u = f  on the unit square,
with a full-Dirichlet boundary, written to mirror the Matlab/Octave code
(fem.m, gen_mesh_rectangle.m, H1_err.m, L2_err.m).

It runs a convergence study for the manufactured solution
    u(x) = exp(-10 |x|^2),   a = 1,  c = 0,
    f    = -Lap u = 20 exp(-10 |x|^2) (2 - 20 |x|^2),
and prints the H1-seminorm and L2 errors together with the observed rates.
Same reference-element assembly and 3-point midpoint quadrature as fem.m.

Requires numpy and scipy.  Run:  python3 fem_convergence.py
"""
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

# ---- data of the manufactured problem -------------------------------------
def u_ex(x):                 # exact solution
    return np.exp(-10.0 * (x[0]**2 + x[1]**2))

def grad_u(x):               # exact gradient (used by the H1 error)
    e = np.exp(-10.0 * (x[0]**2 + x[1]**2))
    return np.array([-20.0 * x[0] * e, -20.0 * x[1] * e])

def f_rhs(x):                # right-hand side  f = -Laplacian(u)
    r2 = x[0]**2 + x[1]**2
    return 20.0 * np.exp(-10.0 * r2) * (2.0 - 20.0 * r2)

# reference-element gradients (columns = grad phihat_i), as in fem.m
Ghat = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]]).T   # shape (2,3)

# ---- uniform NxN triangulation of the unit square -------------------------
def build_mesh(N):
    xs = np.linspace(0.0, 1.0, N + 1)
    coord = np.array([[xs[j], xs[i]] for i in range(N + 1) for j in range(N + 1)])
    elem = []
    for i in range(N):
        for j in range(N):
            a = i * (N + 1) + j
            b = i * (N + 1) + j + 1
            c = (i + 1) * (N + 1) + j
            d = (i + 1) * (N + 1) + j + 1
            elem.append([a, d, c])
            elem.append([d, a, b])
    elem = np.array(elem)
    diri = [k for k in range((N + 1)**2)
            if coord[k, 0] in (0.0, 1.0) or coord[k, 1] in (0.0, 1.0)]
    return coord, elem, np.array(diri)

# ---- assembly and solve (a = 1, c = 0) ------------------------------------
def solve(N):
    coord, elem, diri = build_mesh(N)
    Np = coord.shape[0]
    A = lil_matrix((Np, Np))
    fh = np.zeros(Np)
    for v in elem:
        P = coord[v]                                   # 3x2 vertex coordinates
        B = np.array([P[1] - P[0], P[2] - P[0]]).T     # Jacobian [v2-v1  v3-v1]
        area = abs(np.linalg.det(B)) / 2.0
        Binv = np.linalg.inv(B)
        Aloc = area * (Ghat.T @ (Binv @ Binv.T) @ Ghat)          # stiffness
        m12, m23, m31 = (P[0]+P[1])/2, (P[1]+P[2])/2, (P[2]+P[0])/2
        f12, f23, f31 = f_rhs(m12), f_rhs(m23), f_rhs(m31)
        floc = np.array([(f12+f31)/2, (f12+f23)/2, (f23+f31)/2]) * (area/3.0)
        for a_ in range(3):
            fh[v[a_]] += floc[a_]
            for b_ in range(3):
                A[v[a_], v[b_]] += Aloc[a_, b_]
    for d in diri:                                     # Dirichlet by row replacement
        A.rows[d] = [d]
        A.data[d] = [1.0]
        fh[d] = u_ex(coord[d])
    return coord, elem, spsolve(A.tocsr(), fh)

# ---- errors (3-point midpoint quadrature, as in H1_err.m / L2_err.m) -------
def errors(coord, elem, uh):
    H1sq = L2sq = 0.0
    for v in elem:
        P = coord[v]
        B = np.array([P[1] - P[0], P[2] - P[0]]).T
        area = abs(np.linalg.det(B)) / 2.0
        gradUT = np.linalg.solve(B.T, Ghat @ uh[v])    # constant grad(U) on T
        mids = [(P[0]+P[1])/2, (P[1]+P[2])/2, (P[2]+P[0])/2]
        uef = [(uh[v[0]]+uh[v[1]])/2, (uh[v[1]]+uh[v[2]])/2, (uh[v[2]]+uh[v[0]])/2]
        for m in mids:
            d = gradUT - grad_u(m)
            H1sq += (d @ d) / 3.0 * area
        for m, ue in zip(mids, uef):
            L2sq += (u_ex(m) - ue)**2 / 3.0 * area
    return np.sqrt(H1sq), np.sqrt(L2sq)

# ---- convergence study ----------------------------------------------------
if __name__ == "__main__":
    print(f"{'N':>4} {'Np':>7} {'H1err':>12} {'rate':>5} {'L2err':>12} {'rate':>5}")
    prevH = prevL = None
    for N in [4, 8, 16, 32, 64, 128]:
        coord, elem, uh = solve(N)
        h1, l2 = errors(coord, elem, uh)
        rH = np.log(prevH / h1) / np.log(2) if prevH else float('nan')
        rL = np.log(prevL / l2) / np.log(2) if prevL else float('nan')
        print(f"{N:>4} {coord.shape[0]:>7} {h1:>12.4e} {rH:>5.2f} {l2:>12.4e} {rL:>5.2f}")
        prevH, prevL = h1, l2
