"""
fem_convergence_Lshape.py

Same P1 solver as fem_convergence.py, on the L-shaped domain
    Omega = (-1,1)^2 \\ ([0,1]x[-1,0])
with the classical reentrant-corner solution
    u(r,theta) = r^(2/3) sin(2 theta/3),   theta in [0, 3 pi/2],
which is harmonic (so f = 0) and vanishes on the two edges meeting at the
reentrant corner.  Because u is only in H^{1+2/3-eps}, uniform refinement
gives the SUBOPTIMAL rates
    |u-U|_{H1} ~ h^{2/3} ~ Np^{-1/3},   ||u-U||_{L2} ~ h^{4/3} ~ Np^{-2/3}
(the L2 order is 4/3, not 2, because the dual problem is equally singular).

Requires numpy and scipy.  Run:  python3 fem_convergence_Lshape.py
"""
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

def theta(x):
    t = np.arctan2(x[1], x[0])
    return t + 2*np.pi if t < 0 else t

def u_ex(x):
    r = np.hypot(x[0], x[1])
    return 0.0 if r == 0 else r**(2.0/3.0) * np.sin(2.0/3.0 * theta(x))

def grad_u(x):
    r = np.hypot(x[0], x[1])
    if r == 0:
        return np.array([0.0, 0.0])
    th = theta(x); c = (2.0/3.0) * r**(-1.0/3.0)
    return np.array([-c*np.sin(th/3.0), c*np.cos(th/3.0)])

Ghat = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]]).T

def build_mesh(N):
    """Uniform mesh of the L-shape; keep cells whose centre is inside Omega."""
    h = 1.0 / N
    npts = 2*N + 1
    grid = -np.ones((npts, npts), dtype=int)   # global index of each grid point, -1 if unused
    coords, elem = [], []
    def gid(I, J):                              # lazily number the grid points we use
        if grid[I, J] < 0:
            grid[I, J] = len(coords)
            coords.append([-1.0 + J*h, -1.0 + I*h])
        return grid[I, J]
    for I in range(2*N):
        for J in range(2*N):
            cx, cy = -1.0 + (J+0.5)*h, -1.0 + (I+0.5)*h
            if cx > 0 and cy < 0:               # skip the missing quadrant
                continue
            a, b = gid(I, J),   gid(I, J+1)
            c, d = gid(I+1, J), gid(I+1, J+1)
            elem.append([a, d, c]); elem.append([a, b, d])
    coord = np.array(coords); elem = np.array(elem)
    # boundary vertices = endpoints of edges belonging to exactly one triangle
    from collections import Counter
    ec = Counter()
    for t in elem:
        for e in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])]:
            ec[tuple(sorted(e))] += 1
    bnd = sorted({v for e, k in ec.items() if k == 1 for v in e})
    return coord, elem, np.array(bnd)

def solve(N):
    coord, elem, diri = build_mesh(N)
    Np = coord.shape[0]
    A = lil_matrix((Np, Np)); fh = np.zeros(Np)   # f = 0
    for v in elem:
        P = coord[v]
        B = np.array([P[1]-P[0], P[2]-P[0]]).T
        area = abs(np.linalg.det(B)) / 2.0
        Binv = np.linalg.inv(B)
        Aloc = area * (Ghat.T @ (Binv @ Binv.T) @ Ghat)
        for a_ in range(3):
            for b_ in range(3):
                A[v[a_], v[b_]] += Aloc[a_, b_]
    for d in diri:
        A.rows[d] = [d]; A.data[d] = [1.0]; fh[d] = u_ex(coord[d])
    return coord, elem, spsolve(A.tocsr(), fh)

def errors(coord, elem, uh):
    H1sq = L2sq = 0.0
    for v in elem:
        P = coord[v]; B = np.array([P[1]-P[0], P[2]-P[0]]).T
        area = abs(np.linalg.det(B)) / 2.0
        gradUT = np.linalg.solve(B.T, Ghat @ uh[v])
        mids = [(P[0]+P[1])/2, (P[1]+P[2])/2, (P[2]+P[0])/2]
        uef = [(uh[v[0]]+uh[v[1]])/2, (uh[v[1]]+uh[v[2]])/2, (uh[v[2]]+uh[v[0]])/2]
        for m in mids:
            d = gradUT - grad_u(m); H1sq += (d @ d)/3.0 * area
        for m, ue in zip(mids, uef):
            L2sq += (u_ex(m) - ue)**2 / 3.0 * area
    return np.sqrt(H1sq), np.sqrt(L2sq)

if __name__ == "__main__":
    print(f"{'N':>4} {'Np':>7} {'H1err':>12} {'rate':>5} {'L2err':>12} {'rate':>5}")
    prevH = prevL = None
    for N in [4, 8, 16, 32, 64, 128]:
        coord, elem, uh = solve(N)
        h1, l2 = errors(coord, elem, uh)
        rH = np.log(prevH/h1)/np.log(2) if prevH else float('nan')
        rL = np.log(prevL/l2)/np.log(2) if prevL else float('nan')
        print(f"{N:>4} {coord.shape[0]:>7} {h1:>12.4e} {rH:>5.2f} {l2:>12.4e} {rL:>5.2f}")
        prevH, prevL = h1, l2
