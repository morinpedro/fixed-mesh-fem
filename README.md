# AFEM — Fixed-Mesh P1 Finite Element Code

Minimal Matlab/Octave code for solving a linear second-order elliptic problem
with continuous piecewise-linear (P1) finite elements on a **fixed** triangular
mesh, plus Python scripts that reproduce the same method and run the convergence
studies. It accompanies Section 3 ("Computational Implementation on a Fixed
Mesh") of the short course *Adaptive Finite Element Methods for Elliptic PDEs*.

## The model problem

The code solves

$$-\nabla\cdot(a\,\nabla u) + c\,u = f \quad \text{in } \Omega, \qquad
u = g_D \text{ on } \Gamma_D, \qquad a\,\partial_n u = g_N \text{ on } \Gamma_N,$$

with `a > 0` and `c ≥ 0` piecewise constant. Assembly is done element by element
via the reference triangle and the affine map, the load vector uses the
three-point midpoint quadrature (exact for quadratics), Dirichlet conditions are
imposed by row replacement, and the sparse system is solved with a direct solver.

## Getting the code

```bash
git clone https://github.com/morinpedro/fixed-mesh-fem.git
cd fixed-mesh-fem
```

## Requirements

- **Matlab or GNU Octave** for the `.m` files (Octave ≥ 5 recommended; it is free).
- **Python 3** with **NumPy** and **SciPy** for the `.py` scripts:
  ```bash
  pip install numpy scipy
  ```

## Files

| File | Description |
|------|-------------|
| `fem.m` | Driver: assembles and solves the problem above, then plots the solution. |
| `gen_mesh_rectangle.m` | Uniform triangulation of a rectangle (all-Dirichlet); writes the four mesh files. |
| `gen_mesh_L_shape.m` | Uniform triangulation of the L-shaped domain. |
| `H1_err.m` | `H^1`-seminorm error against a known solution (midpoint quadrature). |
| `L2_err.m` | `L^2` error against a known solution. |
| `fem_convergence.py` | Python convergence study on the unit square (smooth solution). |
| `fem_convergence_Lshape.py` | Python convergence study on the L-shaped domain (reentrant-corner solution). |

The mesh generators write four plain-text files that `fem.m` reads:
`vertex_coordinates.txt` (one `x y` per vertex), `elem_vertices.txt` (three
vertex indices per triangle), `dirichlet.txt` (Dirichlet vertex indices), and
`neumann.txt` (Neumann boundary segments, one `i j` per line; may be absent).

## Quick start (Octave/Matlab)

From this folder, start Octave (or Matlab) and run:

```matlab
gen_mesh_rectangle(16)   % uniform 16x16 mesh of the unit square -> writes the .txt files
fem                      % assemble, solve, and plot the P1 solution
```

The problem data (`coef_a`, `coef_c`, and the functions `fc_f`, `fc_gD`,
`fc_gN`) are set at the top of `fem.m`; edit them to change the problem. The
default corresponds to the manufactured solution `u(x) = exp(-10|x|^2)`, for
which `f = -Δu` and `g_D = u` are already provided.

To measure the error on the current mesh (workspace variables `uh`,
`elem_vertices`, `vertex_coordinates` are left by `fem.m`):

```matlab
u  = @(x) exp(-10*(x(1)^2 + x(2)^2));
gu = @(x) -20*[x(1); x(2)]*exp(-10*(x(1)^2 + x(2)^2));
eH1 = H1_err(elem_vertices, vertex_coordinates, uh, gu)   % energy error
eL2 = L2_err(elem_vertices, vertex_coordinates, uh, u)    % L2 error
```

Refining (`gen_mesh_rectangle(32)`, `64`, ...) and repeating shows the rates
`|u-U|_{H^1} = O(h)` and `||u-U||_{L^2} = O(h^2)` on the square.

## Convergence studies (Python)

The two Python scripts reproduce the identical P1 method and print the error
table and observed rates automatically:

```bash
python3 fem_convergence.py          # smooth solution on the unit square
python3 fem_convergence_Lshape.py   # u = r^(2/3) sin(2θ/3) on the L-shape
```

### Expected output

`fem_convergence.py` (unit square, smooth solution) — the `rate` columns are the
observed orders in `h = 1/N`, converging to `1` (energy) and `2` (`L^2`):

```
   N      Np        H1err  rate        L2err  rate
   4      25   3.9344e-01   nan   3.2858e-02   nan
   8      81   2.1884e-01  0.85   9.8784e-03  1.73
  16     289   1.1214e-01  0.96   2.5768e-03  1.94
  32    1089   5.6416e-02  0.99   6.5137e-04  1.98
  64    4225   2.8251e-02  1.00   1.6330e-04  2.00
 128   16641   1.4131e-02  1.00   4.0854e-05  2.00
```

`fem_convergence_Lshape.py` (L-shape, reentrant corner) — the orders in `h` are
now `2/3` (energy) and `4/3` (`L^2`), i.e. `N_p^{-1/3}` and `N_p^{-2/3}`:

```
   N      Np        H1err  rate        L2err  rate
   4      65   1.8150e-01   nan   1.8391e-02   nan
   8     225   1.1671e-01  0.64   7.2432e-03  1.34
  16     833   7.4587e-02  0.65   2.8706e-03  1.34
  32    3201   4.7433e-02  0.65   1.1420e-03  1.33
  64   12545   3.0064e-02  0.66   4.5497e-04  1.33
 128   49665   1.9013e-02  0.66   1.8129e-04  1.33
```

**Summary:**

- **Unit square** (smooth `u`): energy error `∝ N_p^{-1/2}` (i.e. `O(h)`),
  `L^2` error `∝ N_p^{-1}` (i.e. `O(h^2)`) — optimal.
- **L-shape** (reentrant corner, `u ∈ H^{1+2/3-ε}`): energy error
  `∝ N_p^{-1/3}` and `L^2` error `∝ N_p^{-2/3}` — **suboptimal**, which is the
  motivation for adaptivity in the later lectures.

## Notes

- The mesh generators write their output files into the current working
  directory, so run them from the folder where you will run `fem.m`.

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
