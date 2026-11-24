# Full-tensor diffusion on a non-convex domain — source code

This directory contains the three C99 programs used in the article
*"Anisotropic full-tensor diffusion on a non-convex domain: Finite Elements, Finite
Volumes and a three-field RBF domain decomposition, compared under automatic
residual-based adaptivity"* (`../article/en/`).

All three programs solve the same problem

    -div( K(x) grad u ) = f   in Omega,      u = g_D   on dOmega,

on the L-shaped domain `Omega = [0,1]x[0,2] U [1,2]x[0,1]` (re-entrant corner at `(1,1)`),
with a **full** symmetric positive definite tensor `K` (`K_xy != 0`). Each program runs
its own **automatic residual-based adaptive loop** and compares the numerical solution with
the exact one.

| Directory | Program | Method | Adaptivity |
|---|---|---|---|
| `01_3FDD_Rectangles/` | `tfdd_rect_lshape` | Three-field domain decomposition + local RBF (3FDD–RBF), mesh-free | residual subsampling (ARS): new centres + splitting of rectangles |
| `03_EF_P1/` | `ef_p1_lshape` | P1 Finite Elements (Galerkin) | residual estimator + Dörfler marking + Delaunay remeshing |
| `04_VF_TPFA/` | `vf_tpfa_lshape` | Two-point flux Finite Volumes (Voronoi cells) | same as FE |

---

## 1. Requirements and compilation

- A C99 compiler (`gcc`) and the math library. No other dependency.
- Optional: `gnuplot`, for the live display and the GIF animations.

```sh
cd 01_3FDD_Rectangles && make    # -> tfdd_rect_lshape
cd ../03_EF_P1         && make    # -> ef_p1_lshape
cd ../04_VF_TPFA       && make    # -> vf_tpfa_lshape
```

Each `Makefile` compiles with `gcc -O3 -std=c99 -Wall -Wextra -Wpedantic ... -lm`. Each program is
a single source file. `make clean` removes the executable and the generated files
(`*.dat`, and in most folders also `*.pdf`, `*.gif` and LaTeX by-products).

---

## 2. Test cases

The examples are hard-coded and selected with `ex=`:

| `ex=` | Tensor `K` | Exact solution | Difficulty |
|---|---|---|---|
| 1 | `[[2,-1],[-1,2]]` | `sin(pi x/2) sin(pi y/2) + 0.25 x y` | constant full tensor |
| 2 | `I` for `x<1`, `[[8,-7],[-7,20]]` for `x>1` | `exp(0.5(x+y))` | full tensor, discontinuous at `x=1` (20:1 contrast) |
| 3 | `R(30°) diag(50,1) R(30°)^T` | `exp(-35((x-0.92)^2+(y-0.92)^2)) + 0.2(x+y)` | anisotropy 50:1, peak at the re-entrant corner |
| 0 | `I` (FV only) | same as Ex. 1 | isotropic control case: checks that TPFA is correct when `K` is diagonal |

Without `ex=`, a program runs all its examples one after the other.

---

## 3. Running the programs

Typical use (no graphics, one example):

```sh
./tfdd_rect_lshape gui=0 ex=3
./ef_p1_lshape     gui=0 ex=3
./vf_tpfa_lshape   gui=0 ex=3
```

### Options common to the three programs

| Option | Meaning | Default |
|---|---|---|
| `ex=k` | run only example `k` (`1`, `2`, `3`; `0` also for FV) | all examples |
| `gui=0` / `gui=1` | `gui=1` opens a live gnuplot window updated at each iteration | `1` |
| `delay=ms` | pause between frames of the live display / GIF (milliseconds) | `650` |
| `-o file.gif` | write an animated GIF of the adaptive iterations (needs gnuplot) | — |
| `pts=1` | at the end of each example, print the table `u_h - u_exact` at **every** point | `0` |

`tfdd_rect_lshape` also accepts `out=file.gif` as a synonym of `-o file.gif`.

### Options of `tfdd_rect_lshape` (3FDD–RBF)

| Option | Meaning | Default |
|---|---|---|
| `n=N` | `N x N` regular centres in every subdomain (initial and after each split) | `6` |
| `rbf=...` | local RBF, always augmented with P2 polynomials: `phs` (= `phs5`, r^5), `phs3` (r^3), `tps` (r^4 log r), `gauss`, `mq` (multiquadric), `imq` (inverse MQ), `wendland` | `phs` |
| `kappa=c` | shape parameter of `gauss`, `mq`, `imq`: width `s = kappa * h_k`, `h_k = 1/(sqrt(N_k)-1)` in normalised coordinates | `1` |
| `nb=b` | Legendre modes of the multiplier `lambda_k` per edge (`nc = b-1` modes for the trace `psi`) | `4` |
| `mmax=M` | maximum number of subdomains | `128` |
| `alpha=a`, `alpha=adapt` | Wendland only: support factor, or its automatic residual-based tuning | `18` |

The article uses the defaults: `rbf=phs n=6 nb=4 mmax=128`. The RBF/`n` comparison of the article
(Table 1) is obtained by looping over `rbf=`, `kappa=` and `n=`.

---

## 4. What is printed and written

### Terminal output (messages are in French)

One line per adaptive iteration, then a summary line `Bilan`.

- **3FDD**: `Iter`, `Nb_Dom` (number of subdomains M), `nu(m..M)` (min..max centres per
  subdomain), `N_u` (total centres), `N_Sig` (unknowns on the skeleton Sigma), `e_L2`,
  `e_max` (max error at the centres), `e_Linf_q` (max error at the quadrature points),
  `saut_Sig` (max jump of `u_h` across interfaces), `t(s)` (CPU time of one solve).
- **FE / FV**: `Iter`, `Nb_Nds`/`Nb_Vol` (nodes / volumes), `Nb_Tri` (triangles), `angmin`
  (smallest angle of the mesh, degrees), `obtus` (number of obtuse triangles), `iterK`
  (BiCGSTAB iterations), `e_L2`, `e_max` (max error at the nodes), `eta` (global error
  estimator), `t(s)`.

With `pts=1`, the final table lists for every point: index, `x`, `y`, `u_h`, `u_exact` and
`u_h - u_exact` (3FDD also prints the subdomain index `k`).

### Files (written in the current directory, `k` = example number)

| File | Content |
|---|---|
| `exk_history.dat` | one line per iteration. 3FDD: `iter M N N_Sigma e_L2 e_Linf_q e_energy jump t alpha_min alpha_max e_max`. FE/FV: `iter N N_tri e_L2 e_Linf e_energy eta t e_max angle_min n_obtuse nnz_per_row` |
| `exk_erreur_points.dat` | pointwise error at the final iteration. 3FDD: `k i x y u_h u_exact u_h-u_exact`. FE/FV: `i x y u_h u_exact u_h-u_exact` |
| `exk_init_rects.dat`, `exk_final_rects.dat` | 3FDD: subdomain rectangles (closed polylines, gnuplot format) |
| `exk_init_nodes.dat`, `exk_final_nodes.dat` | 3FDD: centres `x y u_h k` |
| `exk_final_tri.dat`, `exk_final_sol.dat` | FE/FV: final triangles (closed polylines) and nodal solution |
| `live_*.dat` | temporary files of the live display |

Each folder contains a `make_figs.gp` that draws its own figures from these files (`gnuplot make_figs.gp`).

---

## 5. The three methods in the code

### 5.1 `01_3FDD_Rectangles/tfdd_rect_lshape.c` — 3FDD–RBF (mesh-free)

**Discretisation.** The domain is split into rectangles `Omega_k`. On each rectangle, the
solution `u_k` is expanded in a local cardinal RBF basis built on the centres of `Omega_k`.
The Gram matrix `G_k = [Phi_k P_k; P_k^T 0]` is factorised once (LU). The normal flux
`lambda_k` on each edge and the trace `psi` on the skeleton `Sigma` are Legendre expansions
(`nb = 4`, `nc = 3` modes; 2 modes on short edges).

**Assembly** (`assemble_local`). The local stiffness matrix is
`(A_k)_ij = sum_q w_q (K grad phi_i) . grad phi_j`. It uses Gauss 5x5 quadrature on `s x s`
sub-rectangles, with `s = max(2, ceil(sqrt(N_k)/3))`. The program also builds the coupling
blocks `B_k` (centres/multipliers) and `C_k` (multipliers/trace). The local saddle-point
block `K_k = [A_k B_k; B_k^T 0]` is dense but small (at most a few tens of centres); it is
factorised by LU with partial pivoting.

**Global solve** (`solve_rect_3fdd`). The local blocks are condensed into a symmetric
positive definite Schur complement `S psi = r_Sigma` on the skeleton, which is solved by
LU. `backsub_local` then recovers `u_k`, `lambda_k` and the RBF coefficients
independently in each subdomain.

**Initial decomposition** (`init_lshape_2rects`). Minimal decomposition of the L-shape:
`R1 = [0,1]x[0,2]` and `R2 = [1,2]x[0,1]`, each with `n x n` regular centres
(`rect_seed_points`). Edge points are added only where the inf-sup condition requires them
(hanging nodes).

**Adaptivity** (`adapt_rect_and_nodes`, Algorithm 2 of the article).
1. Strong residual `D_k |f + div(K grad u_k)|` on an 8x8 test grid in each rectangle
   (exact second derivatives of the RBF).
2. Level 1: up to 12 new centres where the residual is at least `theta_r` times the global
   maximum (`theta_r = 0.25` for Ex. 1–2, `0.15` for Ex. 3), with a minimum distance.
3. Level 2: a rectangle with more than `n_split = max(44, n^2+8)` centres is cut in two
   along its longest side.
4. 2:1 balancing: a rectangle touching a neighbour more than 2.5 times smaller is cut.
5. Interfaces are rebuilt (non-conforming sub-segments, `build_rect_interfaces`).

### 5.2 `03_EF_P1/ef_p1_lshape.c` — P1 Finite Elements

- **Mesh**: Bowyer–Watson Delaunay triangulation of a point cloud (initial grid `h = 1/6`, 133 nodes); triangles outside the L-shape are removed.
- **Assembly** (`solve_ef_p1`): element stiffness `|T| (K_T grad l_c) . grad l_a` with `K_T` averaged by a 7-point Dunavant rule; Dirichlet nodes are eliminated, so the matrix is symmetric positive definite (about 6.5 non-zeros per row).
- **Solver**: CSR storage, ILU(0)-preconditioned BiCGSTAB, relative tolerance `1e-12`.
- **Estimator** (`estimate_eta`): `eta_T^2 = h_T^2 ||f||^2_T + 1/2 sum_e |e|^2 [[ (K grad u_h).n ]]^2`.
- **Refinement** (`refine_dorfler`): Dörfler marking (`theta = 0.5`), insertion of edge midpoints (minimum distance `0.35 |e|`), new Delaunay triangulation.
- **Mesh quality** (`mesh_quality`): smallest angle and number of obtuse triangles, printed at every iteration.

### 5.3 `04_VF_TPFA/vf_tpfa_lshape.c` — two-point flux Finite Volumes

- **Same mesh, solver, estimator and refinement** as the FE code.
- **Assembly** (`solve_tpfa`): Voronoi control volumes. For each edge `(i,j)` of a triangle with opposite angle `theta_k`, the transmissibility is `tau = 1/2 cot(theta_k) n^T K n`, with `n = (x_j - x_i)/|x_j - x_i|`.
- **Expected behaviour**: only `n^T K n` is used, so the scheme is **inconsistent** for a full tensor. The error stagnates on Ex. 1–3 and converges only on the isotropic control case Ex. 0. This is the point made in the article. Obtuse angles also make some transmissibilities negative.

---

## 6. Reproducing the article results

```sh
for d in 01_3FDD_Rectangles 03_EF_P1 04_VF_TPFA; do (cd $d && make); done
(cd 01_3FDD_Rectangles && ./tfdd_rect_lshape gui=0 rbf=phs n=6)
(cd 03_EF_P1           && ./ef_p1_lshape gui=0)
(cd 04_VF_TPFA         && ./vf_tpfa_lshape gui=0)
```

Each program writes `ex1..3_history.dat` and `ex1..3_erreur_points.dat` in its own folder.
These files contain the values of Tables 2–5 of the article (`../article/en/tables_results.tex`
is generated from such files by `../article/en/make_tables_en.py`; the figures by
`../article/en/make_figs_en.gp`). Expected final values:

| Ex. | 3FDD (M, N, e_L2) | FE (N, e_L2) | FV (N, e_L2) |
|---|---|---|---|
| 1 | 47, 1697, 1.07e-4 | 2347, 7.22e-4 | 3389, 4.23e-2 |
| 2 | 36, 1316, 1.48e-5 | 1567, 6.85e-4 | 484, 2.30e-2 |
| 3 | 113, 4263, 3.27e-4 | 2369, 1.69e-3 | 1100, 1.26e-1 |

A full run of the three examples takes less than one second of CPU per program.
