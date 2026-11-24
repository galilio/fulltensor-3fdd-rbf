# Full-tensor diffusion on a non-convex domain: FE, FV and 3FDD-RBF

Comparison of three methods for `-div(K grad u) = f` with a **full** diffusion tensor on an
L-shaped (non-convex) domain, all driven by automatic residual-based adaptivity:

- **3FDD-RBF**: three-field domain decomposition with local radial basis functions (r^5 + P2), mesh-free;
- **FE P1**: Galerkin finite elements on adaptive Delaunay meshes;
- **FV TPFA**: two-point flux finite volumes on the same meshes.

## Contents

| Path | Content |
|---|---|
| `code/01_3FDD_Rectangles/` | C99 implementation of adaptive 3FDD-RBF (three-field domain decomposition) |
| `code/03_EF_P1/` | C99 implementation of adaptive Galerkin FE P1 |
| `code/04_VF_TPFA/` | C99 implementation of adaptive FV TPFA |
| `code/README.md` | Compilation and usage instructions for all three solvers |

## Quick start (compile and run)

Requires `gcc` and `make` (`gnuplot` is optional, for live display `gui=1` and GIF export `-o`).

```sh
# 1. 3FDD-RBF (Rectangles)
cd code/01_3FDD_Rectangles && make
./tfdd_rect_lshape gui=0 ex=3              # Console mode (Example 3, PHS r^5 + P2)
./tfdd_rect_lshape gui=1 ex=3              # Live gnuplot display
./tfdd_rect_lshape gui=0 ex=3 rbf=wendland # Adaptive Wendland C^4 RBF

# 2. Galerkin FE P1
cd ../03_EF_P1 && make
./ef_p1_lshape gui=0 ex=3

# 3. FV TPFA
cd ../04_VF_TPFA && make
./vf_tpfa_lshape gui=0 ex=0                # Isotropic check (K = I)
./vf_tpfa_lshape gui=0 ex=3                # Full tensor (O(1) stagnation)
```

See [`code/README.md`](code/README.md) for all options (`ex=`, `gui=`, `rbf=`, `-o file.gif`, `pts=1`).

Author: galilio.
