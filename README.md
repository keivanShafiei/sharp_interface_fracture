# Sharp-Interface Variational Fracture Mechanics
## A Mesh-Conforming Finite Element Method Based on Configurational Forces

**Journal:** International Journal for Numerical Methods in Engineering  
**Status:** Preprint / Under Review  

---

## Overview

This repository implements a **sharp-interface finite element method** for quasi-static brittle fracture driven by configurational (material) forces. The method is motivated by the Francfort–Marigo variational principle but employs **local energy descent** — not global minimization — to track thermodynamically admissible crack evolution.

The discrete crack is represented as a union of mesh edges (conforming triangulation), enabling exact evaluation of the fracture surface energy `G_c · H¹(Γ_h)` without phase-field regularization.

---

## Energy Functional

The total potential energy is:

```
E(u, Γ) = ∫_{Ω\Γ} W(∇u) dx  +  G_c · H¹(Γ)  -  W_ext(u)
```

where:
- `W(∇u) = ½ ε : C : ε` is the linear elastic strain energy density (plane strain)
- `G_c` is the critical energy release rate [J/m²]
- `H¹(Γ)` is the 1-D Hausdorff measure (crack length)
- `W_ext` is the external work

The discrete crack set `Γ_h` is a union of element edges of a conforming triangulation `T_h`.

---

## Mathematical Assumptions

The following assumptions underpin the theoretical results. **All must be satisfied for the guarantees to hold.**

| Assumption | Statement | Status in Code |
|---|---|---|
| A1: Linear elasticity | `W = ½ ε:C:ε`, small strains | ✓ Enforced |
| A2: Plane strain | `ε_zz = 0` throughout | ✓ Enforced (flag in Material) |
| A3: Coercivity | `E_h(·,Γ_h)` strongly convex with constant `α > 0` | ✓ Holds for linear elasticity |
| A4: Shape regularity | `θ_min ≥ 22°` on all elements | ⚠ Enforced by remesher; not verified post-remesh |
| A5: Crack irreversibility | `Γ_h` only grows | ✓ CrackPath::add_segment only appends |
| A6: Quasi-static loading | Inertia neglected | ✓ No mass matrix assembled |

---

## What Has Been Proved

- **Energy monotonicity of accepted updates** (Proposition 3.1): Every accepted propagation step satisfies `E_h^{k+1} ≤ E_h^k + ε_tol · |E_h^k|` with `ε_tol = 10⁻⁶`.
- **Finite termination of backtracking** (Lemma 3.1): Backtracking terminates in at most `⌈log₂(β/β_min)⌉` halvings.
- **ε-stationary convergence** (Proposition 3.2): The sequence converges to a point where `|G_tip| ≤ G_c + θ_h`, where `θ_h = 0.01 G_c`.
- **O(h²) energy convergence** for prescribed smooth cracks on aligned meshes.
- **O(h) geometry convergence** for prescribed smooth cracks.

## What Has NOT Been Proved

- **Global minimality**: The algorithm finds local/metastable stationary states, not global minimizers of `E(u,Γ)`.
- **Γ-convergence under adaptive remeshing**: Static Γ-convergence holds for fixed aligned meshes; under adaptive remeshing (crack path not known a priori) this is unresolved.
- **Uniqueness of crack path**: Different load steps or mesh configurations may produce different paths.
- **Branching**: The current code supports only single-tip propagation. No branching logic is implemented.
- **3D extension**: All results are strictly 2D plane strain.

---

## Claimed Performance

| Benchmark | Speedup over AT2 Phase-Field | Condition |
|---|---|---|
| SENT (Mode-I) | 6.1–12.0× | Equivalent COD accuracy (5% tolerance) |
| TPB (Mixed-Mode) | 6.1–12.0× | Equivalent COD accuracy |

⚠ **Important:** The phase-field solver used for comparison is **not included** in this repository. The speedup figures cannot be independently reproduced without a reference AT2 implementation. See `LIMITATIONS.md`.

---

## Repository Structure

```
/src           → C++ source files
/include       → C++ header files
/tests         → Scientific verification tests
/benchmarks    → Benchmark scripts (SENT, TPB, Double-Notch)
/docs          → THEORY.md, VERIFICATION.md, LIMITATIONS.md
README.md
THEORY.md
VERIFICATION.md
LIMITATIONS.md
requirements.txt
reproduce_all_results.sh
```

---

## Requirements

```
C++14 or later
PETSc ≥ 3.18 (with GMRES + BoomerAMG)
Eigen ≥ 3.4
CGAL ≥ 5.4 (for constrained Delaunay remeshing)
Triangle (Shewchuk 1996) for initial mesh generation
GCC ≥ 9 or Clang ≥ 12
```

Reference hardware: Intel Core i7-9700K @ 3.6 GHz, 32 GB RAM, Ubuntu 22.04.2 LTS, GCC 11.3.0 `-O3 -march=native`

---

## How to Reproduce Paper Results

```bash
# Build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4

# Run SENT benchmark
./sharp_fracture --case SENT --mesh meshes/sent_h010.msh

# Run TPB benchmark
./sharp_fracture --case TPB --mesh meshes/tpb_h020.msh

# Run all and compare
bash reproduce_all_results.sh
```

See `VERIFICATION.md` for pass/fail criteria for each benchmark.

---

## Known Limitations Before Publication

See `LIMITATIONS.md` for a complete list. Critical items:

1. The elasticity tensor in `config_force.cpp` uses non-standard parameterization that differs from `material.hpp` by ~1.6% for ν=0.20.
2. Backtracking in `propagate_with_backtracking()` uses a **linear energy approximation** rather than a full FEM re-solve.
3. The Double-Notch benchmark is specified in the paper but not yet implemented in code.
4. OpenMP parallelization of Phase 2 (configurational force evaluation) is declared in the paper but not implemented.
5. The speedup comparison requires a phase-field reference implementation not included here.

---

## Citation

```bibtex
@article{SharpInterface2024,
  title   = {Sharp-Interface Variational Fracture Mechanics: A Mesh-Conforming 
             Finite Element Method Based on Configurational Forces},
  journal = {International Journal for Numerical Methods in Engineering},
  year    = {2024},
  doi     = {10.1002/0000}
}
```
