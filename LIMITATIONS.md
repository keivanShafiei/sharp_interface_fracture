# LIMITATIONS.md — Known Limitations, Assumptions, and Reviewer Risk Areas

This file is required for scientific transparency. Every limitation listed here is either inherent to the method or is a gap between the paper's description and the provided implementation. Reviewers, users, and contributors must read this file before drawing conclusions from simulation results.

---

## I. Theoretical Limitations (Inherent to the Method)

### L-T1: No Global Minimality Guarantee

**What the code does:** Seeks local/metastable stationary points of `E_h(u_h, Γ_h)` via a staggered descent scheme.

**What it does NOT do:** Minimize `E_h` globally over all admissible crack sets `K_h`.

**Implication:** The algorithm can miss lower-energy crack configurations (branching, path selection in non-convex geometries). Results represent one physically plausible metastable state, not the energetically optimal crack path.

**Paper statement:** Explicitly acknowledged in Section 1.3 and Section 3.3.

### L-T2: No Branching Support

**What the code does:** Tracks a single crack tip and propagates it in one direction per step.

**What it does NOT do:** Detect or simulate crack branching, coalescence, or multiple simultaneous tips.

**Implication:** The method cannot reproduce branching phenomena observed in experiments at high loading rates or in geometries with competing defects. The Double-Notch benchmark (path competition) is partially relevant here: the code does not implement it.

### L-T3: Γ-Convergence Under Adaptive Remeshing is Unresolved

**Status:** Static Γ-convergence (fixed mesh, prescribed crack) is established. For evolving cracks with adaptive remeshing, Γ-convergence as h→0 remains an open mathematical problem.

**Implication:** Mesh refinement does not guarantee convergence to the continuum solution for evolving crack paths. Results must be interpreted as mesh-dependent approximations.

### L-T4: Local Descent Does Not Imply Physical Uniqueness

The algorithm is consistent with the **local variational principle** of Larsen (2024) but is not a unique physical predictor. Different initial conditions or load histories may produce different crack paths that are all locally stationary.

---

## II. Implementation Gaps (Discrepancies Between Paper and Code)

### L-I1: Elasticity Tensor Inconsistency (CRITICAL)

**Severity:** Critical — affects all configurational force computations.

**Description:** Two different plane-strain elasticity tensors are used:
- `material.hpp` (stiffness assembly): Standard plane strain, **correct**
- `config_force.cpp` (Eshelby tensor): Non-standard parameterization using `E' = E/(1-ν²)`, **~1.6% error for ν=0.20**

**Effect:** Configurational forces G_tip are computed with a slightly incorrect material model. For ν = 0.20, C_11 is underestimated by 1.56%. The propagation criterion and direction are both affected.

**Fix required:** Replace the elasticity tensor in `config_force.cpp` with the standard formula:
```
λ = Eν / ((1+ν)(1-2ν)),  μ = E / (2(1+ν))
C_11 = λ + 2μ
```

### L-I2: Backtracking Uses Linear Energy Approximation (CRITICAL)

**Severity:** Critical — undermines the energy monotonicity proof.

**Description:** Paper Algorithm 2 requires re-solving the FEM system (Phases 1–5) after each crack extension trial to evaluate the true `E_h^{k+1}`. The implementation in `propagate_with_backtracking()` uses:
```cpp
E_trial.elastic -= force_magnitude * delta_a;   // ΔE ≈ -|G|·Δa  (linear approx)
E_trial.surface += Gc * delta_a;
```

**Effect:** The energy check is based on a first-order Taylor approximation, not the actual post-remesh energy. Proposition 3.1 (energy monotonicity) cannot be rigorously verified from the provided code as written.

**Fix required:** The backtracking loop must perform actual FEM re-assembly after each trial propagation.

### L-I3: Crack Normal Hardcoded (Major)

**Severity:** Major — affects correctness of G_tip for curved cracks.

**Description:** In `main.cpp`, the tip normal `ν̂` in `G_tip = G_mech - G_c·ν̂` is hardcoded as `{0, 1}` (upward) for all steps and both benchmarks. For TPB where the crack curves, this is incorrect once the crack deviates from horizontal.

**Fix required:** Compute ν̂ from the tangent of the last crack segment: if the last segment direction is `(dx, dy)`, then `ν̂ = (-dy, dx)/|...|` (perpendicular, pointing into uncracked region).

### L-I4: Propagation Threshold Missing θ_h (Minor)

**Description:** Paper criterion: `|G| ≥ G_c - θ_h` with `θ_h = 0.01 G_c`. Code uses: `|G| ≥ G_c` (no tolerance). This makes propagation slightly more difficult to trigger and may affect initiation load.

### L-I5: Direction Adjustment is Heuristic (Major)

**Description:** In `create_tpb_config()`:
```cpp
config.adjust_propagation_direction = [](double& dir_x, double& dir_y) {
    if (std::abs(dir_x) > std::abs(dir_y)) dir_x *= 0.5;  // HEURISTIC
    if (dir_y < 0.0) dir_y = -dir_y;                        // HEURISTIC
};
```
And in `create_sent_config()`:
```cpp
dir_y *= 0.1;   // damp vertical component for stability — HEURISTIC
```

**Effect:** These adjustments modify the crack propagation direction without physical justification. They are not described in the paper. This artificially biases crack paths and undermines the claim that "Crack direction aligns with the configurational force vector."

**Fix required:** Remove direction adjustment entirely; use G_tip direction directly.

### L-I6: compute_external_work Returns Zero (Major)

**Description:** `Assembler::compute_external_work()` always returns `0.0`. This means the total energy in the backtracking controller is `E = E_elastic + E_surface` (external work not subtracted). The paper's energy functional includes `-W_ext`.

**Effect:** Energy checks in backtracking are not computing the correct total energy.

### L-I7: Missing Core FEM Implementation (Critical for Reproducibility)

**Description:** Key methods declared in `assembler.hpp` but not provided in the zip:
- `Assembler::assemble()` — global stiffness assembly
- `Assembler::solve_linear()` — linear system solver  
- `Assembler::element_stiffness_CST()` — element stiffness matrix
- `Assembler::compute_tip_force()` — (in Assembler class, distinct from ConfigurationalForce class)

The repository as submitted is **not self-contained** and cannot be compiled or reproduced.

### L-I8: Double-Notch Benchmark Not Implemented (Major)

**Description:** The paper describes three benchmarks (SENT, TPB, Double-Notch). Only SENT and TPB are implemented. The Double-Notch case config is absent from `get_case_registry()`.

### L-I9: Phase-Field Comparison Not Included (Major)

**Description:** The speedup claim (6.1–12.0×) requires a reference AT2 phase-field implementation. No such implementation is provided. The speedup cannot be independently reproduced.

### L-I10: OpenMP Parallelization Not Implemented (Minor)

**Description:** Paper Section 3.2.2 claims Phase 2 (configurational force evaluation) uses "OpenMP-parallel loops." No OpenMP directives appear in the provided source code.

---

## III. Numerical Stability Risks

### L-N1: Mesh Quality Post-Remesh Not Verified

After `remesh_local()`, the code does not check that `θ_min ≥ 20°` is maintained. Poor-quality elements near the crack tip can degrade configurational force accuracy and cause ill-conditioning of K.

**Risk level:** Medium. The remesher is stated to enforce quality constraints, but no assertion in the main loop verifies this.

### L-N2: Condition Number of K Not Monitored

For highly refined meshes or near-degenerate geometries, the stiffness matrix K can become ill-conditioned. No condition number check or warning is implemented. GMRES with BoomerAMG may silently converge to an inaccurate solution.

### L-N3: Field Transfer Accuracy Not Verified

After remeshing, `FieldTransfer::transfer_displacement()` projects u_old → u_new. The projection error `‖u_proj - u_old‖` is not logged or bounded. Poor transfer quality can corrupt the warm-start and cause spurious energy changes in the subsequent step.

### L-N4: Angle Sanity Check is Heuristic

```cpp
if (delta_angle > 75.0) { propagate = false; }  // HEURISTIC
```
The 75° threshold has no theoretical justification. Legitimate branching or sharp turning could be incorrectly suppressed.

---

## IV. What Reviewers May Challenge

| Challenge | Risk Level | Response |
|---|---|---|
| "Energy monotonicity proof relies on full FEM re-solve inside backtracking, not a linear approximation" | HIGH | Acknowledge in paper; fix implementation |
| "The elasticity tensor in the Eshelby computation differs from the stiffness assembly by 1.6%" | HIGH | Fix `config_force.cpp`; run T2 test |
| "The crack normal is hardcoded — results for curved crack paths may be incorrect" | HIGH | Fix and re-run TPB |
| "Direction adjustment artificially biases results" | HIGH | Remove heuristics; quantify impact |
| "No phase-field reference → speedup cannot be verified" | MEDIUM | Include reference implementation or provide detailed comparison protocol |
| "Γ-convergence under remeshing is unresolved but presented ambiguously" | MEDIUM | Already acknowledged; ensure abstract/conclusions are unambiguous |
| "Double-Notch benchmark is missing from code" | MEDIUM | Implement or remove from paper |
| "compute_external_work = 0 means energy balance is incorrect" | HIGH | Fix before claiming thermodynamic consistency |
