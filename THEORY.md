# THEORY.md — Mathematical Summary and Equation-to-Code Mapping

## 1. Continuum Energy Functional

**Paper Eq. (1) / (2):**
```
E(u, Γ) = ∫_{Ω\Γ} W(∇u) dx  +  G_c · H¹(Γ)  -  W_ext(u)
```

**Code location:** `src/main.cpp` lines ~495–510  
```cpp
double elastic_energy  = 0.5 * uᵀKu              // ∫W dx ≈ ½uᵀKu
double surface_energy  = config.Gc * crack.total_length()  // G_c · H¹(Γ_h)
double total_energy    = elastic_energy + surface_energy - accumulated_external_work
```
**Note:** External work is accumulated via trapezoidal rule on reaction forces, not via boundary integral.

---

## 2. Plane-Strain Elasticity Tensor

**Paper Section 2 (plane strain, ε_zz = 0):**
```
C_11 = C_22 = E(1-ν) / ((1+ν)(1-2ν))
C_12 = C_21 = Eν / ((1+ν)(1-2ν))
C_33 = E / (2(1+ν))   [shear modulus G]
```

**Code — `material.hpp` (used in stiffness assembly):**
```cpp
double factor = E / ((1.0 + nu) * (1.0 - 2.0*nu));
D[0] = factor * (1.0 - nu);   // C_11  ✓
D[8] = factor * (1.0-2.0*nu)/2.0;  // C_33 ✓
```

**Code — `config_force.cpp` (used in Eshelby tensor):**
```cpp
double E_prime = E_ / (1.0 - nu_*nu_);          // E' = E/(1-ν²)
double lambda_prime = E_prime * nu_ / (1.0-nu_); // λ' = E'ν/(1-ν)
double mu = E_ / (2.0*(1.0+nu_));
C(0,0) = lambda_prime + 2.0*mu;
```

### ⚠️ DISCREPANCY (Critical)

For ν = 0.20:
- Correct C_11 = 30e9 × 0.8 / (1.2 × 0.6) = **33.33 GPa**
- `config_force.cpp` C_11 = 7.81 + 25.0 = **32.81 GPa** (error ≈ 1.6%)

The formula in `config_force.cpp` is non-standard. It appears to use plane-stress notation (`E'=E/(1-ν²)`) but applies it to a plane-strain context. This creates a **1.6% error** in the Eshelby tensor for ν=0.20 and **grows with ν**.

**Required fix:**
```cpp
// config_force.cpp — replace with standard plane strain:
double lambda = E_ * nu_ / ((1.0+nu_) * (1.0-2.0*nu_));
double mu     = E_ / (2.0*(1.0+nu_));
C(0,0) = lambda + 2.0*mu;
C(1,1) = lambda + 2.0*mu;
C(0,1) = C(1,0) = lambda;
C(2,2) = mu;
```

---

## 3. Configurational (Eshelby) Force

**Paper Eq. B.8:**
```
Σ = W(∇u)·I  −  (∇u)ᵀ · σ     [Eshelby stress tensor, 2×2]
```

**Paper Eq. B.14:**
```
G_I = -∑_{K ∈ Star(I)} ∫_K Σ : ∇N_I dx
```

**Paper Eq. B.18:**
```
G_tip = G_mech − G_c · ν̂     [ν̂ = unit outward normal at tip]
```

**Code — `config_force.cpp`:**

| Equation | Function | File | Notes |
|---|---|---|---|
| B.8 (Eshelby) | `compute_eshelby_tensor()` | `config_force.cpp` | ✓ Correct formula |
| B.14 (J-integral) | `integrate_element_contribution()` | `config_force.cpp` | ✓ Sign correct |
| B.18 (G_tip) | `compute_tip_force()` end | `config_force.cpp` | ✓ Implements `G−G_c·ν` |
| B.20–B.21 (Duffy) | `duffy_integrate_tip_element()` | `config_force.cpp` | ⚠ See note below |

### ⚠️ Duffy Architecture Issue

For P1 (linear) elements, `∇u` is constant per element — there is **no integrable singularity** in the Eshelby tensor at the discrete level. The Duffy transformation is designed to cancel a `r^{-1/2}` singularity in the displacement gradient, which is not present in P1 FEM.

**Consequence:** The Duffy integration is mathematically equivalent to standard Gauss quadrature for P1 elements but adds ~7× computational overhead. The justification in Appendix B.3.1 is physically motivated (continuous fields have singularity) but the discrete implementation does not actually encounter this singularity.

**This does not cause incorrect results** but the claim of singularity treatment via Duffy for P1 elements is misleading and should be clarified.

---

## 4. Crack Tip Normal (ν̂)

**Paper:** ν̂ points **into the uncracked region**, defined as the unit normal to the crack surface at the tip.

**Code — `main.cpp` line ~534:**
```cpp
std::array<double, 2> crack_normal = {0.0, 1.0};  // HARDCODED upward
```

### ⚠️ DISCREPANCY (Major)

The crack normal is hardcoded as `{0, 1}` (pointing upward) for ALL cases including SENT. For SENT (horizontal crack propagating rightward), the correct ν̂ should be the normal to the crack surface, i.e., `{0, ±1}` — which happens to be correct by coincidence for a horizontal crack. However, for TPB where the crack curves, this hardcoded normal will produce errors in G_tip as the crack angle changes.

**Required fix:** Compute ν̂ dynamically from the last crack segment direction.

---

## 5. Propagation Criterion

**Paper (discrete Griffith criterion):**
```
‖G_{i*}‖ ≥ G_c − θ_h,   θ_h = 0.01 G_c
```

**Code — `main.cpp` line ~557:**
```cpp
bool propagate = (G.magnitude >= config.Gc) && std::isfinite(G.magnitude);
```

### ⚠️ DISCREPANCY (Minor)

Code uses `G.magnitude ≥ G_c` (exact equality threshold). Paper uses `G_c − θ_h` with `θ_h = 0.01 G_c`. The tolerance `θ_h` is absent from the propagation check in `main.cpp`.

---

## 6. Backtracking (Energy Monotonicity)

**Paper Algorithm 2, Phase 6:**
```
Accept if: E_h^{k+1} ≤ E_h^k + ε_tol · |E_h^k|,  ε_tol = 10⁻⁶
Backtrack: Δa ← Δa/2
Terminate backtrack if: Δa < h_min
```

**Code — `adaptive_stepping.hpp`:**
```cpp
double tolerance = params_.epsilon_tol * std::abs(E_old.total);  // ε_tol × |E|  ✓
bool energy_acceptable = (delta_E <= tolerance);                  // ✓
delta_a *= params_.reduction_factor;  // reduction_factor = 0.5   ✓
```

**Code — `assembler.cpp` (`propagate_with_backtracking`):**
```cpp
E_trial.elastic -= force_magnitude * delta_a;  // Linear approximation!
E_trial.surface += Gc * delta_a;
```

### ⚠️ DISCREPANCY (Critical)

**Paper Algorithm 2** re-solves the FEM system after each crack extension trial (Phase 1 → Phase 5 → Phase 6 in the loop). The code uses a **first-order Taylor approximation**: `ΔE ≈ (G_c - |G|)·Δa`. This is only valid for infinitesimal steps and does not correctly evaluate `E_h^{k+1}` after remeshing and field projection. The energy monotonicity guarantee (Proposition 3.1) requires the actual post-remesh energy.

---

## 7. Stiffness Assembly (CST Elements)

**Paper Section 3.2 (P1 elements, constant strain triangles):**

The element stiffness matrix is `K_e = t · B^T C B · A_e` where:
- `B` is the 3×6 strain-displacement matrix (constant for CST)
- `C` is the elasticity tensor (3×3 Voigt)
- `A_e` is element area
- `t` is thickness

**Code:** `Assembler::element_stiffness_CST()` is declared in `assembler.hpp` but **its full implementation is not present** in the provided source files. Only the new methods (`propagate_with_backtracking`, energy functions) appear in `assembler.cpp`.

### ⚠️ INCOMPLETE SOURCE

Core FEM routines (`assemble`, `solve_linear`, `compute_tip_force` within Assembler, `element_stiffness_CST`) are declared but not provided. The repository as submitted is **not self-contained**.

---

## 8. Equation → Function Mapping Table

| Paper Equation | Description | File | Function/Line |
|---|---|---|---|
| Eq. (1) | Total energy functional | `main.cpp` | ~495–510 |
| Eq. (2) | Discrete energy | `main.cpp` | ~495–510 |
| Eq. B.8 | Eshelby tensor Σ | `config_force.cpp` | `compute_eshelby_tensor()` |
| Eq. B.14 | Nodal config. force | `config_force.cpp` | `integrate_element_contribution()` |
| Eq. B.18 | Tip force G_tip | `config_force.cpp` | `compute_tip_force()` |
| Eq. B.20–B.21 | Duffy transformation | `config_force.cpp` | `duffy_integrate_tip_element()` |
| Alg. 2, Phase 1 | Newton–Raphson solve | `assembler.hpp` | `solve_linear()` [MISSING IMPL] |
| Alg. 2, Phase 3 | Crack extension Δa=βh | `assembler.cpp` | `propagate_with_backtracking()` |
| Alg. 2, Phase 4 | CDT remeshing | `remesher.cpp` | `remesh_local()` |
| Alg. 2, Phase 5 | L²-projection | `field_transfer.cpp` | `transfer_displacement()` |
| Alg. 2, Phase 6 | Energy check | `adaptive_stepping.hpp` | `attempt_propagation()` |
| Prop. 3.1 | Energy monotonicity | `adaptive_stepping.hpp` | `attempt_propagation()` |
| Table 4 | Material params (TPB) | `main.cpp` | `create_tpb_config()` |
| Table 4 | Material params (SENT) | `main.cpp` | `create_sent_config()` |
| App. C | Load stepping | `adaptive_stepping.hpp` | `AdaptiveStepController` |
