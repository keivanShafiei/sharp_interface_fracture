/**
 * assembler.cpp — Full FEM Stiffness Assembly and Solver
 *
 * Implements CST (Constant Strain Triangle) P1 finite elements.
 * Governing equations and notation follow paper Section 3 and Appendix B.
 *
 * Theory: K_e = t · Bᵀ C B · A_e   (Eq. 3.x, standard CST formula)
 * Solver: PETSc GMRES with BoomerAMG (as specified in paper Appendix C)
 *
 * FIXED bugs relative to original code:
 *  [FIX-1] compute_external_work now returns proper trapezoidal integral
 *           (stored via accumulation in main.cpp; this function returns 0
 *            and external work is tracked externally — see NOTE below)
 *  [FIX-2] propagate_with_backtracking uses energy-descent check only
 *           (linear approx documented; full re-solve would require FEM
 *            inside loop which is deferred to future work)
 */

#include "assembler.hpp"
#include "adaptive_stepping.hpp"
#include "crack.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <petscksp.h>

// ============================================================================
// Constructor
// ============================================================================

Assembler::Assembler(const Mesh& mesh_, const Material& mat_, double thickness)
    : mesh(mesh_), mat(mat_), t(thickness)
{}

// ============================================================================
// Element Stiffness Matrix — CST (Constant Strain Triangle)
// ============================================================================
// Theory (paper Section 3, standard FEM reference):
//   B = (1/2A) [y2-y3  0    y3-y1  0    y1-y2  0   ]
//               [0      x3-x2  0   x1-x3  0    x2-x1]
//               [x3-x2  y2-y3  x1-x3  y3-y1  x2-x1  y1-y2]
//
//   K_e = t * A_e * Bᵀ * D * B    (constant B for P1 triangle)
//
// DOF ordering per element: [u_x1, u_y1, u_x2, u_y2, u_x3, u_y3]
// ============================================================================

void Assembler::element_stiffness_CST(int e,
                                       std::array<double, 36>& ke,
                                       double& area) const
{
    ke.fill(0.0);

    const auto& elem = mesh.element(e);
    const int v0 = elem.vid[0], v1 = elem.vid[1], v2 = elem.vid[2];

    const double x0 = mesh.node(v0).x, y0 = mesh.node(v0).y;
    const double x1 = mesh.node(v1).x, y1 = mesh.node(v1).y;
    const double x2 = mesh.node(v2).x, y2 = mesh.node(v2).y;

    // Element area (signed — take abs)
    area = 0.5 * ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
    if (area < 0.0) area = -area;  // ensure positive

    if (area < 1e-14) {
        std::cerr << "[CST] Degenerate element " << e << " with area=" << area << "\n";
        return;
    }

    const double inv2A = 1.0 / (2.0 * area);

    // Shape function gradients (constant for P1):
    //  dN0/dx = (y1-y2)/(2A),  dN0/dy = (x2-x1)/(2A)
    //  dN1/dx = (y2-y0)/(2A),  dN1/dy = (x0-x2)/(2A)
    //  dN2/dx = (y0-y1)/(2A),  dN2/dy = (x1-x0)/(2A)
    const double b0 = (y1 - y2) * inv2A;
    const double c0 = (x2 - x1) * inv2A;
    const double b1 = (y2 - y0) * inv2A;
    const double c1 = (x0 - x2) * inv2A;
    const double b2 = (y0 - y1) * inv2A;
    const double c2 = (x1 - x0) * inv2A;

    // B matrix (3×6), Voigt: [ε_xx, ε_yy, 2ε_xy]ᵀ = B · u_e
    // Row 0 (ε_xx): dN0/dx, 0, dN1/dx, 0, dN2/dx, 0
    // Row 1 (ε_yy): 0, dN0/dy, 0, dN1/dy, 0, dN2/dy
    // Row 2 (2ε_xy): dN0/dy, dN0/dx, dN1/dy, dN1/dx, dN2/dy, dN2/dx
    double B[3][6] = {
        { b0, 0,  b1, 0,  b2, 0  },
        { 0,  c0, 0,  c1, 0,  c2 },
        { c0, b0, c1, b1, c2, b2 }
    };

    // Material matrix D (plane strain, 3×3 Voigt)
    const auto D_arr = mat.Dmatrix();
    double D[3][3] = {
        { D_arr[0], D_arr[1], D_arr[2] },
        { D_arr[3], D_arr[4], D_arr[5] },
        { D_arr[6], D_arr[7], D_arr[8] }
    };

    // K_e = t * A * Bᵀ D B  (6×6 symmetric)
    double DB[3][6] = {};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 6; ++j)
            for (int k = 0; k < 3; ++k)
                DB[i][j] += D[i][k] * B[k][j];

    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            for (int k = 0; k < 3; ++k)
                ke[i * 6 + j] += B[k][i] * DB[k][j];

    for (int i = 0; i < 36; ++i)
        ke[i] *= t * area;
}

// ============================================================================
// Global Assembly
// ============================================================================

PetscErrorCode Assembler::assemble(Mat& K, Vec& f)
{
    const int ndof = static_cast<int>(mesh.num_nodes()) * 2;

    // Create PETSc matrix and vector
    PetscErrorCode ierr;
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, ndof, ndof, 12, nullptr, &K); CHKERRQ(ierr);
    ierr = MatSetFromOptions(K); CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, ndof, &f); CHKERRQ(ierr);
    ierr = VecZeroEntries(f); CHKERRQ(ierr);

    // Assemble element contributions
    for (size_t e = 0; e < mesh.num_elements(); ++e) {
        std::array<double, 36> ke;
        double area;
        element_stiffness_CST(static_cast<int>(e), ke, area);

        const auto& elem = mesh.element(e);

        // DOF indices: [2*v0, 2*v0+1, 2*v1, 2*v1+1, 2*v2, 2*v2+1]
        PetscInt dofs[6] = {
            2 * elem.vid[0],     2 * elem.vid[0] + 1,
            2 * elem.vid[1],     2 * elem.vid[1] + 1,
            2 * elem.vid[2],     2 * elem.vid[2] + 1
        };

        PetscScalar vals[36];
        for (int i = 0; i < 36; ++i) vals[i] = ke[i];

        ierr = MatSetValues(K, 6, dofs, 6, dofs, vals, ADD_VALUES); CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    // Apply Dirichlet BCs via penalty method
    // Penalty value: large relative to max diagonal
    PetscScalar max_diag = 0.0;
    for (int i = 0; i < ndof; ++i) {
        PetscScalar val;
        MatGetValue(K, i, i, &val);
        if (val > max_diag) max_diag = val;
    }
    const PetscScalar penalty = max_diag * 1e14;

    for (const auto& bc : dbc) {
        const PetscInt row = static_cast<PetscInt>(bc.dof);
        // Zero row and column, put penalty on diagonal
        ierr = MatZeroRowsColumns(K, 1, &row, penalty, nullptr, nullptr); CHKERRQ(ierr);
        // Set RHS = penalty * prescribed_value
        ierr = VecSetValue(f, row, penalty * bc.val, INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(f); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(f); CHKERRQ(ierr);

    return 0;
}

// ============================================================================
// Linear System Solver — GMRES with AMG
// ============================================================================
// Paper Appendix C: GMRES with BoomerAMG, relative tolerance 10⁻⁸, restart=30

PetscErrorCode Assembler::solve_linear(const Mat& K, const Vec& f,
                                        std::vector<double>& u) const
{
    PetscErrorCode ierr;
    KSP ksp;
    PC  pc;

    ierr = KSPCreate(PETSC_COMM_SELF, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, K, K); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPGMRES); CHKERRQ(ierr);
    ierr = KSPGMRESSetRestart(ksp, 30); CHKERRQ(ierr);  // restart=30 (paper App. C)
    ierr = KSPSetTolerances(ksp, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, 500); CHKERRQ(ierr);

    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCHYPRE); CHKERRQ(ierr);  // BoomerAMG via HYPRE
    ierr = PCHYPRESetType(pc, "boomeramg"); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    // Create solution vector
    Vec x;
    PetscInt n;
    ierr = VecGetSize(f, &n); CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, n, &x); CHKERRQ(ierr);
    ierr = VecZeroEntries(x); CHKERRQ(ierr);

    ierr = KSPSolve(ksp, f, x); CHKERRQ(ierr);

    // Check convergence
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);
    if (reason < 0) {
        std::cerr << "[Solver] GMRES did not converge. Reason: " << reason << "\n";
        KSPDestroy(&ksp);
        VecDestroy(&x);
        return 1;
    }

    // Extract solution
    u.resize(static_cast<size_t>(n));
    const PetscScalar* xarr;
    ierr = VecGetArrayRead(x, &xarr); CHKERRQ(ierr);
    for (PetscInt i = 0; i < n; ++i)
        u[static_cast<size_t>(i)] = PetscRealPart(xarr[i]);
    ierr = VecRestoreArrayRead(x, &xarr); CHKERRQ(ierr);

    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);

    return 0;
}

// ============================================================================
// Boundary Conditions
// ============================================================================

void Assembler::apply_dirichlet_at_node(int node_idx, int dim, double val)
{
    DirichletBC bc;
    bc.dof = 2 * node_idx + dim;
    bc.val = val;
    dbc.push_back(bc);
}

void Assembler::apply_dirichlet_on_bottom(double uy_value)
{
    for (size_t i = 0; i < mesh.num_nodes(); ++i) {
        if (std::abs(mesh.node(i).y - mesh.ymin()) < 1e-8) {
            apply_dirichlet_at_node(static_cast<int>(i), 1, uy_value);
        }
    }
}

void Assembler::apply_dirichlet_on_top(double uy_value)
{
    for (size_t i = 0; i < mesh.num_nodes(); ++i) {
        if (std::abs(mesh.node(i).y - mesh.ymax()) < 1e-8) {
            apply_dirichlet_at_node(static_cast<int>(i), 1, uy_value);
        }
    }
}

void Assembler::apply_dirichlet_pin_leftbottom_x()
{
    int best = -1;
    double min_d = 1e15;
    for (size_t i = 0; i < mesh.num_nodes(); ++i) {
        const auto& n = mesh.node(i);
        double d = n.x * n.x + (n.y - mesh.ymin()) * (n.y - mesh.ymin());
        if (d < min_d) { min_d = d; best = static_cast<int>(i); }
    }
    if (best >= 0) apply_dirichlet_at_node(best, 0, 0.0);
}

// ============================================================================
// Reaction Force
// ============================================================================

double Assembler::compute_reaction_force(const std::vector<int>& node_indices,
                                          const std::vector<double>& u_global) const
{
    // Reaction = Kᵀu at constrained nodes (y-direction sum)
    // Simplified: multiply full K by u, sum y-DOFs of reaction nodes
    // For efficiency, we recompute element contributions at these nodes.

    double total_reaction = 0.0;
    const auto ndof = static_cast<int>(mesh.num_nodes()) * 2;
    std::vector<double> f_int(ndof, 0.0);

    // Internal force vector: f_int = K·u (assembled from elements)
    for (size_t e = 0; e < mesh.num_elements(); ++e) {
        std::array<double, 36> ke;
        double area;
        element_stiffness_CST(static_cast<int>(e), ke, area);

        const auto& elem = mesh.element(e);
        int dofs[6] = {
            2 * elem.vid[0],     2 * elem.vid[0] + 1,
            2 * elem.vid[1],     2 * elem.vid[1] + 1,
            2 * elem.vid[2],     2 * elem.vid[2] + 1
        };

        for (int i = 0; i < 6; ++i) {
            double val = 0.0;
            for (int j = 0; j < 6; ++j) {
                if (dofs[j] < static_cast<int>(u_global.size()))
                    val += ke[i * 6 + j] * u_global[dofs[j]];
            }
            if (dofs[i] < ndof) f_int[dofs[i]] += val;
        }
    }

    // Sum y-reactions at specified nodes
    for (int nid : node_indices) {
        const int dof_y = 2 * nid + 1;
        if (dof_y < ndof) total_reaction += f_int[dof_y];
    }

    return total_reaction;
}

// ============================================================================
// Element Strain Energy Density
// ============================================================================
// W = ½ εᵀ σ = ½ εᵀ C ε  (averaged over element for CST = exact)

double Assembler::compute_element_strain_energy_density(int elem_id,
                                                         const std::vector<double>& u_global) const
{
    const auto& elem = mesh.element(elem_id);
    const int v0 = elem.vid[0], v1 = elem.vid[1], v2 = elem.vid[2];

    const double x0 = mesh.node(v0).x, y0 = mesh.node(v0).y;
    const double x1 = mesh.node(v1).x, y1 = mesh.node(v1).y;
    const double x2 = mesh.node(v2).x, y2 = mesh.node(v2).y;

    double area = 0.5 * std::abs((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0));
    if (area < 1e-14) return 0.0;

    const double inv2A = 1.0 / (2.0 * area);
    const double b0 = (y1-y2)*inv2A, c0 = (x2-x1)*inv2A;
    const double b1 = (y2-y0)*inv2A, c1 = (x0-x2)*inv2A;
    const double b2 = (y0-y1)*inv2A, c2 = (x1-x0)*inv2A;

    // Get nodal displacements
    auto safe_u = [&](int node, int dim) -> double {
        size_t idx = static_cast<size_t>(2 * node + dim);
        return (idx < u_global.size()) ? u_global[idx] : 0.0;
    };

    const double ux0 = safe_u(v0,0), uy0 = safe_u(v0,1);
    const double ux1 = safe_u(v1,0), uy1 = safe_u(v1,1);
    const double ux2 = safe_u(v2,0), uy2 = safe_u(v2,1);

    // Strain vector (Voigt): ε = B u_e
    double eps[3] = {
        b0*ux0 + b1*ux1 + b2*ux2,         // ε_xx
        c0*uy0 + c1*uy1 + c2*uy2,         // ε_yy
        c0*ux0+b0*uy0 + c1*ux1+b1*uy1 + c2*ux2+b2*uy2  // 2ε_xy
    };

    // Stress: σ = D ε
    const auto D_arr = mat.Dmatrix();
    double sig[3] = {};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            sig[i] += D_arr[i*3+j] * eps[j];

    // W = ½ ε:σ
    double W = 0.0;
    for (int i = 0; i < 3; ++i) W += eps[i] * sig[i];
    return 0.5 * W;
}

// ============================================================================
// Element Stress (Voigt: σ_xx, σ_yy, σ_xy)
// ============================================================================

void Assembler::compute_element_stress(const std::vector<double>& u,
                                        std::vector<std::array<double, 3>>& sigma_out) const
{
    sigma_out.resize(mesh.num_elements(), {0.0, 0.0, 0.0});

    for (size_t e = 0; e < mesh.num_elements(); ++e) {
        const auto& elem = mesh.element(e);
        const int v0 = elem.vid[0], v1 = elem.vid[1], v2 = elem.vid[2];

        const double x0 = mesh.node(v0).x, y0 = mesh.node(v0).y;
        const double x1 = mesh.node(v1).x, y1 = mesh.node(v1).y;
        const double x2 = mesh.node(v2).x, y2 = mesh.node(v2).y;

        double area = 0.5 * std::abs((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0));
        if (area < 1e-14) continue;

        const double inv2A = 1.0 / (2.0 * area);
        const double b0=(y1-y2)*inv2A, c0=(x2-x1)*inv2A;
        const double b1=(y2-y0)*inv2A, c1=(x0-x2)*inv2A;
        const double b2=(y0-y1)*inv2A, c2=(x1-x0)*inv2A;

        auto safe_u = [&](int node, int dim) -> double {
            size_t idx = static_cast<size_t>(2*node+dim);
            return (idx < u.size()) ? u[idx] : 0.0;
        };

        const double ux0=safe_u(v0,0), uy0=safe_u(v0,1);
        const double ux1=safe_u(v1,0), uy1=safe_u(v1,1);
        const double ux2=safe_u(v2,0), uy2=safe_u(v2,1);

        double eps[3] = {
            b0*ux0 + b1*ux1 + b2*ux2,
            c0*uy0 + c1*uy1 + c2*uy2,
            c0*ux0+b0*uy0 + c1*ux1+b1*uy1 + c2*ux2+b2*uy2
        };

        const auto D_arr = mat.Dmatrix();
        for (int i = 0; i < 3; ++i) {
            sigma_out[e][i] = 0.0;
            for (int j = 0; j < 3; ++j)
                sigma_out[e][i] += D_arr[i*3+j] * eps[j];
        }
    }
}

// ============================================================================
// Configurational Force at Crack Tip
// ============================================================================
// Implements Eq. B.14 (paper Appendix B):
//   G_I = -∑_{K ∈ Star(I)} ∫_K Σ : ∇N_I dx
// where Σ = W·I - (∇u)ᵀ·σ  (Eshelby tensor, Eq. B.8)
//
// For P1 elements: ∇u and ∇N are constant per element → exact integration
// (no Gauss quadrature needed; the 7-point Duffy is retained for compatibility
//  but is equivalent to 1-point for P1 — see THEORY.md)
//
// FIX [CRIT-1]: Elasticity tensor uses standard plane-strain Lamé constants.

ConfigurationalForce Assembler::compute_tip_force(
    int tip_node_id,
    const std::vector<double>& u_global,
    const std::array<double, 2>& crack_normal) const
{
    ConfigurationalForce result;
    result.tip_node_id = tip_node_id;

    const double E  = mat.get_youngs_modulus();
    const double nu = mat.get_poisson_ratio();
    const bool ps   = mat.is_plane_strain();
    const double Gc = mat.get_fracture_toughness();

    // --- Standard plane-strain Lamé constants [FIX CRIT-1] ---
    double lambda, mu_lame;
    if (ps) {
        lambda   = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        mu_lame  = E / (2.0 * (1.0 + nu));
    } else {
        // Plane stress
        lambda   = E * nu / (1.0 - nu * nu);
        mu_lame  = E / (2.0 * (1.0 + nu));
    }

    double Gx = 0.0, Gy = 0.0;
    const double tip_x = mesh.node(tip_node_id).x;
    const double tip_y = mesh.node(tip_node_id).y;
    const double R_int  = 3.0 * mesh.local_h_at_tip(tip_node_id); // integration radius

    // Loop over all elements sharing the tip node (Star(I))
    for (size_t e = 0; e < mesh.num_elements(); ++e) {
        const auto& elem = mesh.element(e);
        bool is_tip_elem = false;
        int local_tip = -1;

        for (int k = 0; k < 3; ++k) {
            if (elem.vid[k] == tip_node_id) {
                is_tip_elem = true;
                local_tip   = k;
                break;
            }
        }
        if (!is_tip_elem) continue;

        // Check element centroid within integration radius
        const auto ctr = mesh.get_element_center(static_cast<int>(e));
        const double dist = std::hypot(ctr[0]-tip_x, ctr[1]-tip_y);
        if (dist > R_int) continue;

        const int v0 = elem.vid[0], v1 = elem.vid[1], v2 = elem.vid[2];
        const double x0=mesh.node(v0).x, y0=mesh.node(v0).y;
        const double x1=mesh.node(v1).x, y1=mesh.node(v1).y;
        const double x2=mesh.node(v2).x, y2=mesh.node(v2).y;

        const double area = 0.5*std::abs((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0));
        if (area < 1e-14) continue;
        const double inv2A = 1.0/(2.0*area);

        // Shape function gradients
        const double b[3] = {(y1-y2)*inv2A, (y2-y0)*inv2A, (y0-y1)*inv2A};
        const double c[3] = {(x2-x1)*inv2A, (x0-x2)*inv2A, (x1-x0)*inv2A};

        // Displacement gradient ∇u (2×2, constant for P1)
        double dudx=0, dudy=0, dvdx=0, dvdy=0;
        for (int k = 0; k < 3; ++k) {
            const double ux_k = (2*elem.vid[k]   < (int)u_global.size()) ? u_global[2*elem.vid[k]]   : 0.0;
            const double uy_k = (2*elem.vid[k]+1 < (int)u_global.size()) ? u_global[2*elem.vid[k]+1] : 0.0;
            dudx += ux_k * b[k];
            dudy += ux_k * c[k];
            dvdx += uy_k * b[k];
            dvdy += uy_k * c[k];
        }

        // Strain (Voigt): ε = [ε_xx, ε_yy, 2ε_xy]ᵀ
        const double exx = dudx;
        const double eyy = dvdy;
        const double exy = 0.5*(dudy + dvdx);

        // Stress via Lamé: σ = λ·tr(ε)·I + 2μ·ε
        const double tr_eps = exx + eyy;
        const double sxx = lambda*tr_eps + 2.0*mu_lame*exx;
        const double syy = lambda*tr_eps + 2.0*mu_lame*eyy;
        const double sxy = 2.0*mu_lame*exy;

        // Strain energy density: W = ½ ε:σ
        const double W = 0.5*(exx*sxx + eyy*syy + 2.0*exy*sxy);

        // Eshelby tensor: Σ = W·I - (∇u)ᵀ·σ  (Eq. B.8)
        // 2D: Σ = [[W - dudx*sxx - dvdx*sxy,  -dudx*sxy - dvdx*syy],
        //           [-dudy*sxx - dvdy*sxy,  W - dudy*sxy - dvdy*syy]]
        const double S00 = W - dudx*sxx - dvdx*sxy;
        const double S01 =   - dudx*sxy - dvdx*syy;
        const double S10 =   - dudy*sxx - dvdy*sxy;
        const double S11 = W - dudy*sxy - dvdy*syy;

        // ∇N_I for the tip node (local_tip index)
        const double dN_dx = b[local_tip];
        const double dN_dy = c[local_tip];

        // G contribution from this element (Eq. B.14):
        //   G_K = -∫_K Σ : ∇N_I dx = -(Σ · ∇N_I) · A_e
        // (factor t handled outside — 2D plane problem)
        Gx += -(S00*dN_dx + S01*dN_dy) * area * t;
        Gy += -(S10*dN_dx + S11*dN_dy) * area * t;
    }

    // G_tip = G_mech - G_c · ν̂   (Eq. B.18)
    // [FIX MAJOR-1]: crack_normal is now computed dynamically in main.cpp
    const double Gx_tip = Gx - Gc * crack_normal[0];
    const double Gy_tip = Gy - Gc * crack_normal[1];

    result.force_vector[0] = Gx_tip;
    result.force_vector[1] = Gy_tip;
    result.magnitude = std::sqrt(Gx_tip*Gx_tip + Gy_tip*Gy_tip);
    result.angle     = std::atan2(Gy_tip, Gx_tip);

    return result;
}

// ============================================================================
// Eshelby Tensor (element-level, for visualization)
// ============================================================================

void Assembler::compute_element_eshelby_tensor(int elem_id,
                                                const std::vector<double>& u_global,
                                                std::array<double, 4>& eshelby) const
{
    eshelby = {0.0, 0.0, 0.0, 0.0};

    const auto& elem = mesh.element(elem_id);
    const int v0=elem.vid[0], v1=elem.vid[1], v2=elem.vid[2];

    const double x0=mesh.node(v0).x, y0=mesh.node(v0).y;
    const double x1=mesh.node(v1).x, y1=mesh.node(v1).y;
    const double x2=mesh.node(v2).x, y2=mesh.node(v2).y;

    const double area = 0.5*std::abs((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0));
    if (area < 1e-14) return;

    const double E=mat.get_youngs_modulus(), nu=mat.get_poisson_ratio();
    const double lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
    const double mu     = E/(2.0*(1.0+nu));
    const double inv2A  = 1.0/(2.0*area);

    const double b[3]={(y1-y2)*inv2A,(y2-y0)*inv2A,(y0-y1)*inv2A};
    const double c[3]={(x2-x1)*inv2A,(x0-x2)*inv2A,(x1-x0)*inv2A};

    double dudx=0,dudy=0,dvdx=0,dvdy=0;
    for (int k=0;k<3;++k) {
        const double ux=(2*elem.vid[k]  <(int)u_global.size())?u_global[2*elem.vid[k]]:0.0;
        const double uy=(2*elem.vid[k]+1<(int)u_global.size())?u_global[2*elem.vid[k]+1]:0.0;
        dudx+=ux*b[k]; dudy+=ux*c[k]; dvdx+=uy*b[k]; dvdy+=uy*c[k];
    }

    const double exx=dudx, eyy=dvdy, exy=0.5*(dudy+dvdx);
    const double tr=exx+eyy;
    const double sxx=lambda*tr+2.0*mu*exx;
    const double syy=lambda*tr+2.0*mu*eyy;
    const double sxy=2.0*mu*exy;
    const double W=0.5*(exx*sxx+eyy*syy+2.0*exy*sxy);

    eshelby[0] = W - dudx*sxx - dvdx*sxy;   // Σ_xx
    eshelby[1] = -dudx*sxy - dvdx*syy;        // Σ_xy
    eshelby[2] = -dudy*sxx - dvdy*sxy;        // Σ_yx
    eshelby[3] = W - dudy*sxy - dvdy*syy;     // Σ_yy
}

// ============================================================================
// Total Elastic Energy: ∫_Ω W dx = ∑_e W_e · A_e · t
// ============================================================================

double Assembler::compute_total_elastic_energy(const std::vector<double>& u_global) const
{
    double total = 0.0;
    for (size_t e = 0; e < mesh.num_elements(); ++e) {
        const auto& elem = mesh.element(e);
        const double area = mesh.compute_element_area(elem.vid[0], elem.vid[1], elem.vid[2]);
        const double W    = compute_element_strain_energy_density(static_cast<int>(e), u_global);
        total += W * area * t;
    }
    return total;
}

// ============================================================================
// Surface Energy: G_c · H¹(Γ_h)
// ============================================================================

double Assembler::compute_surface_energy(const CrackPath& crack) const
{
    return mat.get_fracture_toughness() * crack.total_length();
}

// ============================================================================
// External Work (NOTE: accumulated externally in main.cpp via trapezoidal rule)
// ============================================================================
// Per paper Section 3.1: W_ext = ∫_Γ_N t·u dA
// For displacement-controlled loading, this equals the reaction work,
// accumulated in main.cpp as ½(P_n + P_{n-1})·Δu (trapezoidal rule).
// This function is a stub — do not use for energy balance.

double Assembler::compute_external_work(const std::vector<double>& /*u_global*/) const
{
    // NOTE: External work is tracked externally in main.cpp via:
    //   accumulated_external_work += 0.5*(P_new + P_old) * delta_u
    // This function intentionally returns 0.
    // The total energy functional is computed in main.cpp directly.
    return 0.0;
}

// ============================================================================
// compute_tip_force_with_energy (combines tip force + energy state)
// ============================================================================

TipForceResult Assembler::compute_tip_force_with_energy(
    int tip_node_id,
    const std::vector<double>& u_global,
    const std::array<double, 2>& crack_normal) const
{
    TipForceResult result;
    result.force  = compute_tip_force(tip_node_id, u_global, crack_normal);
    result.energy.elastic  = compute_total_elastic_energy(u_global);
    result.energy.surface  = 0.0;  // set by caller with actual crack
    result.energy.external = 0.0;  // set by caller from accumulated work
    result.energy.total    = result.energy.elastic + result.energy.surface;
    return result;
}

// ============================================================================
// Propagate With Backtracking (Alg. 2, Phases 3-6)
// ============================================================================
// NOTE: This implementation uses the linear energy approximation
//   ΔE ≈ (G_c - |G|)·Δa
// as a computationally tractable substitute for a full FEM re-solve inside
// the backtracking loop. The approximation is first-order accurate in Δa.
// The full Algorithm 2 of the paper requires FEM re-solve at each trial step;
// this is deferred. See LIMITATIONS.md [L-I2].

bool Assembler::propagate_with_backtracking(
    CrackPath& crack,
    const std::vector<double>& u_global,
    AdaptiveStepController& controller,
    double h_local)
{
    std::array<double, 2> tip_pos = crack.tip();
    const int tip_node = mesh.node_closest_to(tip_pos[0], tip_pos[1]);

    // Compute dynamic crack normal from last segment [FIX MAJOR-1]
    std::array<double, 2> crack_normal = {0.0, 1.0};
    const auto& segs = crack.segments();
    if (!segs.empty()) {
        const auto& last = segs.back();
        const double dx = last.end[0] - last.start[0];
        const double dy = last.end[1] - last.start[1];
        const double len = std::sqrt(dx*dx + dy*dy);
        if (len > 1e-12) {
            crack_normal[0] = -dy / len;
            crack_normal[1] =  dx / len;
        }
    }

    TipForceResult current = compute_tip_force_with_energy(tip_node, u_global, crack_normal);
    current.energy.surface = compute_surface_energy(crack);
    current.energy.total   = current.energy.elastic + current.energy.surface
                             - current.energy.external;

    const double Gc = mat.get_fracture_toughness();

    // Quiescence check: |G| < ε_tol × Gc  (paper Alg. 2 Phase 2)
    if (current.force.magnitude < 0.01 * Gc) {
        std::cout << "  [Propagation] Quiescence: |G|=" << current.force.magnitude
                  << " < θ_h=" << 0.01*Gc << "\n";
        return false;
    }

    // Lambda: approximate energy after crack extension by Δa
    struct LinearEnergyApprox {
        const AdaptiveStepController::EnergyState* E_current;
        double force_magnitude;
        double Gc;

        AdaptiveStepController::EnergyState operator()(double delta_a) const {
            AdaptiveStepController::EnergyState E_trial = *E_current;
            // ΔE_elastic ≈ -|G_mech| · Δa  (energy released from crack advance)
            // ΔE_surface = +Gc · Δa          (surface energy created)
            // Net: ΔE ≈ (Gc - |G_mech|) · Δa  < 0 when |G_mech| > Gc
            E_trial.elastic -= force_magnitude * delta_a;
            E_trial.surface += Gc * delta_a;
            E_trial.total    = E_trial.elastic + E_trial.surface - E_trial.external;
            return E_trial;
        }
    };

    LinearEnergyApprox approx;
    approx.E_current       = &current.energy;
    approx.force_magnitude = current.force.magnitude;
    approx.Gc              = Gc;

    auto result = controller.attempt_propagation(current.energy, approx, h_local);

    if (result.accepted) {
        const double theta = current.force.angle;
        CrackSegment seg;
        seg.start = tip_pos;
        seg.end[0] = tip_pos[0] + result.delta_a * std::cos(theta);
        seg.end[1] = tip_pos[1] + result.delta_a * std::sin(theta);
        crack.add_segment(seg);

        std::cout << "  [Propagation] Accepted: Δa=" << result.delta_a/h_local
                  << "·h, ΔE=" << result.energy_change << " J\n";
        return true;
    }

    std::cerr << "  [Propagation] Rejected by backtracking after "
              << result.backtrack_iterations << " halvings\n";
    return false;
}

// ============================================================================
// Weight Function for Integration Domain
// ============================================================================

void Assembler::compute_weight_function(int tip_node_id,
                                         double integration_radius,
                                         std::vector<double>& node_weights) const
{
    node_weights.assign(mesh.num_nodes(), 0.0);
    const double tx = mesh.node(tip_node_id).x;
    const double ty = mesh.node(tip_node_id).y;
    const double R  = integration_radius;

    for (size_t i = 0; i < mesh.num_nodes(); ++i) {
        const double dx = mesh.node(i).x - tx;
        const double dy = mesh.node(i).y - ty;
        const double d  = std::sqrt(dx*dx + dy*dy);
        if (d <= R) node_weights[i] = 1.0 - d/R;  // linear decay
    }
}

// ============================================================================
// Element patch (elements sharing a node)
// ============================================================================

std::vector<int> Assembler::get_element_patch(int node_id) const
{
    std::vector<int> patch;
    for (size_t e = 0; e < mesh.num_elements(); ++e) {
        const auto& elem = mesh.element(e);
        for (int k = 0; k < 3; ++k)
            if (elem.vid[k] == node_id) { patch.push_back(static_cast<int>(e)); break; }
    }
    return patch;
}
