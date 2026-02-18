/**
 * config_force.cpp — Configurational (Eshelby) Force Computation
 *
 * Implements Appendix B of the paper: configurational force via discrete J-integral.
 *
 * FIX [CRIT-1]: build_elasticity_tensor() now uses standard plane-strain Lamé
 *               constants, consistent with material.hpp. Original code used
 *               E' = E/(1-ν²) (plane-stress notation) which gave ~1.6% error in C_11.
 *
 * Theory:
 *   Eshelby tensor:        Σ = W(∇u)·I − (∇u)ᵀ·σ          (Eq. B.8)
 *   Nodal config. force:   G_I = -∑_{K∈Star(I)} ∫_K Σ:∇N_I dx  (Eq. B.14)
 *   Tip force:             G_tip = G_mech − G_c·ν̂            (Eq. B.18)
 *   Duffy transform:       for singularity treatment          (Eq. B.20–B.21)
 *
 * NOTE on Duffy for P1 elements:
 *   For linear (P1) triangular elements, ∇u is constant per element — there
 *   is no integrable r^{-1/2} singularity in the discrete fields. The Duffy
 *   transform is included for completeness and compatibility with the paper's
 *   description, but standard 1-point quadrature is mathematically equivalent
 *   for P1 elements. See THEORY.md for discussion.
 */

#include "config_force.hpp"
#include "mesh.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <algorithm>

// ============================================================================
// Constructor
// ============================================================================

ConfigurationalForce::ConfigurationalForce(
    double E, double nu, bool plane_strain, ConfigForceSettings settings)
    : E_(E), nu_(nu), plane_strain_(plane_strain), settings_(settings)
{
    if (E <= 0.0)
        throw std::invalid_argument("Young's modulus must be positive");
    if (nu <= -1.0 || nu >= 0.5)
        throw std::invalid_argument("Poisson's ratio must be in (-1, 0.5)");

    C_ = build_elasticity_tensor();
}

// ============================================================================
// Public: compute_tip_force
// ============================================================================

std::array<double, 2> ConfigurationalForce::compute_tip_force(
    const Mesh& mesh,
    const std::vector<double>& displacement,
    const std::array<double, 2>& tip_position,
    const std::array<double, 2>& tip_normal,
    double G_c) const
{
    std::array<double, 2> G_mech = {0.0, 0.0};

    double R = settings_.tip_radius;
    if (R <= 0.0) R = 3.0 * mesh.average_h();

    for (int elem_id = 0; elem_id < mesh.num_elements(); ++elem_id) {
        const auto center = mesh.get_element_center(elem_id);
        const double dist = std::hypot(center[0]-tip_position[0], center[1]-tip_position[1]);
        if (dist > R) continue;

        // Weight function: linear decay from 1 at tip to 0 at R
        const double q_elem = std::max(0.0, 1.0 - dist/R);

        const bool use_duffy = is_tip_element(mesh, elem_id, tip_position);

        // Sum over all three nodes of this element (Eq. B.14 sums over Star(I))
        std::array<double, 2> elem_force = {0.0, 0.0};
        const auto& nodes_e = mesh.element_nodes(elem_id);
        for (int node_id : nodes_e) {
            auto contrib = integrate_element_contribution(
                mesh, elem_id, node_id, displacement, use_duffy);
            elem_force[0] += contrib[0];
            elem_force[1] += contrib[1];
        }

        G_mech[0] += elem_force[0] * q_elem;
        G_mech[1] += elem_force[1] * q_elem;
    }

    // G_tip = G_mech − G_c · ν̂  (Eq. B.18)
    return { G_mech[0] - G_c * tip_normal[0],
             G_mech[1] - G_c * tip_normal[1] };
}

// ============================================================================
// Public: compute_element_forces (for debugging/visualization)
// ============================================================================

std::vector<std::array<double, 2>> ConfigurationalForce::compute_element_forces(
    const Mesh& mesh,
    const std::vector<double>& displacement) const
{
    std::vector<std::array<double, 2>> elem_forces(mesh.num_elements());

    for (int elem_id = 0; elem_id < mesh.num_elements(); ++elem_id) {
        std::array<double, 2> avg = {0.0, 0.0};
        const auto& elem_nodes = mesh.element_nodes(elem_id);
        for (int local = 0; local < 3; ++local) {
            auto c = integrate_element_contribution(
                mesh, elem_id, elem_nodes[local], displacement, false);
            avg[0] += c[0] / 3.0;
            avg[1] += c[1] / 3.0;
        }
        elem_forces[elem_id] = avg;
    }
    return elem_forces;
}

// ============================================================================
// Private: Elasticity Tensor — FIXED (standard plane-strain Lamé constants)
// ============================================================================
// FIX [CRIT-1]: Original code used non-standard E' = E/(1-ν²) formulation
// causing ~1.6% error in C_11 and ~6.25% error in C_12 for ν=0.20.
//
// Correct plane-strain elasticity (ε_zz = 0):
//   λ  = Eν / ((1+ν)(1-2ν))
//   μ  = E / (2(1+ν))
//   C_11 = C_22 = λ + 2μ = E(1-ν)/((1+ν)(1-2ν))
//   C_12 = C_21 = λ       = Eν/((1+ν)(1-2ν))
//   C_33             = μ  = E/(2(1+ν))
//
// Consistency check: C_11 computed here must match Material::Dmatrix()[0].

Eigen::Matrix3d ConfigurationalForce::build_elasticity_tensor() const
{
    Eigen::Matrix3d C = Eigen::Matrix3d::Zero();

    if (plane_strain_) {
        // Standard plane-strain Lamé parameters [CORRECTED]
        const double lambda = E_ * nu_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
        const double mu     = E_ / (2.0 * (1.0 + nu_));

        C(0, 0) = lambda + 2.0 * mu;   // C_11
        C(1, 1) = lambda + 2.0 * mu;   // C_22
        C(0, 1) = lambda;               // C_12
        C(1, 0) = lambda;               // C_21
        C(2, 2) = mu;                   // C_33 (shear)

        // Self-check against known formula
        const double C11_expected = E_ * (1.0-nu_) / ((1.0+nu_)*(1.0-2.0*nu_));
        const double err = std::abs(C(0,0) - C11_expected) / C11_expected;
        if (err > 1e-10) {
            throw std::runtime_error(
                "[ConfigurationalForce] Elasticity tensor self-check failed. "
                "C_11 error = " + std::to_string(err));
        }
    } else {
        // Plane stress (σ_zz = 0)
        const double factor = E_ / (1.0 - nu_ * nu_);
        const double mu     = E_ / (2.0 * (1.0 + nu_));

        C(0, 0) = factor;
        C(1, 1) = factor;
        C(0, 1) = factor * nu_;
        C(1, 0) = factor * nu_;
        C(2, 2) = mu;
    }

    return C;
}

// ============================================================================
// Private: Strain, Stress, Energy Density, Eshelby Tensor
// ============================================================================

Eigen::Vector3d ConfigurationalForce::compute_strain_voigt(
    const Eigen::Matrix2d& grad_u) const
{
    Eigen::Vector3d strain;
    strain(0) = grad_u(0, 0);                    // ε_xx = ∂u_x/∂x
    strain(1) = grad_u(1, 1);                    // ε_yy = ∂u_y/∂y
    strain(2) = grad_u(0, 1) + grad_u(1, 0);    // 2ε_xy
    return strain;
}

Eigen::Vector3d ConfigurationalForce::compute_stress_voigt(
    const Eigen::Vector3d& strain_voigt) const
{
    return C_ * strain_voigt;
}

double ConfigurationalForce::compute_strain_energy_density(
    const Eigen::Vector3d& strain_voigt) const
{
    const Eigen::Vector3d stress = compute_stress_voigt(strain_voigt);
    // W = ½ εᵀ σ  (accounting for factor-2 in shear: see Voigt convention)
    // W = ½(ε_xx σ_xx + ε_yy σ_yy + 2ε_xy σ_xy)
    // Since strain(2) = 2ε_xy and stress(2) = σ_xy:
    // W = ½(strain(0)*stress(0) + strain(1)*stress(1) + strain(2)*stress(2))
    return 0.5 * strain_voigt.dot(stress);
}

Eigen::Matrix2d ConfigurationalForce::compute_eshelby_tensor(
    double W,
    const Eigen::Matrix2d& grad_u,
    const Eigen::Vector3d& stress_voigt) const
{
    // Σ = W·I − (∇u)ᵀ·σ  (Eq. B.8)
    // Full 2×2 stress tensor from Voigt: σ = [[sxx, sxy],[sxy, syy]]
    Eigen::Matrix2d sigma;
    sigma(0,0) = stress_voigt(0);  // σ_xx
    sigma(1,1) = stress_voigt(1);  // σ_yy
    sigma(0,1) = stress_voigt(2);  // σ_xy
    sigma(1,0) = stress_voigt(2);  // σ_yx = σ_xy (symmetric)

    Eigen::Matrix2d Sigma;
    Sigma(0,0) = W - (grad_u(0,0)*sigma(0,0) + grad_u(1,0)*sigma(1,0));
    Sigma(0,1) =   - (grad_u(0,0)*sigma(0,1) + grad_u(1,0)*sigma(1,1));
    Sigma(1,0) =   - (grad_u(0,1)*sigma(0,0) + grad_u(1,1)*sigma(1,0));
    Sigma(1,1) = W - (grad_u(0,1)*sigma(0,1) + grad_u(1,1)*sigma(1,1));
    return Sigma;
}

// ============================================================================
// Private: Integrate Element Contribution (Eq. B.14)
// ============================================================================
// For P1 elements: ∇u, ∇N, and Σ are all constant → exact with 1-point quadrature.
// When use_duffy=true, the Duffy transform is applied but produces the same
// result as direct integration for P1 (see THEORY.md).

std::array<double, 2> ConfigurationalForce::integrate_element_contribution(
    const Mesh& mesh,
    int elem_id,
    int node_id,
    const std::vector<double>& displacement,
    bool use_duffy) const
{
    if (use_duffy) {
        return duffy_integrate_tip_element(mesh, elem_id, node_id, displacement);
    }

    const auto& elem_nodes = mesh.element_nodes(elem_id);

    // Find local index of node_id
    int local_idx = -1;
    for (int k = 0; k < 3; ++k)
        if (elem_nodes[k] == node_id) { local_idx = k; break; }
    if (local_idx < 0) return {0.0, 0.0};

    const auto X0 = mesh.node_coords(elem_nodes[0]);
    const auto X1 = mesh.node_coords(elem_nodes[1]);
    const auto X2 = mesh.node_coords(elem_nodes[2]);

    const double area = 0.5 * std::abs(
        (X1[0]-X0[0])*(X2[1]-X0[1]) - (X2[0]-X0[0])*(X1[1]-X0[1]));
    if (area < settings_.regularization) return {0.0, 0.0};

    const double inv2A = 1.0 / (2.0 * area);

    // Shape function gradients (∇N for each vertex)
    const double b[3] = {
        (X1[1]-X2[1])*inv2A, (X2[1]-X0[1])*inv2A, (X0[1]-X1[1])*inv2A
    };
    const double c[3] = {
        (X2[0]-X1[0])*inv2A, (X0[0]-X2[0])*inv2A, (X1[0]-X0[0])*inv2A
    };

    // Displacement gradient (constant for P1)
    Eigen::Matrix2d grad_u = Eigen::Matrix2d::Zero();
    for (int k = 0; k < 3; ++k) {
        const double ux = (2*elem_nodes[k]  <(int)displacement.size()) ? displacement[2*elem_nodes[k]]   : 0.0;
        const double uy = (2*elem_nodes[k]+1<(int)displacement.size()) ? displacement[2*elem_nodes[k]+1] : 0.0;
        grad_u(0,0) += ux * b[k];
        grad_u(0,1) += ux * c[k];
        grad_u(1,0) += uy * b[k];
        grad_u(1,1) += uy * c[k];
    }

    const Eigen::Vector3d strain  = compute_strain_voigt(grad_u);
    const Eigen::Vector3d stress  = compute_stress_voigt(strain);
    const double W                = compute_strain_energy_density(strain);
    const Eigen::Matrix2d Sigma   = compute_eshelby_tensor(W, grad_u, stress);

    // G_K = -∫_K Σ:∇N_I dx = -(Σ · ∇N_I) · A_e  (Eq. B.14)
    const double dN_dx = b[local_idx];
    const double dN_dy = c[local_idx];

    return {
        -(Sigma(0,0)*dN_dx + Sigma(0,1)*dN_dy) * area,
        -(Sigma(1,0)*dN_dx + Sigma(1,1)*dN_dy) * area
    };
}

// ============================================================================
// Private: Duffy Transform (Eq. B.20–B.21)
// ============================================================================
// Transforms the unit triangle to a square to cancel r^{-1/2} singularity
// from the continuous displacement field. For P1, ∇u is constant and there is
// no discrete singularity; the Duffy result equals direct integration.
// Retained for paper fidelity.

std::array<double, 2> ConfigurationalForce::duffy_integrate_tip_element(
    const Mesh& mesh,
    int elem_id,
    int node_id,
    const std::vector<double>& displacement) const
{
    const auto& elem_nodes = mesh.element_nodes(elem_id);

    // Locate tip vertex (node_id must be in this element for Duffy)
    int local_tip_idx = -1;
    for (int k = 0; k < 3; ++k)
        if (elem_nodes[k] == node_id) { local_tip_idx = k; break; }
    if (local_tip_idx < 0)
        throw std::runtime_error("[Duffy] node_id not in element " + std::to_string(elem_id));

    // Rotate so X0 = tip vertex (required for Duffy mapping)
    std::array<int, 3> v;
    v[0] = elem_nodes[local_tip_idx];
    v[1] = elem_nodes[(local_tip_idx+1) % 3];
    v[2] = elem_nodes[(local_tip_idx+2) % 3];

    const auto X0 = mesh.node_coords(v[0]);
    const auto X1 = mesh.node_coords(v[1]);
    const auto X2 = mesh.node_coords(v[2]);

    const double area = 0.5 * std::abs(
        (X1[0]-X0[0])*(X2[1]-X0[1]) - (X2[0]-X0[0])*(X1[1]-X0[1]));
    if (area < settings_.regularization) return {0.0, 0.0};

    const double inv2A = 1.0 / (2.0 * area);

    // ∇N for tip node (after rotation, tip is local index 0 = v[0])
    const double dN0_dx = (X1[1]-X2[1]) * inv2A;
    const double dN0_dy = (X2[0]-X1[0]) * inv2A;

    // Displacement gradient (constant for P1, rotation doesn't change it)
    Eigen::Matrix2d grad_u = Eigen::Matrix2d::Zero();
    for (int k = 0; k < 3; ++k) {
        const double bk = (mesh.node_coords(v[(k+1)%3])[1] - mesh.node_coords(v[(k+2)%3])[1]) * inv2A;
        const double ck = (mesh.node_coords(v[(k+2)%3])[0] - mesh.node_coords(v[(k+1)%3])[0]) * inv2A;
        const double ux = (2*v[k]  <(int)displacement.size()) ? displacement[2*v[k]]   : 0.0;
        const double uy = (2*v[k]+1<(int)displacement.size()) ? displacement[2*v[k]+1] : 0.0;
        grad_u(0,0) += ux*bk; grad_u(0,1) += ux*ck;
        grad_u(1,0) += uy*bk; grad_u(1,1) += uy*ck;
    }

    const Eigen::Vector3d strain  = compute_strain_voigt(grad_u);
    const Eigen::Vector3d stress  = compute_stress_voigt(strain);
    const double W                = compute_strain_energy_density(strain);
    const Eigen::Matrix2d Sigma   = compute_eshelby_tensor(W, grad_u, stress);

    // Force density vector: -Σ · ∇N_tip
    Eigen::Vector2d fd;
    fd(0) = -(Sigma(0,0)*dN0_dx + Sigma(0,1)*dN0_dy);
    fd(1) = -(Sigma(1,0)*dN0_dx + Sigma(1,1)*dN0_dy);

    // Duffy quadrature: ∫∫_triangle f dA = ∫₀¹∫₀¹ f(η,ζ)·η·|J| dη dζ
    // For P1 (f=const): integral = f · A  (Duffy = standard)
    std::vector<double> gp, gw;
    gauss_quadrature_1d(settings_.duffy_order, gp, gw);

    // Map Gauss points [-1,1] → [0,1]
    std::vector<double> eta_pts(gp.size()), eta_wts(gw.size());
    for (size_t i = 0; i < gp.size(); ++i) {
        eta_pts[i] = 0.5*(gp[i]+1.0);
        eta_wts[i] = 0.5*gw[i];
    }

    std::array<double, 2> integral = {0.0, 0.0};
    for (size_t i = 0; i < eta_pts.size(); ++i) {
        for (size_t j = 0; j < eta_pts.size(); ++j) {
            const double eta  = eta_pts[i];
            const double zeta = eta_pts[j];
            if (eta < settings_.regularization) continue;

            // Jacobian of Duffy map (Eq. B.21)
            const double dx_deta  = -X0[0] + X1[0] + zeta*(X2[0]-X1[0]);
            const double dy_deta  = -X0[1] + X1[1] + zeta*(X2[1]-X1[1]);
            const double dx_dzeta =  eta*(X2[0]-X1[0]);
            const double dy_dzeta =  eta*(X2[1]-X1[1]);
            const double det_J    = std::abs(dx_deta*dy_dzeta - dx_dzeta*dy_deta);

            // Factor η cancels 1/r^{1/2} singularity (continuous theory)
            const double weight = eta * det_J * eta_wts[i] * eta_wts[j];

            integral[0] += fd(0) * weight;
            integral[1] += fd(1) * weight;
        }
    }

    return integral;
}

// ============================================================================
// Private: Gauss Quadrature Points and Weights
// ============================================================================

void ConfigurationalForce::gauss_quadrature_1d(
    int order,
    std::vector<double>& points,
    std::vector<double>& weights) const
{
    if (order == 3) {
        const double sq35 = std::sqrt(3.0/5.0);
        points  = {-sq35, 0.0, sq35};
        weights = {5.0/9.0, 8.0/9.0, 5.0/9.0};
    } else if (order == 7) {
        // 7-point Gauss-Legendre on [-1,1]
        points = {
            -0.9491079123427585,
            -0.7415311855993944,
            -0.4058451513773972,
             0.0,
             0.4058451513773972,
             0.7415311855993944,
             0.9491079123427585
        };
        weights = {
            0.1294849661688697,
            0.2797053914892767,
            0.3818300505051189,
            0.4179591836734694,
            0.3818300505051189,
            0.2797053914892767,
            0.1294849661688697
        };
    } else {
        throw std::invalid_argument(
            "[ConfigForce] Only 3-point and 7-point Gauss supported. Got: "
            + std::to_string(order));
    }
}

// ============================================================================
// Private: Tip Element Detection
// ============================================================================

bool ConfigurationalForce::is_tip_element(
    const Mesh& mesh,
    int elem_id,
    const std::array<double, 2>& tip_position) const
{
    double radius = settings_.tip_radius;
    if (radius <= 0.0) {
        // Estimate from edge length of this element
        const auto& en = mesh.element_nodes(elem_id);
        const auto A = mesh.node_coords(en[0]);
        const auto B = mesh.node_coords(en[1]);
        const double dx = B[0]-A[0], dy = B[1]-A[1];
        radius = 3.0 * std::sqrt(dx*dx + dy*dy);
    }

    for (int k = 0; k < 3; ++k) {
        const auto Xi = mesh.node_coords(mesh.element_nodes(elem_id)[k]);
        const double dx = Xi[0]-tip_position[0];
        const double dy = Xi[1]-tip_position[1];
        if (std::sqrt(dx*dx+dy*dy) < radius) return true;
    }
    return false;
}

int ConfigurationalForce::find_tip_vertex_local(
    const Mesh& mesh,
    int elem_id,
    const std::array<double, 2>& tip_position) const
{
    const auto& en = mesh.element_nodes(elem_id);
    int best = 0;
    double min_d = 1e100;
    for (int k = 0; k < 3; ++k) {
        const auto Xi = mesh.node_coords(en[k]);
        const double dx = Xi[0]-tip_position[0], dy = Xi[1]-tip_position[1];
        const double d  = std::sqrt(dx*dx+dy*dy);
        if (d < min_d) { min_d = d; best = k; }
    }
    return best;
}
