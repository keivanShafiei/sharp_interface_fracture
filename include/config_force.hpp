#pragma once

#include <Eigen/Dense>
#include <array>
#include <vector>

// Forward declarations
class Mesh;

/**
 * @brief Settings for configurational force computation
 * 
 * Based on Appendix B, Section B.3.1:
 * - duffy_order: Number of Gauss points in each direction (η,ζ) for Duffy transformation
 *                Default 7 achieves O(h^1.8) convergence (Fig. B.1 in paper)
 * - standard_order: Gauss order for elements away from crack tip (3 points sufficient for P1)
 * - tip_radius: Radius within which elements are treated as "crack-tip elements" 
 *               requiring Duffy transformation (typically 3h per paper)
 */
struct ConfigForceSettings {
    int duffy_order = 7;        // 7-point Gauss in each direction (η,ζ)
    int standard_order = 3;     // 3-point Gauss for non-tip elements
    double tip_radius = 0.0;    // Automatically set to 3*h_local if zero
    double regularization = 1e-12; // Numerical zero for singularity checks
};

/**
 * @brief Configurational force computation via discrete J-integral
 * 
 * Implements Appendix B: "Mechanics of Configurational Forces and Numerical Implementation"
 * 
 * Core equations:
 * - Eq. B.8:  Eshelby tensor Σ = W(∇u)I - (∇u)ᵀσ
 * - Eq. B.14: G_I = -∑_{K∈Star(I)} ∫_K Σ : ∇N_I dx
 * - Eq. B.18: G_tip = G_mech - G_c ν̂
 * - Eq. B.20-B.21: Duffy transformation for singularity treatment
 * 
 * Material assumptions (Section 2, Table 1):
 * - Plane strain (ε_zz = 0)
 * - Linear isotropic elasticity
 * - Small strain kinematics
 */
class ConfigurationalForce {
public:
    /**
     * @brief Construct configurational force evaluator
     * @param E Young's modulus [Pa]
     * @param nu Poisson's ratio [-]
     * @param plane_strain true for plane strain, false for plane stress
     * @param settings Numerical integration settings
     */
    ConfigurationalForce(double E, double nu, bool plane_strain,
                        ConfigForceSettings settings = ConfigForceSettings());

    /**
     * @brief Compute configurational force at crack tip
     * 
     * Implements Eq. B.18: G_tip = G_mech - G_c ν̂
     * 
     * @param mesh Computational mesh
     * @param displacement Global displacement field [u_x, u_y, u_x, u_y, ...]
     * @param tip_position Crack tip coordinates [x, y]
     * @param tip_normal Unit normal at crack tip pointing into uncracked region
     * @param G_c Critical energy release rate [J/m²]
     * @return Configurational force vector [G_x, G_y] [N/m] in 2D (force per unit thickness)
     */
    std::array<double, 2> compute_tip_force(
        const Mesh& mesh,
        const std::vector<double>& displacement,
        const std::array<double, 2>& tip_position,
        const std::array<double, 2>& tip_normal,
        double G_c) const;

    /**
     * @brief Compute element-wise configurational forces (for debugging/visualization)
     * 
     * Returns G_K for each element K (useful for checking patch independence)
     * 
     * @param mesh Computational mesh
     * @param displacement Global displacement field
     * @return Vector of element configurational forces [G_x, G_y] per element
     */
    std::vector<std::array<double, 2>> compute_element_forces(
        const Mesh& mesh,
        const std::vector<double>& displacement) const;

private:
    // Material properties (Table 1)
    double E_;           // Young's modulus
    double nu_;          // Poisson's ratio
    bool plane_strain_;  // true = plane strain (ε_zz=0), false = plane stress (σ_zz=0)
    
    // Derived material constants
    Eigen::Matrix3d C_;  // Voigt notation: 3×3 for 2D (xx, yy, xy components)
    
    // Numerical settings
    ConfigForceSettings settings_;

    /**
     * @brief Build elasticity tensor in Voigt notation
     * 
     * For plane strain (Eq. 2.1 in paper):
     * E' = E/(1-ν²)
     * C_11 = C_22 = E'(1-ν)
     * C_12 = E'ν
     * C_33 = E'(1-2ν)/2 = G
     * 
     * @return 3×3 matrix in Voigt notation [σ_xx, σ_yy, σ_xy]ᵀ = C [ε_xx, ε_yy, 2ε_xy]ᵀ
     */
    Eigen::Matrix3d build_elasticity_tensor() const;

    /**
     * @brief Compute strain tensor from displacement gradient
     * 
     * ε = ½(∇u + (∇u)ᵀ) (small strain assumption, Section 2)
     * 
     * @param grad_u 2×2 displacement gradient [∂u_i/∂x_j]
     * @return 3×1 strain in Voigt notation [ε_xx, ε_yy, 2ε_xy]ᵀ
     */
    Eigen::Vector3d compute_strain_voigt(const Eigen::Matrix2d& grad_u) const;

    /**
     * @brief Compute stress from strain via constitutive law
     * 
     * σ = C : ε (Eq. B.8 context)
     * 
     * @param strain_voigt Strain in Voigt notation
     * @return Stress in Voigt notation [σ_xx, σ_yy, σ_xy]ᵀ
     */
    Eigen::Vector3d compute_stress_voigt(const Eigen::Vector3d& strain_voigt) const;

    /**
     * @brief Compute strain energy density
     * 
     * W = ½ ε : C : ε = ½ σ : ε (Eq. B.8)
     * 
     * @param strain_voigt Strain in Voigt notation
     * @return Scalar strain energy density [J/m³]
     */
    double compute_strain_energy_density(const Eigen::Vector3d& strain_voigt) const;

    /**
     * @brief Compute Eshelby stress tensor
     * 
     * Σ = W·I - (∇u)ᵀ·σ (Eq. B.8)
     * 
     * Returns 2×2 Eshelby tensor (full tensor, not Voigt)
     * 
     * @param W Strain energy density
     * @param grad_u Displacement gradient 2×2
     * @param stress_voigt Stress in Voigt notation
     * @return 2×2 Eshelby stress tensor
     */
    Eigen::Matrix2d compute_eshelby_tensor(
        double W,
        const Eigen::Matrix2d& grad_u,
        const Eigen::Vector3d& stress_voigt) const;

    /**
     * @brief Integrate configurational force over a single element
     * 
     * G_K = -∫_K Σ : ∇N_I dx (Eq. B.14)
     * 
     * @param mesh Mesh
     * @param elem_id Element index
     * @param node_id Global node index I
     * @param displacement Global displacement field
     * @param use_duffy If true, apply Duffy transformation (Eq. B.20-B.21)
     * @return Contribution to configurational force from this element [G_x, G_y]
     */
    std::array<double, 2> integrate_element_contribution(
        const Mesh& mesh,
        int elem_id,
        int node_id,
        const std::vector<double>& displacement,
        bool use_duffy) const;

    /**
     * @brief Apply Duffy transformation for crack-tip element
     * 
     * Implements Eq. B.20-B.21:
     * ξ(η,ζ) = (1-η)X_tip + η(1-ζ)X_2 + ηζ X_3
     * ∫_K f dx = ∫₀¹∫₀¹ f(ξ) · η · |det J| dηdζ
     * 
     * The factor η cancels the r^{-1/2} singularity (Appendix B.3.1)
     * 
     * @param mesh Mesh
     * @param elem_id Element index (must be crack-tip element with tip at vertex 0)
     * @param node_id Node for which to compute force
     * @param displacement Global displacement field
     * @return Duffy-transformed integral contribution
     */
    std::array<double, 2> duffy_integrate_tip_element(
        const Mesh& mesh,
        int elem_id,
        int node_id,
        const std::vector<double>& displacement) const;

    /**
     * @brief Gauss quadrature points and weights in 1D
     * 
     * @param order Number of points (3 or 7 supported)
     * @param points Output quadrature points in [-1,1]
     * @param weights Output quadrature weights
     */
    void gauss_quadrature_1d(int order,
                            std::vector<double>& points,
                            std::vector<double>& weights) const;

    /**
     * @brief Check if element is a crack-tip element requiring Duffy transformation
     * 
     * An element is "tip element" if one of its vertices is within tip_radius of tip_position
     * 
     * @param mesh Mesh
     * @param elem_id Element index
     * @param tip_position Crack tip [x,y]
     * @return true if element requires Duffy treatment
     */
    bool is_tip_element(const Mesh& mesh,
                       int elem_id,
                       const std::array<double, 2>& tip_position) const;

    /**
     * @brief Find vertex index within element that is closest to tip
     * 
     * Used to determine which vertex should be mapped to origin in Duffy transformation
     * 
     * @param mesh Mesh
     * @param elem_id Element index
     * @param tip_position Crack tip [x,y]
     * @return Local vertex index (0, 1, or 2 for triangle)
     */
    int find_tip_vertex_local(const Mesh& mesh,
                             int elem_id,
                             const std::array<double, 2>& tip_position) const;
};
