#ifndef ASSEMBLER_HPP
#define ASSEMBLER_HPP

#include "mesh.hpp"
#include "material.hpp"
#include "crack.hpp"  // اضافه شد
#include "adaptive_stepping.hpp"  // اضافه شد
#include <petsc.h>
#include <vector>
#include <array>
#include <cmath>

struct DirichletBC {
    int dof;
    double val;
};

struct ConfigurationalForce {
    std::array<double, 2> force_vector;
    double magnitude;
    double angle;
    int tip_node_id;
    
    ConfigurationalForce() : force_vector(), magnitude(0.0), angle(0.0), tip_node_id(-1) {
        force_vector[0] = 0.0;
        force_vector[1] = 0.0;
    }
};

// اضافه شد: ساختار برای نتیجه force + energy
struct TipForceResult {
    ConfigurationalForce force;
    AdaptiveStepController::EnergyState energy;
};

class Assembler {
public:
    Assembler(const Mesh& mesh_, const Material& mat_, double thickness);

    // Existing methods
    PetscErrorCode assemble(Mat& K, Vec& f);
    PetscErrorCode solve_linear(const Mat& K, const Vec& f, std::vector<double>& u) const;
    void compute_element_stress(const std::vector<double>& u, std::vector<std::array<double,3> >& sigma_out) const;

    void apply_dirichlet_on_bottom(double uy_value);
    void apply_dirichlet_on_top(double uy_value);
    void apply_dirichlet_pin_leftbottom_x();
    void apply_dirichlet_at_node(int node_idx, int dim, double val);
    void clear_boundary_conditions() { dbc.clear(); }  // Clear BCs before new step
    
    double compute_reaction_force(const std::vector<int>& node_indices, const std::vector<double>& u_global) const;
    
    ConfigurationalForce compute_tip_force(int tip_node_id, 
                                           const std::vector<double>& u_global,
                                           const std::array<double, 2>& crack_normal) const;
    
    void compute_element_eshelby_tensor(int elem_id,
                                        const std::vector<double>& u_global,
                                        std::array<double, 4>& eshelby_tensor) const;
    
    void compute_weight_function(int tip_node_id,
                                 double integration_radius,
                                 std::vector<double>& node_weights) const;
    
    std::vector<int> get_element_patch(int node_id) const;
    
    double compute_element_strain_energy_density(int elem_id,
                                                  const std::vector<double>& u_global) const;

    // ✅ NEW: Enhanced methods for backtracking support
    TipForceResult compute_tip_force_with_energy(
        int tip_node_id,
        const std::vector<double>& u_global,
        const std::array<double, 2>& crack_normal) const;
    
    double compute_total_elastic_energy(const std::vector<double>& u_global) const;
    double compute_surface_energy(const CrackPath& crack) const;
    double compute_external_work(const std::vector<double>& u_global) const;
    
    bool propagate_with_backtracking(
        CrackPath& crack,
        const std::vector<double>& u_global,
        AdaptiveStepController& controller,
        double h_local);

private:
    void element_stiffness_CST(int e, std::array<double,36>& ke, double& area) const;
    
    const Mesh& mesh;
    const Material& mat;
    double t;
    std::vector<DirichletBC> dbc;
};

#endif // ASSEMBLER_HPP
