// Material.hpp extension - Add this to your Material class

#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <array>
#include <cmath>

class Material {
public:
    Material(double E_, double nu_, bool plane_strain_ = true, double Gc_ = 100.0)
        : E(E_), nu(nu_), plane_strain(plane_strain_), Gc(Gc_) {}
    
    // Existing method - returns constitutive matrix D
    std::array<double, 9> Dmatrix() const {
        std::array<double, 9> D = {0};
        
        if (plane_strain) {
            double factor = E / ((1.0 + nu) * (1.0 - 2.0*nu));
            D[0] = factor * (1.0 - nu);  // D11
            D[1] = factor * nu;           // D12
            D[2] = 0.0;                   // D13 (shear coupling)
            D[3] = factor * nu;           // D21
            D[4] = factor * (1.0 - nu);  // D22
            D[5] = 0.0;                   // D23
            D[6] = 0.0;                   // D31
            D[7] = 0.0;                   // D32
            D[8] = factor * (1.0 - 2.0*nu) / 2.0;  // D33 (shear)
        } else {
            // Plane stress
            double factor = E / (1.0 - nu*nu);
            D[0] = factor;               // D11
            D[1] = factor * nu;          // D12
            D[2] = 0.0;
            D[3] = factor * nu;          // D21
            D[4] = factor;               // D22
            D[5] = 0.0;
            D[6] = 0.0;
            D[7] = 0.0;
            D[8] = factor * (1.0 - nu) / 2.0;  // D33
        }
        
        return D;
    }
    
    // NEW: Accessor for fracture toughness (needed for configurational forces)
    double get_fracture_toughness() const { return Gc; }
    
    // NEW: Compute stress intensity factor from energy release rate
    // K_I = sqrt(E' * G) where E' = E/(1-ν²) for plane strain
    double compute_stress_intensity_factor(double G_computed) const {
        double E_prime = plane_strain ? E / (1.0 - nu*nu) : E;
        return std::sqrt(E_prime * G_computed);
    }
    
    // NEW: Compute critical load for LEFM
    // For SENT: P_crit = sqrt(E' * Gc * a0 / π) / f(a0/W)
    double compute_critical_load(double crack_length, double width, 
                                  double thickness = 1.0) const {
        double E_prime = plane_strain ? E / (1.0 - nu*nu) : E;
        double a_over_W = crack_length / width;
        
        // Geometric correction factor for SENT (Tada handbook)
        double f = 1.12 - 0.23*a_over_W + 10.6*std::pow(a_over_W, 2) 
                   - 21.7*std::pow(a_over_W, 3) + 30.4*std::pow(a_over_W, 4);
        
        return std::sqrt(E_prime * Gc * crack_length / M_PI) * thickness / f;
    }
    
    double get_youngs_modulus() const { return E; }
    double get_poisson_ratio() const { return nu; }
    bool is_plane_strain() const { return plane_strain; }
    
private:
    double E;             // Young's modulus (Pa)
    double nu;            // Poisson's ratio
    bool plane_strain;    // true = plane strain, false = plane stress
    double Gc;            // Fracture toughness (J/m²)
};

#endif // MATERIAL_HPP
