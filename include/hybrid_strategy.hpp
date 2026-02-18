#ifndef HYBRID_STRATEGY_HPP
#define HYBRID_STRATEGY_HPP

#include "mesh.hpp"
#include "crack.hpp"
#include <vector>
#include <array>
#include <cmath>
#include <limits>

class HybridStrategy {
public:
    struct HybridParams {
        double l_c_factor;
        double h_pf_factor;
        double init_extension_ratio;
        double path_hausdorff_tol;
        bool verbose;
        
        // C++11 compatible constructor
        HybridParams()
            : l_c_factor(5.0)
            , h_pf_factor(1.25)
            , init_extension_ratio(0.1)
            , path_hausdorff_tol(0.05)
            , verbose(true)
        {}
    };
    
    struct PathTransferResult {
        bool success;
        std::vector<std::array<double, 2> > sharp_path_nodes;
        double hausdorff_distance;
        double energy_difference;
        double displacement_l2_error;
        
        PathTransferResult()
            : success(false)
            , hausdorff_distance(0)
            , energy_difference(0)
            , displacement_l2_error(0)
        {}
    };
    
    explicit HybridStrategy(const HybridParams& params = HybridParams())
        : params_(params) {}
    
    bool should_use_hybrid(double K_II_over_K_I, bool user_request = false) const {
        return (K_II_over_K_I > 0.5) || user_request;
    }
    
    std::vector<double> phase_field_initialization(
        const Mesh& mesh,
        double h_bulk,
        double initial_crack_length,
        const std::string& bc_type
    ) {
        if (params_.verbose) {
            std::cout << "[Hybrid Phase 1] Coarse phase-field initialization\n";
            std::cout << "  ℓc = " << params_.l_c_factor * h_bulk << " m\n";
        }
        
        std::vector<double> damage_field(mesh.num_nodes(), 0.0);
        std::cout << "  ⚠️  Phase-field solver not implemented.\n";
        return damage_field;
    }

private:
    HybridParams params_;
    
    double point_to_segment_distance(
        const std::array<double, 2>& p,
        const std::array<double, 2>& a,
        const std::array<double, 2>& b
    ) const {
        double dx = b[0] - a[0];
        double dy = b[1] - a[1];
        double len_sq = dx*dx + dy*dy;
        
        if (len_sq < 1e-12) {
            double dpx = p[0] - a[0];
            double dpy = p[1] - a[1];
            return std::sqrt(dpx*dpx + dpy*dpy);
        }
        
        double t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len_sq;
        t = std::max(0.0, std::min(1.0, t));
        
        double proj_x = a[0] + t * dx;
        double proj_y = a[1] + t * dy;
        
        double dpx = p[0] - proj_x;
        double dpy = p[1] - proj_y;
        
        return std::sqrt(dpx*dpx + dpy*dpy);
    }
};

#endif // HYBRID_STRATEGY_HPP
