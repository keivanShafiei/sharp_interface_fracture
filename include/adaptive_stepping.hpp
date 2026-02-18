#ifndef ADAPTIVE_STEPPING_HPP
#define ADAPTIVE_STEPPING_HPP

#include <vector>
#include <array>
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>

/**
 * @brief Adaptive step size control with energy monotonicity enforcement
 * Compatible with C++11/14
 */
class AdaptiveStepController {
public:
    struct StepResult {
        bool accepted;
        double delta_a;
        double energy_change;
        int backtrack_iterations;
        double energy_residual;
        
        StepResult() : accepted(false), delta_a(0), energy_change(0),
                      backtrack_iterations(0), energy_residual(0) {}
    };
    
    struct EnergyState {
        double elastic;
        double surface;
        double external;
        double total;
        double dissipated;
        
        EnergyState() : elastic(0), surface(0), external(0), total(0), dissipated(0) {}
    };
    
    struct BacktrackingParams {
        double epsilon_tol;
        double beta_default;
        double beta_min;
        double beta_max;
        double reduction_factor;
        int max_iterations;
        double h_min_factor;
        bool enforce_strict;
        bool verbose;
        
        // Constructor with default values (C++11 compatible)
        BacktrackingParams() 
            : epsilon_tol(1e-6)
            , beta_default(1.0)
            , beta_min(0.5)
            , beta_max(2.0)
            , reduction_factor(0.5)
            , max_iterations(5)
            , h_min_factor(0.1)
            , enforce_strict(true)
            , verbose(false)
        {}
    };

    explicit AdaptiveStepController(const BacktrackingParams& params = BacktrackingParams())
        : params_(params), total_backtracks_(0), total_steps_(0) {}

    template<typename ComputeEnergyFunc>
    StepResult attempt_propagation(
        const EnergyState& E_old,
        ComputeEnergyFunc compute_new_state,
        double h_local
    ) {
        StepResult result;
        result.backtrack_iterations = 0;
        
        double delta_a = params_.beta_default * h_local;
        double h_min = params_.h_min_factor * h_local;
        
        for (int iter = 0; iter < params_.max_iterations; ++iter) {
            EnergyState E_new = compute_new_state(delta_a);
            double delta_E = E_new.total - E_old.total;
            double tolerance = params_.epsilon_tol * std::abs(E_old.total);
            
            bool energy_acceptable = (delta_E <= tolerance);
            
            if (!params_.enforce_strict) {
                energy_acceptable = (delta_E <= 10.0 * tolerance);
            }
            
            if (energy_acceptable) {
                result.accepted = true;
                result.delta_a = delta_a;
                result.energy_change = delta_E;
                result.energy_residual = std::abs(delta_E) / std::abs(E_old.total);
                
                if (params_.verbose && result.backtrack_iterations > 0) {
                    std::cout << "  [Backtrack] Accepted after " 
                              << result.backtrack_iterations << " halvings\n";
                }
                
                energy_history_.push_back(E_new);
                total_steps_++;
                return result;
            }
            
            delta_a *= params_.reduction_factor;
            result.backtrack_iterations++;
            total_backtracks_++;
            
            if (params_.verbose) {
                std::cout << "  [Backtrack] Iteration " << iter + 1 
                          << ": ΔE = " << delta_E << "\n";
            }
            
            if (delta_a < h_min) {
                result.accepted = false;
                result.delta_a = 0.0;
                result.energy_change = delta_E;
                return result;
            }
        }
        
        result.accepted = false;
        result.delta_a = 0.0;
        return result;
    }
    
    double compute_energy_residual(const EnergyState& state) const {
        double expected_total = state.elastic + state.surface;
        double residual = std::abs(state.external - expected_total);
        return residual / std::abs(state.external);
    }
    
    double backtracking_frequency() const {
        return total_steps_ > 0 ? 
            static_cast<double>(total_backtracks_) / total_steps_ : 0.0;
    }
    
    int total_backtracks() const { return total_backtracks_; }
    int total_steps() const { return total_steps_; }
    
    const std::vector<EnergyState>& history() const { return energy_history_; }
    
    void export_history(const std::string& filename) const {
        std::ofstream out(filename.c_str());
        out << "# Step, E_elastic, E_surface, E_external, E_total, Delta_E\n";
        
        double E_prev = 0.0;
        for (size_t i = 0; i < energy_history_.size(); ++i) {
            const EnergyState& E = energy_history_[i];
            double delta_E = (i > 0) ? (E.total - E_prev) : 0.0;
            
            out << i << ", "
                << E.elastic << ", "
                << E.surface << ", "
                << E.external << ", "
                << E.total << ", "
                << delta_E << "\n";
            
            E_prev = E.total;
        }
        out.close();
    }

    BacktrackingParams params_;  // Make public for testing

private:
    std::vector<EnergyState> energy_history_;
    int total_backtracks_;
    int total_steps_;
};

#endif // ADAPTIVE_STEPPING_HPP
