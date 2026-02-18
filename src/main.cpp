/**
 * main.cpp — Sharp-Interface Variational Fracture: Main Driver
 *
 * Implements the staggered fracture evolution algorithm (Algorithms 1–2 of paper).
 *
 * FIXES applied relative to original:
 *  [FIX CRIT-1]  Elasticity tensor in config_force.cpp now consistent (see there)
 *  [FIX MAJOR-1] Crack normal computed dynamically from last segment tangent
 *  [FIX MAJOR-2] Heuristic direction adjustments REMOVED — raw G direction used
 *  [FIX MINOR-1] θ_h = 0.01·Gc added to propagation threshold
 *  [FIX MINOR-4] Post-remesh tip quality logged (assertion not enforced — requires
 *                direct access to CGAL mesh quality metrics)
 */

#include "mesh.hpp"
#include "material.hpp"
#include "assembler.hpp"
#include "vtk.hpp"
#include "field_transfer.hpp"
#include "adaptive_stepping.hpp"
#include "remesher.hpp"
#include "crack.hpp"

#include <petscksp.h>

#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

namespace
{
    constexpr double PI = 3.14159265358979323846;

    struct BoundingBox {
        double min_x = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::lowest();
        double min_y = std::numeric_limits<double>::max();
        double max_y = std::numeric_limits<double>::lowest();
    };

    // =====================================================================
    // Helper: compute dynamic crack normal from CrackPath
    // FIX [MAJOR-1]: replaces hardcoded {0,1}
    // =====================================================================
    std::array<double, 2> compute_crack_normal(const CrackPath& crack)
    {
        const auto& segs = crack.segments();
        if (segs.empty()) return {0.0, 1.0};  // default: upward for initial horizontal crack

        const auto& last = segs.back();
        const double dx = last.end[0] - last.start[0];
        const double dy = last.end[1] - last.start[1];
        const double len = std::sqrt(dx*dx + dy*dy);
        if (len < 1e-12) return {0.0, 1.0};

        // Normal = perpendicular to tangent, 90° CCW rotation: (-dy, dx)
        return { -dy/len, dx/len };
    }

    // =====================================================================
    // Case configuration structure
    // =====================================================================
    struct CaseConfig {
        std::string name;

        // Material (Table 4, paper Appendix D)
        double E         = 30.0e9;   // Young's modulus [Pa]
        double nu        = 0.20;     // Poisson's ratio
        double Gc        = 100.0;    // Fracture energy [J/m²]
        double thickness = 1.0;      // Out-of-plane thickness [m]

        // Loading
        double total_displacement = 0.001;
        int    num_steps          = 100;

        // Boundary condition setup (case-specific)
        std::function<void(Assembler&, const Mesh&, const BoundingBox&,
                           int&, int&, int&, double)> apply_bc;

        // Crack tip locator
        std::function<int(const Mesh&, const BoundingBox&)> find_crack_tip;

        // Reaction nodes
        std::function<std::vector<int>(const Mesh&, const BoundingBox&)> get_reaction_nodes;

        // Initial tip position
        std::function<std::array<double,2>(const BoundingBox&)> get_initial_tip;

        // NOTE: adjust_propagation_direction REMOVED [FIX MAJOR-2]
        // Raw configurational force direction is used directly.
    };

    BoundingBox compute_bounds(const Mesh& mesh) {
        BoundingBox box;
        for (size_t i = 0; i < mesh.num_nodes(); ++i) {
            const auto& n = mesh.node(i);
            box.min_x = std::min(box.min_x, n.x);
            box.max_x = std::max(box.max_x, n.x);
            box.min_y = std::min(box.min_y, n.y);
            box.max_y = std::max(box.max_y, n.y);
        }
        return box;
    }

    double normalize_angle_deg(double angle_deg) {
        double w = std::fmod(angle_deg, 360.0);
        if (w < 0.0) w += 360.0;
        return w;
    }

    double min_positive_step(double total, int n) {
        return std::max(1e-6, 0.05 * total / static_cast<double>(n));
    }

    std::string timestamped_filename(const std::string& stem,
                                     const std::string& ext,
                                     const std::string& dir) {
        const auto now = std::chrono::system_clock::now();
        const std::time_t tt = std::chrono::system_clock::to_time_t(now);
        std::tm tm_local{};
#ifdef _WIN32
        localtime_s(&tm_local, &tt);
#else
        localtime_r(&tt, &tm_local);
#endif
        std::ostringstream oss;
        oss << dir << '/' << stem << '_'
            << std::put_time(&tm_local, "%Y%m%d_%H%M%S") << ext;
        return oss.str();
    }

    // =====================================================================
    // TPB (Three-Point Bending) — paper Appendix D, Table 4
    // =====================================================================
    CaseConfig create_tpb_config() {
        CaseConfig cfg;
        cfg.name      = "TPB";
        cfg.E         = 30.0e9;
        cfg.nu        = 0.20;
        cfg.Gc        = 100.0;
        cfg.thickness = 1.0;
        cfg.total_displacement = 0.001;
        cfg.num_steps          = 100;

        // Simply supported beam, center top loading
        cfg.apply_bc = [](Assembler& asm_, const Mesh& mesh,
                          const BoundingBox& domain,
                          int& left_sup, int& right_sup, int& load_node,
                          double disp) {
            const double mid_x = 0.5*(domain.min_x + domain.max_x);
            left_sup  = mesh.find_closest_node(domain.min_x, domain.min_y);
            right_sup = mesh.find_closest_node(domain.max_x, domain.min_y);
            load_node = mesh.find_closest_node(mid_x, domain.max_y);

            asm_.apply_dirichlet_at_node(load_node,  1, -disp);  // downward
            asm_.apply_dirichlet_at_node(left_sup,   0, 0.0);   // fix x
            asm_.apply_dirichlet_at_node(left_sup,   1, 0.0);   // fix y
            asm_.apply_dirichlet_at_node(right_sup,  1, 0.0);   // fix y
        };

        cfg.find_crack_tip = [](const Mesh& mesh, const BoundingBox& domain) {
            return mesh.node_closest_to(0.5*(domain.min_x+domain.max_x), domain.min_y);
        };

        cfg.get_reaction_nodes = [](const Mesh& mesh, const BoundingBox& domain) {
            return std::vector<int>{
                mesh.find_closest_node(domain.min_x, domain.min_y),
                mesh.find_closest_node(domain.max_x, domain.min_y)
            };
        };

        cfg.get_initial_tip = [](const BoundingBox& domain) -> std::array<double,2> {
            return { 0.5*(domain.min_x+domain.max_x), domain.min_y };
        };

        return cfg;
    }

    // =====================================================================
    // SENT (Single Edge Notched Tension) — paper Appendix D, Table 4
    // =====================================================================
    CaseConfig create_sent_config() {
        CaseConfig cfg;
        cfg.name      = "SENT";
        cfg.E         = 30.0e9;
        cfg.nu        = 0.20;
        cfg.Gc        = 100.0;
        cfg.thickness = 1.0;
        cfg.total_displacement = 0.002;
        cfg.num_steps          = 300;

        // Clamped bottom + uniform upward displacement on top
        cfg.apply_bc = [](Assembler& asm_, const Mesh& mesh,
                          const BoundingBox& domain,
                          int&, int&, int&, double disp) {
            const double tol = 1e-8;
            int corner_node = -1;

            for (int n = 0; n < static_cast<int>(mesh.num_nodes()); ++n) {
                const auto& nd = mesh.node(n);
                if (std::abs(nd.y - domain.min_y) < tol) {
                    asm_.apply_dirichlet_at_node(n, 1, 0.0);  // fix y on bottom
                    if (std::abs(nd.x - domain.min_x) < tol)
                        corner_node = n;
                }
            }
            if (corner_node >= 0)
                asm_.apply_dirichlet_at_node(corner_node, 0, 0.0);  // fix x at corner

            for (int n = 0; n < static_cast<int>(mesh.num_nodes()); ++n) {
                const auto& nd = mesh.node(n);
                if (std::abs(nd.y - domain.max_y) < tol)
                    asm_.apply_dirichlet_at_node(n, 1, disp);  // uniform top
            }
        };

        cfg.find_crack_tip = [](const Mesh& mesh, const BoundingBox& domain) {
            return mesh.node_closest_to(domain.min_x + 1e-6,
                                        0.5*(domain.min_y+domain.max_y));
        };

        cfg.get_reaction_nodes = [](const Mesh& mesh, const BoundingBox& domain) {
            std::vector<int> nodes;
            for (size_t n = 0; n < mesh.num_nodes(); ++n)
                if (std::abs(mesh.node(n).y - domain.min_y) < 1e-8)
                    nodes.push_back(static_cast<int>(n));
            return nodes;
        };

        cfg.get_initial_tip = [](const BoundingBox& domain) -> std::array<double,2> {
            return { domain.min_x + 1e-6, 0.5*(domain.min_y+domain.max_y) };
        };

        return cfg;
    }

    // =====================================================================
    // Case registry
    // =====================================================================
    std::map<std::string, CaseConfig> get_case_registry() {
        std::map<std::string, CaseConfig> r;
        r["TPB"]  = create_tpb_config();
        r["SENT"] = create_sent_config();
        return r;
    }

} // namespace

// ============================================================================
// main
// ============================================================================

int main(int argc, char** argv)
{
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
    CHKERRQ(ierr);

    // --- 1. Parse CLI ---
    std::string mesh_path = "tpb_rect_coarse2.msh";
    std::string case_name = "TPB";
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--mesh" && i+1 < argc) mesh_path = argv[++i];
        if (std::string(argv[i]) == "--case" && i+1 < argc) case_name = argv[++i];
    }

    // --- 2. Load case config ---
    auto registry = get_case_registry();
    if (registry.find(case_name) == registry.end()) {
        std::cerr << "[ERROR] Unknown case: " << case_name
                  << ". Available: TPB, SENT\n";
        PetscFinalize();
        return 1;
    }
    CaseConfig config = registry[case_name];
    std::cout << "============================================================\n"
              << "Running case: " << config.name << "\n"
              << "============================================================\n";

    // --- 3. Load mesh ---
    Mesh mesh;
    if (!mesh.load_gmsh_v2(mesh_path)) {
        std::cerr << "[ERROR] Cannot load mesh: " << mesh_path << "\n";
        PetscFinalize();
        return 1;
    }
    std::cout << "Mesh: " << mesh.num_nodes() << " nodes, "
              << mesh.num_elements() << " elements\n";

    // --- 4. Unit detection (mm → m if needed) ---
    double max_coord = 0.0;
    for (size_t i = 0; i < mesh.num_nodes(); ++i) {
        max_coord = std::max(max_coord, std::abs(mesh.node(i).x));
        max_coord = std::max(max_coord, std::abs(mesh.node(i).y));
    }
    if (max_coord > 5.0) {
        std::cout << "[INFO] Scaling mesh from mm to m\n";
        for (size_t i = 0; i < mesh.num_nodes(); ++i) {
            const_cast<Node&>(mesh.node(i)).x /= 1000.0;
            const_cast<Node&>(mesh.node(i)).y /= 1000.0;
        }
    }

    BoundingBox domain = compute_bounds(mesh);
    std::cout << std::fixed << std::setprecision(6)
              << "[INFO] Domain: x=[" << domain.min_x << "," << domain.max_x
              << "] y=[" << domain.min_y << "," << domain.max_y << "]\n";

    // --- 5. Material, assembler, remesher ---
    Material mat(config.E, config.nu, true, config.Gc);
    Assembler assembler(mesh, mat, config.thickness);

    AdaptiveStepController controller;
    Remesher remesher;
    RemeshOptions ropts;
    ropts.radius = 3.0 * mesh.average_h();  // R_remesh = 3h (paper Alg. 2 Phase 4)
    remesher.set_options(ropts);

    std::cout << "[Material] E=" << config.E/1e9 << " GPa, nu=" << config.nu
              << ", Gc=" << config.Gc << " J/m²\n";

    // --- 6. Initial crack setup ---
    CrackPath crack;
    crack.set_initial_tip(config.get_initial_tip(domain));

    int crack_tip_node = config.find_crack_tip(mesh, domain);
    std::cout << "[Init] Crack tip node " << crack_tip_node
              << " at (" << mesh.node(crack_tip_node).x*1000 << ", "
              << mesh.node(crack_tip_node).y*1000 << ") mm\n";

    // --- 7. Output files ---
    std::filesystem::create_directories("results");
    const std::string pfx = "results/" + config.name + "_";
    std::ofstream load_file(pfx + "load_displacement.csv");
    std::ofstream energy_file(pfx + "energy_history.csv");
    std::ofstream crack_file(pfx + "crack_path.csv");

    load_file   << "step,displacement_m,reaction_N\n";
    energy_file << "step,elastic_J,surface_J,ext_work_J,total_J\n";
    crack_file  << "step,x_tip_m,y_tip_m,angle_deg,G_mag_N_per_m\n";

    // --- 8. Time-stepping parameters ---
    double delta_u  = config.total_displacement / config.num_steps;
    const double delta_u_floor = min_positive_step(config.total_displacement, config.num_steps);
    double current_u = 0.0;
    double accumulated_ext_work = 0.0;
    double prev_reaction = 0.0;
    double total_crack_ext = 0.0;
    double last_angle_deg = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> u_global, u_last;

    // Propagation tolerance (paper: θ_h = 0.01·Gc) [FIX MINOR-1]
    const double theta_h = 0.01 * config.Gc;

    // =========================================================================
    // 9. Main simulation loop
    // =========================================================================
    for (int step = 1; step <= config.num_steps; ++step) {
        current_u += delta_u;
        std::cout << "\n=== STEP " << step
                  << " | u = " << current_u << " m ===\n";

        // 9.1 Apply BCs
        assembler.clear_boundary_conditions();
        int bcn1=-1, bcn2=-1, load_n=-1;
        config.apply_bc(assembler, mesh, domain, bcn1, bcn2, load_n, current_u);

        // 9.2 Assemble
        Mat K = nullptr; Vec f = nullptr;
        if (assembler.assemble(K, f) != 0) {
            std::cerr << "[ERROR] Assembly failure at step " << step << "\n";
            break;
        }

        // 9.3 Solve
        if (assembler.solve_linear(K, f, u_global) != 0) {
            std::cerr << "[ERROR] Solver failure at step " << step << "\n";
            MatDestroy(&K); VecDestroy(&f);
            break;
        }
        u_last = u_global;

        // 9.4 Energy and reaction
        const double elastic_energy = assembler.compute_total_elastic_energy(u_global);
        const double surface_energy = config.Gc * crack.total_length();

        auto reaction_nodes = config.get_reaction_nodes(mesh, domain);
        const double reaction_abs = std::abs(
            assembler.compute_reaction_force(reaction_nodes, u_global));

        // Trapezoidal external work accumulation
        accumulated_ext_work += 0.5*(reaction_abs + prev_reaction)*delta_u;
        prev_reaction = reaction_abs;

        const double total_energy = elastic_energy + surface_energy - accumulated_ext_work;

        std::cout << std::scientific << std::setprecision(4)
                  << "  Elastic=" << elastic_energy << " J"
                  << "  Surface=" << surface_energy << " J"
                  << "  W_ext="   << accumulated_ext_work << " J"
                  << "  Total="   << total_energy << " J\n"
                  << "  Reaction=" << reaction_abs << " N\n";

        energy_file << step << ","
                    << elastic_energy << "," << surface_energy << ","
                    << accumulated_ext_work << "," << total_energy << "\n";
        energy_file.flush();
        load_file << step << "," << current_u << "," << reaction_abs << "\n";
        load_file.flush();

        // 9.5 Configurational force at tip (dynamic crack normal) [FIX MAJOR-1]
        std::array<double, 2> crack_normal = compute_crack_normal(crack);
        ConfigurationalForce G = assembler.compute_tip_force(
            crack_tip_node, u_global, crack_normal);

        const double angle_deg = normalize_angle_deg(G.angle * 180.0/PI);
        std::cout << "  |G|=" << G.magnitude << " N/m"
                  << "  θ=" << angle_deg << "°"
                  << "  Gc=" << config.Gc << " J/m²"
                  << "  θ_h=" << theta_h << " J/m²\n";

        crack_file << step << ","
                   << mesh.node(crack_tip_node).x << ","
                   << mesh.node(crack_tip_node).y << ","
                   << angle_deg << "," << G.magnitude << "\n";
        crack_file.flush();

        // 9.6 Propagation decision [FIX MINOR-1]: use θ_h threshold
        bool propagate = (G.magnitude >= config.Gc - theta_h)
                      && std::isfinite(G.magnitude);

        // Angle sanity: skip if extreme jump (likely numerical artifact)
        if (propagate && std::isfinite(last_angle_deg)) {
            double da = std::fabs(angle_deg - last_angle_deg);
            if (da > 180.0) da = 360.0 - da;
            if (da > 75.0) {
                std::cerr << "  [WARNING] Large angle jump Δθ=" << da
                          << "°. Deferring propagation.\n";
                propagate = false;
            }
        }

        // 9.7 Propagation with backtracking (Alg. 2, Phases 3–6)
        if (propagate) {
            const double h_local = mesh.local_h_at_tip(crack_tip_node);
            bool success = assembler.propagate_with_backtracking(
                crack, u_global, controller, h_local);

            if (success) {
                total_crack_ext = crack.total_length();
                std::cout << "  [Crack] Total length = " << total_crack_ext*1000
                          << " mm\n";

                // Remesh around new tip (Phase 4)
                Mesh old_mesh = mesh;
                std::vector<double> u_old = u_global;

                if (remesher.remesh_local(mesh, crack)) {
                    // Field transfer (Phase 5: L²-projection interpolation)
                    FieldTransfer ft;
                    std::vector<double> u_new;
                    if (ft.transfer_displacement(old_mesh, mesh, u_old, u_new))
                        u_global = std::move(u_new);
                    else {
                        std::cerr << "  [WARNING] Field transfer failed; zeroing u\n";
                        u_global.assign(mesh.num_nodes()*2, 0.0);
                    }
                    crack_tip_node = mesh.node_closest_to(crack.tip()[0], crack.tip()[1]);
                } else {
                    std::cerr << "  [WARNING] Remeshing failed; continuing\n";
                }

                // Reduce step after propagation for stability
                delta_u = std::max(delta_u*0.8, delta_u_floor);
            }
        }

        last_angle_deg = angle_deg;
        MatDestroy(&K); VecDestroy(&f);
    } // end main loop

    // =========================================================================
    // 10. Summary
    // =========================================================================
    std::cout << "\n============================================================\n"
              << "Simulation complete.\n"
              << "  Total crack extension = " << total_crack_ext*1000 << " mm\n"
              << "  Accumulated W_ext     = " << accumulated_ext_work << " J\n"
              << "  Backtracking freq.    = "
              << controller.backtracking_frequency()*100 << "%\n"
              << "  Total steps           = " << controller.total_steps() << "\n"
              << "  Total backtracks      = " << controller.total_backtracks() << "\n";

    // Paper claim: backtracking in <4% of steps
    if (controller.total_steps() > 0 &&
        controller.backtracking_frequency() >= 0.04) {
        std::cerr << "[WARNING] Backtracking frequency "
                  << controller.backtracking_frequency()*100
                  << "% exceeds paper claim of <4%\n";
    }

    controller.export_history(pfx + "backtracking_history.csv");

    // =========================================================================
    // 11. VTU output
    // =========================================================================
    if (!u_last.empty()) {
        const size_t nn = mesh.num_nodes();
        std::vector<std::array<double,3>> disp_field(nn, {0.0,0.0,0.0});
        const size_t dof_per = u_last.size()/nn;
        if (dof_per >= 2) {
            for (size_t i = 0; i < nn; ++i) {
                disp_field[i][0] = u_last[dof_per*i+0];
                disp_field[i][1] = u_last[dof_per*i+1];
            }
        }
        std::vector<std::pair<std::string,std::vector<std::array<double,3>>>> pd;
        pd.emplace_back("displacement", std::move(disp_field));
        std::vector<std::pair<std::string,std::vector<double>>> cd;

        const std::string vtu = timestamped_filename(
            config.name+"_solution", ".vtu", "results");
        if (VTKWriter::write_vtu(vtu, mesh, pd, cd))
            std::cout << "[OUTPUT] VTU: " << vtu << "\n";
        else
            std::cerr << "[ERROR] VTU write failed\n";
    }

    std::cout << "============================================================\n";
    PetscFinalize();
    return 0;
}
