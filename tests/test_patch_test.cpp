/**
 * test_patch_test.cpp — Patch Test for CST Element Correctness
 *
 * Verifies that the stiffness assembly produces correct displacements
 * for a known analytical solution (constant stress state).
 *
 * Test setup:
 *   - Square domain [0,1]×[0,1], plane strain
 *   - E = 1 Pa, ν = 0.3, t = 1 m, no crack
 *   - Applied: u_x = ε₀·x everywhere (uniform tension)
 *   - Expected σ_xx = E(1-ν)/((1+ν)(1-2ν)) · ε₀ (constant)
 *
 * Acceptance: |σ_computed - σ_expected| / σ_expected < 1e-8 (machine precision for CST)
 */

#include "assembler.hpp"
#include "material.hpp"
#include "mesh.hpp"
#include "crack.hpp"
#include <cmath>
#include <iostream>
#include <petscksp.h>

// Create a simple 2×2 structured triangular mesh on [0,1]×[0,1]
// 9 nodes, 8 elements (2 triangles per square cell)
static Mesh make_patch_mesh()
{
    Mesh m;
    // 3×3 nodes
    for (int j = 0; j <= 2; ++j)
        for (int i = 0; i <= 2; ++i)
            m.add_node(0.5*i, 0.5*j);

    //   6-7-8
    //   |/|/|
    //   3-4-5
    //   |/|/|
    //   0-1-2
    auto idx = [](int i, int j){ return j*3+i; };

    // Lower-left triangles (BL) + upper-right (UR) for each cell
    for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < 2; ++i) {
            int n0=idx(i,j), n1=idx(i+1,j), n2=idx(i,j+1), n3=idx(i+1,j+1);
            m.add_element(n0, n1, n2);
            m.add_element(n1, n3, n2);
        }
    }
    return m;
}

int main(int argc, char** argv)
{
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
    CHKERRQ(ierr);

    std::cout << "=== PATCH TEST ===\n";

    Mesh mesh = make_patch_mesh();
    const double E_test  = 1.0;
    const double nu_test = 0.3;
    const double eps0    = 0.01;   // applied strain

    Material mat(E_test, nu_test, true, 100.0);
    Assembler asm_(mesh, mat, 1.0);

    // Apply uniform tension: u_x = eps0*x, u_y = 0 everywhere (prescribed)
    for (int n = 0; n < static_cast<int>(mesh.num_nodes()); ++n) {
        const double x = mesh.node(n).x;
        const double y = mesh.node(n).y;
        asm_.apply_dirichlet_at_node(n, 0, eps0*x);  // u_x = ε₀·x
        asm_.apply_dirichlet_at_node(n, 1, 0.0);     // u_y = 0
    }

    Mat K; Vec f;
    if (asm_.assemble(K, f) != 0) {
        std::cerr << "[FAIL] Assembly error\n";
        PetscFinalize();
        return 1;
    }

    std::vector<double> u;
    if (asm_.solve_linear(K, f, u) != 0) {
        std::cerr << "[FAIL] Solver error\n";
        MatDestroy(&K); VecDestroy(&f);
        PetscFinalize();
        return 1;
    }

    MatDestroy(&K); VecDestroy(&f);

    // Compute element stresses
    std::vector<std::array<double,3>> sigma;
    asm_.compute_element_stress(u, sigma);

    // Expected σ_xx for plane strain under ε_xx = ε₀, ε_yy = 0
    // σ_xx = (λ + 2μ)·ε₀  where λ=Eν/((1+ν)(1-2ν)), μ=E/(2(1+ν))
    const double lambda = E_test * nu_test / ((1.0+nu_test)*(1.0-2.0*nu_test));
    const double mu     = E_test / (2.0*(1.0+nu_test));
    const double sxx_expected = (lambda + 2.0*mu) * eps0;

    double max_err = 0.0;
    for (size_t e = 0; e < sigma.size(); ++e) {
        const double err = std::abs(sigma[e][0] - sxx_expected);
        if (err > max_err) max_err = err;
    }
    const double rel_err = max_err / std::abs(sxx_expected);

    std::cout << "  σ_xx expected = " << sxx_expected << " Pa\n";
    std::cout << "  σ_xx max error = " << max_err << " Pa  (rel: " << rel_err << ")\n";

    const bool pass = rel_err < 1e-8;
    std::cout << (pass ? "[PASS] Patch test\n" : "[FAIL] Patch test\n");

    // Check σ_yy (should be λ·ε₀)
    const double syy_expected = lambda * eps0;
    double max_err_y = 0.0;
    for (auto& s : sigma) max_err_y = std::max(max_err_y, std::abs(s[1]-syy_expected));
    std::cout << "  σ_yy expected = " << syy_expected << " Pa, max err = " << max_err_y << "\n";

    // Check σ_xy ≈ 0 (no shear expected)
    double max_shear = 0.0;
    for (auto& s : sigma) max_shear = std::max(max_shear, std::abs(s[2]));
    std::cout << "  σ_xy max = " << max_shear << " (should be ~0)\n";
    const bool pass_all = pass && max_err_y/std::max(std::abs(syy_expected),1e-15) < 1e-6
                                && max_shear < 1e-10;

    PetscFinalize();
    return pass_all ? 0 : 1;
}
