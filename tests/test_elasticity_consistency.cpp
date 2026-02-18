/**
 * test_elasticity_consistency.cpp
 *
 * Verifies that the elasticity tensor in ConfigurationalForce::build_elasticity_tensor()
 * is consistent with Material::Dmatrix() for plane strain.
 *
 * This test was FAILING before the fix (1.56% error in C_11, 6.25% in C_12).
 * After FIX [CRIT-1] both tensors must agree to machine precision.
 *
 * Acceptance: |C_ij_configforce - C_ij_material| / C_ij_material < 1e-8
 */

#include "config_force.hpp"
#include "material.hpp"
#include "mesh.hpp"
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

// Helper: expose C_ via a subclass (for testing only)
class ConfigForceTest : public ConfigurationalForce {
public:
    ConfigForceTest(double E, double nu, bool ps)
        : ConfigurationalForce(E, nu, ps) {}
    const Eigen::Matrix3d& get_C() const { return C_; }
};

static bool nearly_equal(double a, double b, double tol=1e-8) {
    const double denom = std::max(std::abs(b), 1e-15);
    return std::abs(a-b)/denom < tol;
}

int main()
{
    std::cout << "=== ELASTICITY TENSOR CONSISTENCY TEST ===\n\n";
    bool all_pass = true;

    // Test parameters: paper benchmark values
    const double E  = 30.0e9;
    const double nu = 0.20;

    Material mat(E, nu, true, 100.0);
    auto D_arr = mat.Dmatrix();

    // Material::Dmatrix() entries (plane strain)
    // D[0]=C11, D[1]=C12, D[4]=C22, D[8]=C33
    const double C11_material = D_arr[0];
    const double C12_material = D_arr[1];
    const double C22_material = D_arr[4];
    const double C33_material = D_arr[8];

    // Analytical expectations
    const double lambda   = E * nu / ((1.0+nu)*(1.0-2.0*nu));
    const double mu       = E / (2.0*(1.0+nu));
    const double C11_anal = lambda + 2.0*mu;
    const double C12_anal = lambda;
    const double C33_anal = mu;

    std::cout << "Reference (analytical plane strain):\n";
    std::cout << "  C11 = " << C11_anal/1e9 << " GPa  (expected: 33.333)\n";
    std::cout << "  C12 = " << C12_anal/1e9 << " GPa  (expected:  8.333)\n";
    std::cout << "  C33 = " << C33_anal/1e9 << " GPa  (expected: 12.500)\n\n";

    // Check Material::Dmatrix
    bool m_pass = nearly_equal(C11_material, C11_anal)
               && nearly_equal(C12_material, C12_anal)
               && nearly_equal(C33_material, C33_anal);
    std::cout << "Material::Dmatrix check: " << (m_pass ? "PASS" : "FAIL") << "\n";
    if (!m_pass) {
        std::cerr << "  C11=" << C11_material/1e9 << " C12=" << C12_material/1e9
                  << " C33=" << C33_material/1e9 << "\n";
        all_pass = false;
    }

    // ConfigurationalForce tensor (FIXED)
    ConfigForceTest cft(E, nu, true);
    const auto& C = cft.get_C();

    const double C11_cf = C(0,0);
    const double C12_cf = C(0,1);
    const double C22_cf = C(1,1);
    const double C33_cf = C(2,2);

    std::cout << "\nConfigurationalForce::C_ (post-fix):\n";
    std::cout << "  C11 = " << C11_cf/1e9 << " GPa\n";
    std::cout << "  C12 = " << C12_cf/1e9 << " GPa\n";
    std::cout << "  C22 = " << C22_cf/1e9 << " GPa\n";
    std::cout << "  C33 = " << C33_cf/1e9 << " GPa\n";

    bool c_pass = nearly_equal(C11_cf, C11_material)
               && nearly_equal(C12_cf, C12_material)
               && nearly_equal(C22_cf, C22_material)
               && nearly_equal(C33_cf, C33_material);

    std::cout << "\nConsistency check (ConfigForce vs Material): "
              << (c_pass ? "PASS" : "FAIL") << "\n";

    if (!c_pass) {
        std::cerr << "  Error C11: " << 100*(C11_cf-C11_material)/C11_material << "%\n";
        std::cerr << "  Error C12: " << 100*(C12_cf-C12_material)/std::abs(C12_material) << "%\n";
        all_pass = false;
    }

    // Symmetry check: C must be symmetric
    bool sym_pass = std::abs(C(0,1)-C(1,0))<1e-10
                 && std::abs(C(0,2)-C(2,0))<1e-10
                 && std::abs(C(1,2)-C(2,1))<1e-10;
    std::cout << "Symmetry check: " << (sym_pass ? "PASS" : "FAIL") << "\n";
    if (!sym_pass) all_pass = false;

    // Positive definiteness: all eigenvalues > 0
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);
    const double min_eig = es.eigenvalues().minCoeff();
    bool pd_pass = min_eig > 0.0;
    std::cout << "Positive definiteness (min eigenvalue=" << min_eig/1e9
              << " GPa): " << (pd_pass ? "PASS" : "FAIL") << "\n";
    if (!pd_pass) all_pass = false;

    // Also test ν = 0.4 (closer to 0.5, where the error was larger)
    {
        std::cout << "\n--- Test with nu=0.40 ---\n";
        const double nu2    = 0.40;
        const double lam2   = E * nu2 / ((1.0+nu2)*(1.0-2.0*nu2));
        const double mu2    = E / (2.0*(1.0+nu2));
        const double C11_a2 = lam2 + 2.0*mu2;

        ConfigForceTest cft2(E, nu2, true);
        const double C11_c2 = cft2.get_C()(0,0);
        const double err2   = 100*std::abs(C11_c2-C11_a2)/C11_a2;
        bool p2 = nearly_equal(C11_c2, C11_a2);
        std::cout << "  C11(nu=0.40): material=" << C11_a2/1e9
                  << " GPa, configforce=" << C11_c2/1e9
                  << " GPa, err=" << err2 << "%  " << (p2 ? "PASS" : "FAIL") << "\n";
        if (!p2) all_pass = false;
    }

    std::cout << "\n=== OVERALL: " << (all_pass ? "PASS" : "FAIL") << " ===\n";
    return all_pass ? 0 : 1;
}
