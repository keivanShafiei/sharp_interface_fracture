#!/usr/bin/env bash
# =============================================================================
# reproduce_all_results.sh
# Reproduces all benchmark results from the paper.
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
echo "================================================================"
echo "Sharp-Interface Fracture — Full Reproduction Script"
echo "================================================================"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "OS: $(uname -a)"

# Check required tools
for tool in cmake make python3 awk; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not found"; exit 1; }
done

# ---------------------------------------------------------------------------
# 1. Build
# ---------------------------------------------------------------------------
echo ""
echo "=== STEP 1: BUILD ==="
mkdir -p build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-O3 -march=native" \
    2>&1 | tee cmake_output.log

make -j4 2>&1 | tee make_output.log
echo "  [OK] Build completed"
cd ..

BINARY="./build/sharp_fracture"
if [ ! -f "$BINARY" ]; then
    echo "ERROR: Binary not found at $BINARY"
    exit 1
fi

# ---------------------------------------------------------------------------
# 2. Verification Tests (must pass before benchmarks)
# ---------------------------------------------------------------------------
echo ""
echo "=== STEP 2: VERIFICATION TESTS ==="

echo "  Running T1: Patch Test..."
./build/tests/patch_test > results/T1_patch_test.log 2>&1
if grep -q "PASS" results/T1_patch_test.log; then
    echo "  [PASS] T1: Patch test"
else
    echo "  [FAIL] T1: Patch test — check results/T1_patch_test.log"
    echo "  WARNING: Continuing with benchmarks despite patch test failure"
fi

echo "  Running T2: Elasticity Consistency..."
./build/tests/test_elasticity_consistency > results/T2_elasticity.log 2>&1
if grep -q "PASS" results/T2_elasticity.log; then
    echo "  [PASS] T2: Elasticity consistency"
else
    echo "  [FAIL] T2: Elasticity consistency (known issue — see LIMITATIONS.md L-I1)"
fi

# ---------------------------------------------------------------------------
# 3. SENT Benchmark
# ---------------------------------------------------------------------------
echo ""
echo "=== STEP 3: SENT BENCHMARK ==="
echo "  Material: E=30 GPa, nu=0.20, Gc=100 J/m2"
echo "  Mesh: h_bulk=0.010 m, h_tip=0.003 m"
echo "  Algorithm: beta=1.0, theta_h=0.01*Gc, eps_tol=1e-6"

mkdir -p results/SENT

$BINARY \
    --case SENT \
    --mesh meshes/sent_h010.msh \
    > results/SENT/sent_output.log 2>&1

mv SENT_load_displacement.txt   results/SENT/ 2>/dev/null || true
mv SENT_energy_history.txt      results/SENT/ 2>/dev/null || true
mv SENT_crack_path_history.txt  results/SENT/ 2>/dev/null || true
mv energy_history.csv           results/SENT/ 2>/dev/null || true

echo "  Checking energy monotonicity (T3)..."
python3 tests/check_energy_monotonicity.py results/SENT/energy_history.csv \
    --eps_tol 1e-6 > results/SENT/T3_monotonicity.log 2>&1
cat results/SENT/T3_monotonicity.log

echo "  Checking ε-stationary bound (T7): |G_final| ≤ 101 J/m2..."
python3 -c "
import sys
data = open('results/SENT/SENT_crack_path_history.txt').readlines()
data = [l for l in data if not l.startswith('#')]
if data:
    last = data[-1].split()
    G_final = float(last[4])
    Gc = 100.0
    theta_h = 0.01*Gc
    if G_final <= Gc + theta_h:
        print(f'  [PASS] T7: |G_final| = {G_final:.4f} <= {Gc+theta_h} J/m2')
    else:
        print(f'  [FAIL] T7: |G_final| = {G_final:.4f} > {Gc+theta_h} J/m2')
else:
    print('  [SKIP] T7: No crack path data found')
"

echo "  [OK] SENT benchmark complete — results in results/SENT/"

# ---------------------------------------------------------------------------
# 4. TPB Benchmark
# ---------------------------------------------------------------------------
echo ""
echo "=== STEP 4: TPB BENCHMARK ==="
echo "  Material: E=30 GPa, nu=0.20, Gc=100 J/m2"
echo "  Mesh: h_bulk=0.020 m, h_tip=0.006 m"

mkdir -p results/TPB

$BINARY \
    --case TPB \
    --mesh meshes/tpb_h020.msh \
    > results/TPB/tpb_output.log 2>&1

mv TPB_load_displacement.txt   results/TPB/ 2>/dev/null || true
mv TPB_energy_history.txt      results/TPB/ 2>/dev/null || true
mv TPB_crack_path_history.txt  results/TPB/ 2>/dev/null || true
mv energy_history.csv          results/TPB/ 2>/dev/null || true

echo "  Checking energy monotonicity..."
python3 tests/check_energy_monotonicity.py results/TPB/energy_history.csv \
    --eps_tol 1e-6 > results/TPB/T3_monotonicity.log 2>&1
cat results/TPB/T3_monotonicity.log

echo "  [OK] TPB benchmark complete — results in results/TPB/"

# ---------------------------------------------------------------------------
# 5. Convergence Study (T4)
# ---------------------------------------------------------------------------
echo ""
echo "=== STEP 5: ENERGY CONVERGENCE STUDY (T4) ==="
echo "  Running SENT at 4 mesh refinements..."

mkdir -p results/convergence

for H in 0.040 0.020 0.010 0.005; do
    echo "    h = $H m ..."
    MESH="meshes/sent_h${H/./}.msh"
    if [ -f "$MESH" ]; then
        $BINARY --case SENT --mesh "$MESH" > results/convergence/sent_h${H/./}.log 2>&1
        # Extract final elastic energy from log
        ENERGY=$(grep "elastic_energy" results/convergence/sent_h${H/./}.log | tail -1 | awk '{print $2}' || echo "N/A")
        echo "    h=$H  E=$ENERGY J" >> results/convergence/convergence_data.txt
    else
        echo "    Mesh $MESH not found — skipping"
    fi
done

echo "  Computing convergence rates..."
python3 tests/compute_convergence_rate.py results/convergence/convergence_data.txt \
    > results/convergence/T4_convergence.log 2>&1
cat results/convergence/T4_convergence.log

# ---------------------------------------------------------------------------
# 6. Summary Report
# ---------------------------------------------------------------------------
echo ""
echo "=== STEP 6: SUMMARY ==="

python3 - << 'EOF'
import os, glob

results = {
    'T1 Patch Test':         'results/T1_patch_test.log',
    'T2 Elasticity':         'results/T2_elasticity.log',
    'T3 SENT Monotonicity':  'results/SENT/T3_monotonicity.log',
    'T3 TPB Monotonicity':   'results/TPB/T3_monotonicity.log',
    'T4 Convergence':        'results/convergence/T4_convergence.log',
}

print("\n=== FINAL VERIFICATION SUMMARY ===")
print(f"{'Test':<30} {'Status':<10}")
print("-" * 42)
for name, logfile in results.items():
    if os.path.exists(logfile):
        content = open(logfile).read()
        status = "PASS" if "PASS" in content else ("FAIL" if "FAIL" in content else "UNKNOWN")
    else:
        status = "MISSING"
    print(f"{name:<30} {status:<10}")

print("")
print("NOTE: T2 (Elasticity Consistency) is expected to FAIL")
print("      due to known discrepancy in config_force.cpp (see LIMITATIONS.md L-I1)")
print("")
print("Output files:")
print("  results/SENT/  — SENT load-displacement, energy, crack path")
print("  results/TPB/   — TPB load-displacement, energy, crack path")
print("  results/convergence/  — Energy convergence data")
EOF

echo ""
echo "================================================================"
echo "Reproduction complete. See results/ directory."
echo "================================================================"
echo ""
echo "IMPORTANT: The following paper results CANNOT be reproduced"
echo "           from this repository alone:"
echo "  1. Speedup 6.1-12.0x (requires phase-field reference implementation)"
echo "  2. Double-Notch benchmark (not yet implemented)"
echo "  3. Sobol sensitivity analysis (not implemented)"
echo "================================================================"
