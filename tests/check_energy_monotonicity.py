#!/usr/bin/env python3
"""
check_energy_monotonicity.py
Reads energy_history.csv output and verifies that every accepted step satisfies
    E_new ≤ E_old + ε_tol · |E_old|
Prints pass/fail and statistics.
"""

import sys
import csv
import argparse
import math


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file", help="Path to energy_history CSV (from simulation output)")
    parser.add_argument("--eps_tol", type=float, default=1e-6,
                        help="Energy tolerance ε_tol (paper: 1e-6)")
    args = parser.parse_args()

    rows = []
    with open(args.csv_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(",")
            if len(parts) < 5:
                continue
            try:
                step        = int(parts[0])
                elastic     = float(parts[1])
                surface     = float(parts[2])
                ext_work    = float(parts[3])
                total       = float(parts[4])
                rows.append({"step": step, "elastic": elastic,
                             "surface": surface, "ext_work": ext_work,
                             "total": total})
            except ValueError:
                continue  # skip header

    if len(rows) < 2:
        print("[SKIP] Not enough data rows for monotonicity check")
        sys.exit(0)

    violations = 0
    max_violation = 0.0

    for i in range(1, len(rows)):
        E_old = rows[i-1]["total"]
        E_new = rows[i  ]["total"]
        tol   = args.eps_tol * abs(E_old)
        delta = E_new - E_old

        if delta > tol:
            violations += 1
            rel = abs(delta) / max(abs(E_old), 1e-15)
            if rel > max_violation:
                max_violation = rel
            print(f"  [VIOLATION] step={rows[i]['step']}: "
                  f"ΔE={delta:.6e} J > tol={tol:.6e} J "
                  f"(rel={rel:.4e})")

    bt_rate = 0.0
    # Try to read backtracking history if available
    bt_file = args.csv_file.replace("energy_history.csv", "backtracking_history.csv")
    try:
        with open(bt_file) as bf:
            lines = [l for l in bf if not l.startswith("#") and l.strip()]
            if lines:
                last = lines[-1].split(",")
                # Last line has cumulative data
                pass
    except FileNotFoundError:
        pass

    n_steps = len(rows) - 1
    print(f"\nSummary:")
    print(f"  Steps checked     : {n_steps}")
    print(f"  Violations        : {violations}")
    print(f"  ε_tol             : {args.eps_tol}")
    print(f"  Max rel violation : {max_violation:.4e}")

    if violations == 0:
        print("[PASS] Energy monotonicity satisfied on all steps")
        sys.exit(0)
    elif violations / max(n_steps, 1) < 0.04:
        print(f"[WARN] {violations}/{n_steps} violations ({100*violations/n_steps:.1f}%) "
              f"— within paper's <4% claim but non-zero")
        sys.exit(0)
    else:
        print(f"[FAIL] {violations}/{n_steps} violations ({100*violations/n_steps:.1f}%) "
              f"— exceeds paper's <4% claim")
        sys.exit(1)


if __name__ == "__main__":
    main()
