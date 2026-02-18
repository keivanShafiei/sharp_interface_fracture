#!/usr/bin/env python3
"""
compute_convergence_rate.py
Reads convergence data (h, E_h) and computes log-log slope.
Expected rate for CST elastic energy: ≥ 1.8 (paper claims O(h²)).
"""

import sys
import math
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_file", help="Two-column file: h  energy")
    parser.add_argument("--expected_rate", type=float, default=1.8,
                        help="Minimum expected convergence rate")
    args = parser.parse_args()

    data = []
    with open(args.data_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    h = float(parts[0])
                    E = float(parts[1])
                    data.append((h, E))
                except ValueError:
                    continue

    if len(data) < 2:
        print("[SKIP] Not enough data points for convergence analysis")
        sys.exit(0)

    data.sort(key=lambda x: x[0])  # sort by h ascending

    # Reference: finest mesh energy (or analytical if provided)
    E_ref = data[0][1]

    print("h [m]       |E_h - E_ref|   Rate")
    print("-" * 42)

    rates = []
    prev_h = None
    prev_err = None

    for h, E in data:
        err = abs(E - E_ref)
        rate_str = "---"
        if prev_h is not None and prev_err is not None and err > 0 and prev_err > 0:
            rate = math.log(prev_err/err) / math.log(prev_h/h)
            rates.append(rate)
            rate_str = f"{rate:.3f}"
        print(f"{h:.5f}     {err:.4e}       {rate_str}")
        prev_h, prev_err = h, err

    if rates:
        avg_rate = sum(rates) / len(rates)
        print(f"\nAverage convergence rate: {avg_rate:.3f}")
        if avg_rate >= args.expected_rate:
            print(f"[PASS] Rate {avg_rate:.3f} ≥ expected {args.expected_rate}")
            sys.exit(0)
        else:
            print(f"[FAIL] Rate {avg_rate:.3f} < expected {args.expected_rate}")
            sys.exit(1)
    else:
        print("[SKIP] Not enough points to compute rate")
        sys.exit(0)


if __name__ == "__main__":
    main()
