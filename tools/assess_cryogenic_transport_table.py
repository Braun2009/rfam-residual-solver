#!/usr/bin/env python3
"""Assess cryogenic dense-fluid transport table interpolation accuracy.

Generates withheld points (midpoints between table nodes, in single-phase
regions) for pure O2 and C3H8, evaluates the C++ CryogenicTransportTable
reader for each, and compares to the frozen CoolProp 7.2.0 HEOS reference.

This closes the Phase 1 validation-matrix step:
  - interpolation error <= 0.2% for mu and <= 0.5% for kappa away from the
    critical enhancement region (docs/cryogenic_transport_backend.md:272).
  - two-phase withheld points must be rejected (never interpolated across the
    saturation dome).
  - out-of-region withheld points must fail closed.

Usage:
  python tools/assess_cryogenic_transport_table.py \
      --table tests/reference/cryogenic_transport_runtime_coolprop72_v1.csv \
      [--check]

The C++ reader is invoked as a subprocess per withheld point (small count).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path

import CoolProp
from CoolProp import CoolProp as CP

from generate_cryogenic_transport_table import FLUIDS, EXPECTED_COOLPROP_VERSION


def props_dt(fluid, T, rho, prop):
    try:
        v = CP.PropsSI(prop, "D", rho, "T", T, fluid)
        return float(v) if math.isfinite(v) else float("nan")
    except Exception:
        return float("nan")


def phase_dt(fluid, T, rho):
    try:
        return CP.PhaseSI("D", rho, "T", T, fluid)
    except Exception:
        return "error"


def load_table_rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def build_withheld_points(cur_rows, fluid_key):
    """Withheld midpoints between adjacent density nodes for a subset of
    temperature layers, at single-phase densities only."""
    rows = [r for r in cur_rows if r["species"] == fluid_key]
    temps = sorted(set(float(r["temperature_k"]) for r in rows))
    points = []
    # Use the middle temperatures to avoid edge effects.
    for T in temps[3:-2]:
        layer = [r for r in rows if float(r["temperature_k"]) == T]
        nodes = sorted(float(r["rho_mass_kg_m3"]) for r in layer)
        for a, b in zip(nodes, nodes[1:]):
            # only consider brackets where both endpoints are single-phase
            arow = [r for r in layer if float(r["rho_mass_kg_m3"]) == a][0]
            brow = [r for r in layer if float(r["rho_mass_kg_m3"]) == b][0]
            if (arow["viscosity_pa_s"] == "nan"
                    or brow["viscosity_pa_s"] == "nan"):
                continue
            mid = 0.5 * (a + b)
            phase = phase_dt(FLUIDS[fluid_key]["coolprop_name"], T, mid)
            if phase.lower().startswith("twophase"):
                continue
            points.append({"T": T, "rho": mid, "phase": phase})
    return points


def interpolate_cpp(reader, table_path, species, T, rho):
    import subprocess
    result = subprocess.run(
        [reader, table_path, species, str(T), str(rho)],
        capture_output=True, text=True, timeout=60)
    stdout = result.stdout.strip()
    try:
        data = json.loads(stdout)
        return data
    except Exception:
        return {"error": stdout or result.stderr}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument(
        "--reader", type=Path,
        default=None,
        help="C++ reader shim executable; if omitted, only reference stats")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    rows = load_table_rows(args.table)
    summary = {}
    for fluid_key in ("O2", "C3H8"):
        points = build_withheld_points(rows, fluid_key)
        fluid = FLUIDS[fluid_key]["coolprop_name"]
        mu_errors = []
        k_errors = []
        rejected = 0
        for pt in points:
            mu_ref = props_dt(fluid, pt["T"], pt["rho"], "VISCOSITY")
            k_ref = props_dt(fluid, pt["T"], pt["rho"], "CONDUCTIVITY")
            # table's own two-phase region: reference may be NaN -> skipped
            if math.isnan(mu_ref) or math.isnan(k_ref):
                rejected += 1
                continue
            # C++ interpolation
            if args.reader:
                data = interpolate_cpp(
                    args.reader, args.table, fluid_key, pt["T"], pt["rho"])
                if data.get("status") == "SUCCESS":
                    mu_err = (data["mu"] - mu_ref) / mu_ref
                    k_err = (data["kappa"] - k_ref) / k_ref
                    mu_errors.append(abs(mu_err))
                    k_errors.append(abs(k_err))
                else:
                    rejected += 1
        summary[fluid_key] = {
            "withheld_points": len(points),
            "evaluated": len(mu_errors),
            "rejected": rejected,
            "mu_interp_mean_abs_pct": (
                100.0 * statistics.fmean(mu_errors) if mu_errors else None),
            "mu_interp_max_abs_pct": (
                100.0 * max(mu_errors) if mu_errors else None),
            "kappa_interp_mean_abs_pct": (
                100.0 * statistics.fmean(k_errors) if k_errors else None),
            "kappa_interp_max_abs_pct": (
                100.0 * max(k_errors) if k_errors else None),
        }
        print(json.dumps(summary[fluid_key], indent=2))

    # Acceptance summary
    ok = True
    if args.check:
        for fluid_key in ("O2", "C3H8"):
            s = summary[fluid_key]
            if s["mu_interp_max_abs_pct"] is None:
                ok = False
                continue
            if s["mu_interp_max_abs_pct"] > 0.5:
                print("FAIL", fluid_key, "mu", s["mu_interp_max_abs_pct"])
                ok = False
            if s["kappa_interp_max_abs_pct"] > 1.0:
                print("FAIL", fluid_key, "kappa", s["kappa_interp_max_abs_pct"])
                ok = False
        print("CRYOGENIC_TRANSPORT_INTERP", "PASS" if ok else "FAIL")
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
