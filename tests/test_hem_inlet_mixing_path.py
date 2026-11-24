#!/usr/bin/env python3
"""Validate the bounded RFAM LOX/propane inlet-mixing diagnostic."""

from __future__ import annotations

import argparse
import csv
import io
import math
import subprocess
import sys
from pathlib import Path


EXPECTED_COLUMNS = [
    "pressure_pa",
    "kij",
    "fuel_mass_fraction",
    "o2_mole_fraction",
    "target_h_j_kg",
    "status",
    "phase_status",
    "temperature_k",
    "region",
    "vapor_mass_fraction",
    "enthalpy_relative_residual",
    "valid_temperature_samples",
    "invalid_temperature_samples",
    "phase_topology_validated",
    "production_ready",
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    args = parser.parse_args()

    probe = str(args.probe.resolve())
    completed = subprocess.run(
        [probe],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        print("HEM_INLET_MIXING_PATH_TEST_FAIL probe-returned-nonzero")
        print(completed.stdout)
        print(completed.stderr, file=sys.stderr)
        return 1

    refined = subprocess.run(
        [probe, "--scan-intervals", "1024"],
        check=False,
        capture_output=True,
        text=True,
    )
    if refined.returncode != 0:
        print("HEM_INLET_MIXING_PATH_TEST_FAIL refined-probe-returned-nonzero")
        print(refined.stdout)
        print(refined.stderr, file=sys.stderr)
        return 1

    finest = subprocess.run(
        [probe, "--scan-intervals", "2048"],
        check=False,
        capture_output=True,
        text=True,
    )
    if finest.returncode != 0:
        print("HEM_INLET_MIXING_PATH_TEST_FAIL finest-probe-returned-nonzero")
        print(finest.stdout)
        print(finest.stderr, file=sys.stderr)
        return 1

    lines = completed.stdout.splitlines()
    if not lines or lines[0] != "# rfam.hem-inlet-mixing-path-diagnostic/v1":
        print("HEM_INLET_MIXING_PATH_TEST_FAIL schema-marker")
        return 1
    expected_summary = (
        "# HEM_INLET_MIXING_PATH_DIAGNOSTIC_BLOCKED rows=225 "
        "successes=220 failures=5 unbracketed=3 two_phase=126 "
        "scan_intervals=512 "
        "phase_topology_validated=0 production_ready=0"
    )
    if lines[-1] != expected_summary:
        print("HEM_INLET_MIXING_PATH_TEST_FAIL blocked-marker")
        return 1

    rows = list(csv.DictReader(io.StringIO("\n".join(lines[1:-1]))))
    refined_lines = refined.stdout.splitlines()
    refined_summary = expected_summary.replace(
        "scan_intervals=512", "scan_intervals=1024"
    )
    if not refined_lines or refined_lines[-1] != refined_summary:
        print("HEM_INLET_MIXING_PATH_TEST_FAIL refined-blocked-marker")
        return 1
    refined_rows = list(
        csv.DictReader(io.StringIO("\n".join(refined_lines[1:-1])))
    )
    finest_lines = finest.stdout.splitlines()
    finest_summary = expected_summary.replace(
        "scan_intervals=512", "scan_intervals=2048"
    )
    if not finest_lines or finest_lines[-1] != finest_summary:
        print("HEM_INLET_MIXING_PATH_TEST_FAIL finest-blocked-marker")
        return 1
    finest_rows = list(
        csv.DictReader(io.StringIO("\n".join(finest_lines[1:-1])))
    )
    checks = 0

    def check(condition: bool, label: str) -> None:
        nonlocal checks
        checks += 1
        if not condition:
            raise ValueError(label)

    try:
        check(list(rows[0]) == EXPECTED_COLUMNS, "exact-columns")
        check(len(rows) == 225, "row-count")
        check(
            {float(row["pressure_pa"]) for row in rows}
            == {4.0e6, 5.0e6, 6.0e6},
            "pressure-grid",
        )
        check(
            {float(row["kij"]) for row in rows}
            == {-0.05, 0.0, 0.055, 0.10, 0.15},
            "kij-grid",
        )
        check(
            all(row["phase_topology_validated"] == "0" for row in rows),
            "topology-remains-unvalidated",
        )
        check(
            all(row["production_ready"] == "0" for row in rows),
            "production-remains-disabled",
        )
        check(len(refined_rows) == len(rows), "refined-row-count")
        check(len(finest_rows) == len(rows), "finest-row-count")
        for base, fine, finest_row in zip(rows, refined_rows, finest_rows):
            check(
                base["pressure_pa"] == fine["pressure_pa"]
                == finest_row["pressure_pa"]
                and base["kij"] == fine["kij"]
                == finest_row["kij"]
                and base["fuel_mass_fraction"] == fine["fuel_mass_fraction"]
                == finest_row["fuel_mass_fraction"]
                and base["status"] == fine["status"]
                == finest_row["status"]
                and base["phase_status"] == fine["phase_status"]
                == finest_row["phase_status"]
                and base["region"] == fine["region"]
                == finest_row["region"],
                "512-to-1024-to-2048-state-classification",
            )
            check(
                base["temperature_k"] == fine["temperature_k"]
                == finest_row["temperature_k"]
                and base["vapor_mass_fraction"]
                == fine["vapor_mass_fraction"]
                == finest_row["vapor_mass_fraction"]
                and base["enthalpy_relative_residual"]
                == fine["enthalpy_relative_residual"]
                == finest_row["enthalpy_relative_residual"],
                "512-to-1024-to-2048-result",
            )
        check(
            any(row["status"] != "SUCCESS" for row in rows),
            "failure-region-is-visible",
        )
        check(
            any(row["region"] == "TWO_PHASE" for row in rows),
            "two-phase-region-is-visible",
        )
        global_fraction = 1.05 / 3.85
        global_rows = [
            row
            for row in rows
            if math.isclose(
                float(row["fuel_mass_fraction"]),
                global_fraction,
                rel_tol=0.0,
                abs_tol=1.0e-15,
            )
        ]
        check(len(global_rows) == 15, "global-mixture-at-all-pressure-kij-pairs")
        check(
            all(
                math.isfinite(float(row["o2_mole_fraction"]))
                and 0.0 < float(row["o2_mole_fraction"]) < 1.0
                for row in global_rows
            ),
            "global-mole-basis-is-finite",
        )
        global_6mpa = {
            float(row["kij"]): row["status"]
            for row in global_rows
            if float(row["pressure_pa"]) == 6.0e6
        }
        check(
            global_6mpa
            == {
                -0.05: "SUCCESS",
                0.0: "SUCCESS",
                0.055: "SUCCESS",
                0.10: "SUCCESS",
                0.15: "SUCCESS",
            },
            "global-6mpa-kij-dependent-ph-closure",
        )
        check(
            all(
                row["status"] == "SUCCESS"
                and row["region"] == "TWO_PHASE"
                for row in global_rows
                if float(row["pressure_pa"]) in {4.0e6, 5.0e6}
            ),
            "global-4mpa-5mpa-two-phase-diagnostic",
        )
        check(
            "phase_topology_validated=0 production_ready=0" in lines[-1],
            "summary-retains-blocked-scope",
        )
    except (KeyError, ValueError) as exc:
        print(f"HEM_INLET_MIXING_PATH_TEST_FAIL reason={exc}")
        return 1

    print(
        f"HEM_INLET_MIXING_PATH_TEST_PASS checks={checks} "
        f"rows={len(rows)} refined_rows={len(refined_rows)} "
        f"finest_rows={len(finest_rows)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
