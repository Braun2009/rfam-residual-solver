#!/usr/bin/env python3
"""Focused regression for the six-species R2 equilibrium screen."""

from __future__ import annotations

import argparse
import csv
import io
import subprocess
import sys
from pathlib import Path

try:
    import cantera as ct
except ModuleNotFoundError:
    ct = None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    args = parser.parse_args()
    if ct is None:
        print("R2_SIX_SPECIES_EQUILIBRIUM_TEST_SKIPPED cantera-unavailable")
        return 0

    probe = args.root.resolve() / "tools/r2_six_species_equilibrium_probe.py"
    completed = subprocess.run(
        [sys.executable, str(probe), "--root", str(args.root.resolve())],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        print("R2_SIX_SPECIES_EQUILIBRIUM_TEST_FAIL probe-returned-nonzero")
        print(completed.stdout)
        print(completed.stderr, file=sys.stderr)
        return 1

    lines = completed.stdout.splitlines()
    if len(lines) < 4 or lines[0] != "# rfam.r2-six-species-equilibrium-screen/v1":
        print("R2_SIX_SPECIES_EQUILIBRIUM_TEST_FAIL schema-marker")
        return 1
    rows = list(csv.DictReader(io.StringIO("\n".join(lines[1:-1]))))
    checks = 0

    def check(condition: bool, label: str) -> None:
        nonlocal checks
        checks += 1
        if not condition:
            raise ValueError(label)

    try:
        check(len(rows) == 6, "case-count")
        check(
            {row["partition"] for row in rows}
            == {"protocol-training", "protocol-holdout"},
            "partition-set",
        )
        gaps = [float(row["temperature_gap_k"]) for row in rows]
        check(min(gaps) > 300.0, "all-temperature-gaps-exceed-300K")
        check(max(gaps) > 1500.0, "maximum-temperature-gap-exceeds-1500K")
        check(
            max(
                float(row["reference_terminal_temperature_residual_k"])
                for row in rows
            )
            <= 1.0e-3,
            "detailed-reference-terminal-is-equilibrium",
        )
        check(
            max(
                float(row["maximum_major_species_mass_fraction_gap"])
                for row in rows
            )
            > 0.15,
            "major-species-structural-gap-visible",
        )
        check(
            all(
                row["fit_approved"] == "0"
                and row["physical_validation_approved"] == "0"
                and row["production_ready"] == "0"
                for row in rows
            ),
            "approvals-remain-false",
        )
        check(
            lines[-1].startswith("# R2_SIX_SPECIES_EQUILIBRIUM_BLOCKED cases=6 ")
            and lines[-1].endswith(
                "fit_approved=0 physical_validation_approved=0 production_ready=0"
            ),
            "blocked-summary",
        )
    except (KeyError, ValueError) as error:
        print(f"R2_SIX_SPECIES_EQUILIBRIUM_TEST_FAIL reason={error}")
        return 1

    print(f"R2_SIX_SPECIES_EQUILIBRIUM_TEST_PASS checks={checks} cases={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
