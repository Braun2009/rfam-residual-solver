#!/usr/bin/env python3
"""Validate the bounded PR binary phase-topology risk screen."""

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
    "temperature_k",
    "kij",
    "composition_intervals",
    "valid_samples",
    "invalid_samples",
    "hull_vertices",
    "tie_segments",
    "lle_segments",
    "vle_segments",
    "other_segments",
    "maximum_contact_clusters",
    "three_phase_tangent_pairs",
    "three_phase_tie_segments",
    "maximum_three_phase_shared_phase_error",
    "maximum_three_phase_slope_error",
    "maximum_three_phase_intercept_error",
    "maximum_skipped_samples",
    "refined_tie_segments",
    "unrefined_tie_segments",
    "tp_supported_tie_segments",
    "tp_disagreements",
    "tp_flash_failures",
    "tp_non_two_phase",
    "tp_different_splits",
    "tp_failure_status",
    "tp_stability_attempted",
    "tp_stability_converged",
    "tp_stability_failed",
    "tp_flash_iterations",
    "tp_maximum_fugacity_residual",
    "tp_maximum_material_residual",
    "tp_liquid_minimum_tpd",
    "tp_vapor_minimum_tpd",
    "maximum_tp_endpoint_error",
    "widest_tie_span",
    "widest_tie_left_x_o2",
    "widest_tie_right_x_o2",
    "widest_tie_left_z",
    "widest_tie_right_z",
    "widest_tie_fugacity_residual",
    "tie_left_stability_status",
    "tie_left_stability_attempted",
    "tie_left_stability_converged",
    "tie_left_stability_failed",
    "tie_left_minimum_tpd",
    "tie_right_stability_status",
    "tie_right_stability_attempted",
    "tie_right_stability_converged",
    "tie_right_stability_failed",
    "tie_right_minimum_tpd",
    "possible_three_phase",
    "production_ready",
]


def run(probe: Path, intervals: int) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(probe.resolve()), "--composition-intervals", str(intervals)],
        capture_output=True,
        text=True,
        check=False,
    )


def run_single(
    probe: Path,
    pressure: float,
    temperature: float,
    kij: float,
    intervals: int,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            str(probe.resolve()),
            "--state",
            repr(pressure),
            repr(temperature),
            repr(kij),
            "--composition-intervals",
            str(intervals),
        ],
        capture_output=True,
        text=True,
        check=False,
    )


def parse(result: subprocess.CompletedProcess[str], intervals: int) -> tuple[list[dict[str, str]], str]:
    if result.returncode != 0:
        raise ValueError(
            f"probe-{intervals}-returned-{result.returncode}:"
            f"{result.stdout}{result.stderr}"
        )
    lines = result.stdout.splitlines()
    if not lines or lines[0] != "# rfam.pr-binary-phase-topology-screen/v1":
        raise ValueError(f"schema-marker-{intervals}")
    if len(lines) < 4:
        raise ValueError(f"short-output-{intervals}")
    rows = list(csv.DictReader(io.StringIO("\n".join(lines[1:-1]))))
    if list(rows[0]) != EXPECTED_COLUMNS:
        raise ValueError(f"exact-columns-{intervals}")
    return rows, lines[-1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    args = parser.parse_args()

    checks = 0

    def check(condition: bool, label: str) -> None:
        nonlocal checks
        checks += 1
        if not condition:
            raise ValueError(label)

    try:
        coarse_rows, coarse_summary = parse(run(args.probe, 1000), 1000)
        fine_rows, fine_summary = parse(run(args.probe, 2000), 2000)
        check(len(coarse_rows) == 270, "coarse-row-count")
        check(len(fine_rows) == 270, "fine-row-count")
        check(
            {float(row["pressure_pa"]) for row in fine_rows}
            == {4.0e6, 5.0e6, 6.0e6},
            "pressure-grid",
        )
        check(
            {float(row["kij"]) for row in fine_rows}
            == {-0.05, 0.0, 0.055, 0.10, 0.15},
            "kij-grid",
        )
        check(
            {float(row["temperature_k"]) for row in fine_rows}
            == {
                110.0, 120.0, 130.0, 140.0, 150.0, 160.0,
                170.0, 180.0, 190.0, 200.0, 210.0, 220.0,
                230.0, 240.0, 250.0, 260.0, 280.0, 300.0,
            },
            "temperature-grid",
        )
        check(
            all(row["production_ready"] == "0" for row in fine_rows),
            "production-remains-disabled",
        )
        check(
            all(int(row["valid_samples"]) == 2001 for row in fine_rows),
            "all-fine-samples-valid",
        )
        check(
            all(int(row["invalid_samples"]) == 0 for row in fine_rows),
            "no-invalid-fine-samples",
        )
        check(
            any(int(row["tie_segments"]) > 0 for row in fine_rows),
            "tie-risk-is-visible",
        )
        check(
            all(
                int(row["tie_segments"])
                == int(row["lle_segments"])
                + int(row["vle_segments"])
                + int(row["other_segments"])
                for row in fine_rows
            ),
            "tie-classification-closes",
        )
        check(
            all(
                int(row["tie_segments"])
                == int(row["tp_supported_tie_segments"])
                + int(row["tp_disagreements"])
                + int(row["three_phase_tie_segments"])
                for row in fine_rows
            ),
            "tp-cross-check-closes",
        )
        check(
            all(
                int(row["tie_segments"])
                == int(row["refined_tie_segments"])
                + int(row["unrefined_tie_segments"])
                for row in fine_rows
            ),
            "tie-refinement-accounting-closes",
        )
        check(
            any(int(row["refined_tie_segments"]) > 0 for row in fine_rows),
            "refined-tie-evidence-is-visible",
        )
        check(
            all(
                int(row["tp_supported_tie_segments"])
                == int(row["tie_segments"])
                - int(row["three_phase_tie_segments"])
                for row in fine_rows
            ),
            "all-nonthree-phase-fine-ties-are-tp-supported",
        )
        check(
            all(int(row["tp_disagreements"]) == 0 for row in fine_rows),
            "zero-fine-tp-disagreements",
        )
        check(
            all(
                int(row["tp_flash_failures"]) == 0
                and int(row["tp_non_two_phase"]) == 0
                and int(row["tp_different_splits"]) == 0
                and row["tp_failure_status"] == "NONE"
                for row in fine_rows
            ),
            "zero-fine-tp-failure-classifications",
        )
        check(
            all(
                (row["possible_three_phase"] == "1")
                == (int(row["three_phase_tie_segments"]) > 0)
                for row in fine_rows
            ),
            "three-phase-flag-is-attributable",
        )
        check(
            all(
                int(row["tp_disagreements"])
                == int(row["tp_flash_failures"])
                + int(row["tp_non_two_phase"])
                + int(row["tp_different_splits"])
                for row in fine_rows
            ),
            "tp-disagreement-classification-closes",
        )
        check(
            all(
                (int(row["tp_flash_failures"]) == 0
                 and row["tp_failure_status"] == "NONE")
                or (int(row["tp_flash_failures"]) > 0
                    and row["tp_failure_status"] != "NONE")
                for row in fine_rows
            ),
            "tp-failure-status-is-attributable",
        )
        check(
            all(
                (int(row["tp_flash_failures"]) == 0
                 and int(row["tp_stability_attempted"]) == 0
                 and int(row["tp_stability_converged"]) == 0
                 and int(row["tp_stability_failed"]) == 0)
                or (int(row["tp_flash_failures"]) > 0
                    and int(row["tp_stability_attempted"])
                    == int(row["tp_stability_converged"])
                    + int(row["tp_stability_failed"]))
                for row in fine_rows
            ),
            "tp-stability-trial-accounting-closes",
        )
        check(
            all(
                int(row["tp_flash_iterations"]) >= 0
                for row in fine_rows
            ),
            "tp-flash-iteration-count-is-nonnegative",
        )
        check(
            "phase_topology_validated=0 production_ready=0" in fine_summary,
            "fine-summary-retains-blocked-scope",
        )
        check(
            "phase_topology_validated=0 production_ready=0" in coarse_summary,
            "coarse-summary-retains-blocked-scope",
        )

        coarse_by_state = {
            (row["pressure_pa"], row["temperature_k"], row["kij"]): row
            for row in coarse_rows
        }
        fine_by_state = {
            (row["pressure_pa"], row["temperature_k"], row["kij"]): row
            for row in fine_rows
        }
        check(coarse_by_state.keys() == fine_by_state.keys(), "same-state-grid")
        single_rows, single_summary = parse(
            run_single(args.probe, 4.0e6, 200.0, 0.055, 2000), 2000
        )
        check(len(single_rows) == 1, "single-state-row-count")
        check(
            single_rows[0] == fine_by_state[("4000000", "200", "0.055")],
            "single-state-matches-matrix",
        )
        check(
            "rows=1" in single_summary
            and "phase_topology_validated=0 production_ready=0" in single_summary,
            "single-state-summary-retains-blocked-scope",
        )
        triple_states = (
            (5.0e6, 156.37717703, 0.055),
            (5.0e6, 155.6846784, 0.10),
            (5.0e6, 155.1625280618, 0.15),
        )
        for pressure, temperature, kij in triple_states:
            triple_rows, triple_summary = parse(
                run_single(args.probe, pressure, temperature, kij, 20000),
                20000,
            )
            check(len(triple_rows) == 1, f"three-phase-row-count:{kij}")
            triple = triple_rows[0]
            check(
                int(triple["tie_segments"]) == 2
                and int(triple["maximum_contact_clusters"]) == 2,
                f"adjacent-two-tie-topology:{kij}",
            )
            check(
                int(triple["three_phase_tangent_pairs"]) == 1
                and int(triple["three_phase_tie_segments"]) == 2
                and triple["possible_three_phase"] == "1",
                f"shared-middle-phase-detected:{kij}",
            )
            check(
                float(triple["maximum_three_phase_shared_phase_error"])
                <= 1.0e-7
                and float(triple["maximum_three_phase_slope_error"])
                <= 1.0e-8
                and float(triple["maximum_three_phase_intercept_error"])
                <= 1.0e-8,
                f"three-phase-common-tangent-closure:{kij}",
            )
            check(
                int(triple["tp_supported_tie_segments"]) == 0
                and int(triple["tp_disagreements"]) == 0
                and triple["tp_failure_status"] == "NONE",
                f"three-phase-ties-skip-two-phase-admission:{kij}",
            )
            check(
                "possible_three_phase_rows=1" in triple_summary
                and "phase_topology_validated=0 production_ready=0"
                in triple_summary,
                f"three-phase-summary-blocked:{kij}",
            )
        malformed = subprocess.run(
            [
                str(args.probe.resolve()),
                "--state",
                "4000000",
                "200",
                "0.055junk",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        check(malformed.returncode == 2, "malformed-state-return-code")
        check(
            "PR_BINARY_PHASE_TOPOLOGY_SCREEN_FAIL reason=invalid-state"
            in malformed.stderr,
            "malformed-state-marker",
        )
        check(
            "rfam.pr-binary-phase-topology-screen" not in malformed.stdout,
            "malformed-state-no-output",
        )
        robust_tie_rows = 0
        for identity, coarse in coarse_by_state.items():
            fine = fine_by_state[identity]
            coarse_has_tie = int(coarse["tie_segments"]) > 0
            fine_has_tie = int(fine["tie_segments"]) > 0
            check(coarse_has_tie == fine_has_tie, f"tie-classification:{identity}")
            if fine_has_tie:
                robust_tie_rows += 1
                coarse_span = float(coarse["widest_tie_span"])
                fine_span = float(fine["widest_tie_span"])
                check(
                    int(fine["refined_tie_segments"])
                    == int(fine["tie_segments"])
                    and int(fine["unrefined_tie_segments"]) == 0,
                    f"all-ties-refined:{identity}",
                )
                if int(fine["three_phase_tie_segments"]) == 0:
                    check(
                        float(fine["maximum_tp_endpoint_error"]) <= 1.0e-7,
                        f"tp-endpoint-agreement:{identity}",
                    )
                    check(
                        float(fine["tp_maximum_fugacity_residual"]) <= 1.0e-9
                        and float(fine["tp_maximum_material_residual"])
                        <= 1.0e-10,
                        f"tp-residual-closure:{identity}",
                    )
                    check(
                        float(fine["tp_liquid_minimum_tpd"]) >= -5.0e-8
                        and float(fine["tp_vapor_minimum_tpd"]) >= -5.0e-8,
                        f"tp-post-split-stability:{identity}",
                    )
                else:
                    check(
                        int(fine["three_phase_tangent_pairs"]) > 0
                        or int(fine["maximum_contact_clusters"]) >= 3,
                        f"three-phase-evidence-visible:{identity}",
                    )
                check(
                    math.isclose(coarse_span, fine_span, rel_tol=0.0, abs_tol=0.004),
                    f"tie-span-grid-convergence:{identity}",
                )
                check(
                    0.0 <= float(fine["widest_tie_left_x_o2"])
                    < float(fine["widest_tie_right_x_o2"]) <= 1.0,
                    f"tie-endpoint-composition:{identity}",
                )
                check(
                    float(fine["widest_tie_left_z"]) > 0.0
                    and float(fine["widest_tie_right_z"]) > 0.0,
                    f"tie-endpoint-z:{identity}",
                )
                check(
                    math.isfinite(float(fine["widest_tie_fugacity_residual"])),
                    f"tie-fugacity-residual-finite:{identity}",
                )
                check(
                    fine["tie_left_stability_status"] == "SUCCESS"
                    and fine["tie_right_stability_status"] == "SUCCESS",
                    f"tie-phase-stability-status:{identity}",
                )
                for side in ("left", "right"):
                    attempted = int(fine[f"tie_{side}_stability_attempted"])
                    converged = int(fine[f"tie_{side}_stability_converged"])
                    failed = int(fine[f"tie_{side}_stability_failed"])
                    check(
                        attempted == 6 and converged == 6 and failed == 0,
                        f"tie-{side}-stability-accounting:{identity}",
                    )
                    check(
                        float(fine[f"tie_{side}_minimum_tpd"]) >= -1.0e-8,
                        f"tie-{side}-minimum-tpd:{identity}",
                    )
        check(robust_tie_rows > 0, "robust-tie-row-count")
    except (KeyError, ValueError) as error:
        print(f"PR_BINARY_PHASE_TOPOLOGY_TEST_FAIL reason={error}")
        return 1

    print(
        f"PR_BINARY_PHASE_TOPOLOGY_TEST_PASS checks={checks} "
        f"rows={len(fine_rows)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
