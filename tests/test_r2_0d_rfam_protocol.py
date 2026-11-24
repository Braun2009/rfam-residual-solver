#!/usr/bin/env python3
"""Focused regressions for the offline RFAM R2 0-D runner protocol."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path


class RegressionFailure(AssertionError):
    pass


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RegressionFailure(label)


def run_protocol(
    root: Path,
    runner: Path,
    output: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    command = [
        sys.executable,
        str(root / "tools/run_r2_0d_rfam_protocol.py"),
        "--root",
        str(root),
        "--contract",
        str(root / "tests/reference/r2_0d_calibration_contract_v1.json"),
        "--reference",
        str(root / "tests/reference/r2_0d_cantera_reference_v1.json"),
        "--runner",
        str(runner),
    ]
    if output is not None:
        command.extend(["--output", str(output)])
    return subprocess.run(command, check=False, text=True, capture_output=True)


def valid_runner_command(root: Path, runner: Path) -> list[str]:
    reference = json.loads(
        (root / "tests/reference/r2_0d_cantera_reference_v1.json").read_text(
            encoding="utf-8"
        )
    )
    case = reference["cases"][0]
    return scheduled_runner_command(
        runner, reference, case, reference["sample_times_s"], 0.0
    )


def scheduled_runner_command(
    runner: Path,
    reference: dict,
    case: dict,
    sample_times: list[float],
    max_internal_step: float,
    solver_tolerance_profile: str | None = None,
) -> list[str]:
    command = [
        str(runner),
        case["id"],
        format(case["temperature_k"], ".17g"),
        format(case["pressure_pa"], ".17g"),
        str(len(sample_times)),
        format(max_internal_step, ".17g"),
        *[
            format(case["initial_major_species_mass_fractions"][species], ".17g")
            for species in reference["major_species"]
        ],
        *[format(value, ".17g") for value in sample_times],
    ]
    if solver_tolerance_profile is not None:
        command.append(
            "--solver-tolerance-profile=" + solver_tolerance_profile)
    return command


def run_scheduled_case(
    runner: Path,
    reference: dict,
    case: dict,
    sample_times: list[float],
    max_internal_step: float,
    solver_tolerance_profile: str | None = None,
) -> tuple[dict[str, str], list[dict[str, str]]]:
    result = subprocess.run(
        scheduled_runner_command(
            runner,
            reference,
            case,
            sample_times,
            max_internal_step,
            solver_tolerance_profile,
        ),
        check=False,
        text=True,
        capture_output=True,
    )
    require(result.returncode == 0, "scheduled-runner-rejected")
    require(not result.stderr.strip(), "scheduled-runner-stderr")
    require(
        result.stdout.splitlines()[-1] == "R2_0D_RFAM_RUNNER_OK",
        "scheduled-runner-marker-missing",
    )
    headers: dict[str, str] = {}
    samples: list[dict[str, str]] = []
    for line in result.stdout.splitlines()[:-1]:
        if line.startswith("sample="):
            samples.append(dict(token.split("=", 1) for token in line.split("|")))
        else:
            key, value = line.split("=", 1)
            headers[key] = value
    require(len(samples) == len(sample_times), "scheduled-sample-count-mismatch")
    return headers, samples


def write_approval_tamper_wrapper(path: Path, runner: Path) -> None:
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import subprocess, sys\n"
        f"runner = {str(runner)!r}\n"
        "result = subprocess.run([runner, *sys.argv[1:]], "
        "text=True, capture_output=True)\n"
        "sys.stderr.write(result.stderr)\n"
        "sys.stdout.write(result.stdout.replace("
        "'physical_validation_approved=0', "
        "'physical_validation_approved=1'))\n"
        "raise SystemExit(result.returncode)\n",
        encoding="ascii",
    )
    path.chmod(path.stat().st_mode | 0o111)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    args = parser.parse_args(argv)
    root = args.root.resolve()
    runner = args.runner.resolve()
    checks = 0
    try:
        require(runner.is_file(), "runner-missing")
        checks += 1
        with tempfile.TemporaryDirectory(prefix="rfam-r2-0d-") as directory:
            temporary = Path(directory)
            output = temporary / "trajectory.json"
            canonical = run_protocol(root, runner, output)
            require(canonical.returncode == 0, "canonical-protocol-rejected")
            require(
                "R2_0D_RFAM_PROTOCOL_BLOCKED cases=6 samples=84 checks=816"
                in canonical.stdout,
                "canonical-blocked-marker-missing",
            )
            require("R2_0D_RFAM_PROTOCOL_FAIL" not in canonical.stdout + canonical.stderr,
                    "canonical-fail-marker-present")
            require(output.is_file(), "canonical-output-missing")
            checks += 4

            document = json.loads(output.read_text(encoding="utf-8"))
            require(
                document["schema"] == "rfam.r2-0d-rfam-protocol-output/v1",
                "output-schema-mismatch",
            )
            require(
                document["status"]
                == "implementation-protocol-only-no-fit-or-physical-validation",
                "output-status-mismatch",
            )
            require(len(document["cases"]) == 6, "output-case-count-mismatch")
            require(
                sum(len(case["samples"]) for case in document["cases"]) == 84,
                "output-sample-count-mismatch",
            )
            require(
                all(value is False for value in document["approvals"].values()),
                "output-approval-enabled",
            )
            require(
                not any(
                    document["comparison"][key]
                    for key in (
                        "cantera_numeric_fit_score_generated",
                        "ignition_delay_fit_metric_enabled",
                        "trajectory_fit_metrics_enabled",
                        "heat_release_fit_metric_enabled",
                    )
                ),
                "output-fit-metric-enabled",
            )
            require(
                all(
                    int(case["runner_headers"]["depletion_event_count"]) >= 0
                    for case in document["cases"]
                ),
                "output-depletion-count-invalid",
            )
            checks += 7

            reference = json.loads(
                (root / "tests/reference/r2_0d_cantera_reference_v1.json")
                .read_text(encoding="utf-8")
            )
            case = reference["cases"][0]
            canonical_times = reference["sample_times_s"]
            sparse_times = [0.0, canonical_times[-1]]
            refined_times = sorted(
                set(
                    canonical_times
                    + [
                        2.0e-9,
                        4.0e-9,
                        6.0e-9,
                        8.0e-9,
                        2.0e-8,
                        5.0e-8,
                        2.0e-7,
                        5.0e-7,
                        2.0e-6,
                        5.0e-6,
                        2.0e-5,
                        5.0e-5,
                    ]
                )
            )
            grid_runs = [
                run_scheduled_case(runner, reference, case, times, 0.0)
                for times in (canonical_times, sparse_times, refined_times)
            ]
            for times, (headers, samples) in zip(
                (canonical_times, sparse_times, refined_times), grid_runs
            ):
                actual_times = [float(sample["time_s"]) for sample in samples]
                require(
                    all(
                        later > earlier
                        for earlier, later in zip(actual_times, actual_times[1:])
                    ),
                    "scheduled-time-not-monotonic",
                )
                require(
                    actual_times == times,
                    "scheduled-time-value-mismatch",
                )
                require(
                    float(headers["accepted_physical_time_s"]) == times[-1],
                    "accepted-physical-time-mismatch",
                )
                require(
                    float(headers["implicit_time_extension_s"]) == 0.0,
                    "implicit-time-extension-present",
                )
                require(
                    int(headers["depletion_event_count"]) == 1,
                    "depletion-event-count-mismatch",
                )
                checks += 5

            canonical_final = grid_runs[0][1][-1]
            for _, samples in grid_runs[1:]:
                final = samples[-1]
                require(
                    abs(
                        float(final["temperature_k"])
                        - float(canonical_final["temperature_k"])
                    )
                    <= 1.0e-8,
                    "output-grid-temperature-not-converged",
                )
                canonical_y = [
                    float(value)
                    for value in canonical_final["mass_fractions"].split(",")
                ]
                final_y = [
                    float(value) for value in final["mass_fractions"].split(",")
                ]
                require(
                    max(abs(actual - expected) for actual, expected in zip(final_y, canonical_y))
                    <= 1.0e-13,
                    "output-grid-composition-not-converged",
                )
                checks += 2

            short_times = [0.0, 1.0e-8, 3.0e-8, 1.0e-7]
            step_runs = [
                run_scheduled_case(runner, reference, case, short_times, max_step)
                for max_step in (0.0, 1.0e-9, 5.0e-10)
            ]
            adaptive_final = step_runs[0][1][-1]
            capped_finals = [run[1][-1] for run in step_runs[1:]]
            for final in capped_finals:
                require(
                    abs(
                        float(final["temperature_k"])
                        - float(adaptive_final["temperature_k"])
                    )
                    <= 1.0e-6,
                    "internal-step-temperature-not-converged",
                )
                adaptive_y = [
                    float(value)
                    for value in adaptive_final["mass_fractions"].split(",")
                ]
                capped_y = [
                    float(value) for value in final["mass_fractions"].split(",")
                ]
                require(
                    max(abs(actual - expected) for actual, expected in zip(capped_y, adaptive_y))
                    <= 1.0e-10,
                    "internal-step-composition-not-converged",
                )
                checks += 2

            strict_headers, strict_samples = run_scheduled_case(
                runner,
                reference,
                case,
                [0.0, 1.0e-14],
                0.0,
                "strict-equivalence-v1",
            )
            require(
                strict_headers["solver_tolerance_profile"]
                == "strict-equivalence-v1",
                "strict-profile-identity-mismatch",
            )
            require(
                float(strict_headers["solver_relative_tolerance"]) == 1.0e-12
                and float(strict_headers["solver_species_absolute_tolerance"])
                == 1.0e-16
                and float(strict_headers["solver_enthalpy_absolute_tolerance"])
                == 1.0e-8,
                "strict-profile-tolerance-mismatch",
            )
            for sample in strict_samples:
                require(
                    "temperature_rate_k_s" in sample
                    and "mass_fraction_rates_s_inv" in sample,
                    "strict-derivative-field-missing",
                )
                temperature_rate = float(sample["temperature_rate_k_s"])
                mass_fraction_rates = [
                    float(value)
                    for value in sample["mass_fraction_rates_s_inv"].split(",")
                ]
                require(
                    math.isfinite(temperature_rate)
                    and len(mass_fraction_rates) == 6
                    and all(math.isfinite(value) for value in mass_fraction_rates),
                    "strict-derivative-value-invalid",
                )
                require(
                    abs(sum(mass_fraction_rates)) <= 1.0e-8,
                    "strict-mass-fraction-rate-not-closed",
                )
            require(
                "temperature_rate_k_s" not in grid_runs[0][1][0]
                and "mass_fraction_rates_s_inv" not in grid_runs[0][1][0],
                "default-protocol-derivative-field-present",
            )
            checks += 8

            invalid_profile = scheduled_runner_command(
                runner,
                reference,
                case,
                [0.0, 1.0e-14],
                0.0,
                "strict-equivalence-v2",
            )
            rejected_profile = subprocess.run(
                invalid_profile, check=False, text=True, capture_output=True)
            require(rejected_profile.returncode == 2, "unknown-profile-accepted")
            require(
                "argument-count-mismatch" in rejected_profile.stderr,
                "unknown-profile-reason-missing",
            )
            require(
                "R2_0D_RFAM_RUNNER_OK"
                not in rejected_profile.stdout + rejected_profile.stderr,
                "unknown-profile-success-marker-present",
            )
            checks += 3

            missing_arguments = subprocess.run(
                [str(runner)], check=False, text=True, capture_output=True
            )
            require(missing_arguments.returncode == 2, "missing-arguments-accepted")
            require(
                "R2_0D_RFAM_RUNNER_FAIL arguments-required"
                in missing_arguments.stderr,
                "missing-arguments-marker-missing",
            )
            require(
                "R2_0D_RFAM_RUNNER_OK"
                not in missing_arguments.stdout + missing_arguments.stderr,
                "missing-arguments-success-marker-present",
            )
            checks += 3

            invalid_composition = valid_runner_command(root, runner)
            invalid_composition[6] = "0.5"
            rejected_composition = subprocess.run(
                invalid_composition, check=False, text=True, capture_output=True
            )
            require(rejected_composition.returncode == 2, "invalid-composition-accepted")
            require(
                "mass-fraction-sum-mismatch" in rejected_composition.stderr,
                "invalid-composition-reason-missing",
            )
            require(
                "R2_0D_RFAM_RUNNER_OK"
                not in rejected_composition.stdout + rejected_composition.stderr,
                "invalid-composition-success-marker-present",
            )
            checks += 3

            wrapper = temporary / "approval_tamper.py"
            write_approval_tamper_wrapper(wrapper, runner)
            tampered = run_protocol(root, wrapper)
            require(tampered.returncode == 2, "approval-tamper-accepted")
            require(
                "runner-header-value-mismatch:train_phi0p8_diluted:"
                "physical_validation_approved"
                in tampered.stderr,
                "approval-tamper-reason-missing",
            )
            require(
                "R2_0D_RFAM_PROTOCOL_BLOCKED"
                not in tampered.stdout + tampered.stderr,
                "approval-tamper-blocked-marker-present",
            )
            require("Traceback" not in tampered.stderr, "approval-tamper-traceback")
            checks += 4
    except (RegressionFailure, OSError, ValueError) as error:
        print(f"R2_0D_RFAM_PROTOCOL_TEST_FAIL {error}", file=sys.stderr)
        return 1
    except Exception as error:
        print(
            "R2_0D_RFAM_PROTOCOL_TEST_FAIL unexpected:"
            + type(error).__name__,
            file=sys.stderr,
        )
        return 1

    print(f"R2_0D_RFAM_PROTOCOL_TEST_PASS checks={checks}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
