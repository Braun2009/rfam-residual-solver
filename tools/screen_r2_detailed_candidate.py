#!/usr/bin/env python3
"""Screen the current RFAM two-step trajectory against a Cantera reference.

This tool deliberately emits diagnostics rather than pass/fail physics. The
input detailed mechanism and its data rights must be reviewed separately before
any fit, publication, or production claim.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

from run_r2_0d_rfam_protocol import (
    ProtocolError,
    finite_number,
    load_json,
    parse_runner_output,
    sha256,
)


SPECIES = ["C3H8", "O2", "CO", "CO2", "H2O", "N2"]
EXPECTED_REFERENCE_SCHEMA = "rfam.r2-0d-cantera-reference/v1"
EXPECTED_REFERENCE_STATUS = "local-cantera-reference-not-fit-eligible"


def stable(value: float) -> float:
    if not math.isfinite(value):
        raise ProtocolError("nonfinite-screening-result")
    return float(format(value, ".12g"))


def require_hex_sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise ProtocolError(f"invalid-sha256:{label}")
    if any(character not in "0123456789abcdef" for character in value):
        raise ProtocolError(f"invalid-sha256:{label}")
    return value


def validate_reference(document: dict[str, Any]) -> tuple[list[float], list[dict[str, Any]]]:
    if document.get("schema") != EXPECTED_REFERENCE_SCHEMA:
        raise ProtocolError("reference-schema-mismatch")
    if document.get("status") != EXPECTED_REFERENCE_STATUS:
        raise ProtocolError("reference-status-mismatch")
    if document.get("major_species") != SPECIES:
        raise ProtocolError("reference-species-order-mismatch")
    require_hex_sha256(document.get("mechanism_sha256"), "mechanism")
    require_hex_sha256(document.get("generator_sha256"), "generator")

    sample_times_raw = document.get("sample_times_s")
    if not isinstance(sample_times_raw, list) or len(sample_times_raw) < 2:
        raise ProtocolError("reference-sample-times-invalid")
    sample_times = [
        finite_number(value, f"reference.sample_times[{index}]")
        for index, value in enumerate(sample_times_raw)
    ]
    if sample_times[0] != 0.0 or any(
        right <= left for left, right in zip(sample_times, sample_times[1:])
    ):
        raise ProtocolError("reference-sample-times-not-increasing")

    cases = document.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ProtocolError("reference-cases-invalid")
    case_ids: set[str] = set()
    for case_index, case in enumerate(cases):
        if not isinstance(case, dict):
            raise ProtocolError(f"reference-case-not-object:{case_index}")
        case_id = case.get("id")
        if not isinstance(case_id, str) or not case_id or case_id in case_ids:
            raise ProtocolError(f"reference-case-id-invalid:{case_index}")
        case_ids.add(case_id)
        finite_number(case.get("temperature_k"), f"{case_id}.temperature")
        finite_number(case.get("pressure_pa"), f"{case_id}.pressure")
        ignition_delay = finite_number(
            case.get("ignition_delay_s"), f"{case_id}.ignition-delay")
        if ignition_delay <= 0.0:
            raise ProtocolError(f"reference-ignition-delay-invalid:{case_id}")
        initial = case.get("initial_major_species_mass_fractions")
        if not isinstance(initial, dict) or set(initial) != set(SPECIES):
            raise ProtocolError(f"reference-initial-composition-invalid:{case_id}")
        initial_values = [
            finite_number(initial[name], f"{case_id}.initial.{name}")
            for name in SPECIES
        ]
        if any(value < 0.0 for value in initial_values) or abs(sum(initial_values) - 1.0) > 2.0e-9:
            raise ProtocolError(f"reference-initial-composition-sum:{case_id}")
        samples = case.get("samples")
        if not isinstance(samples, list) or len(samples) != len(sample_times):
            raise ProtocolError(f"reference-sample-count-mismatch:{case_id}")
        for sample_index, (sample, expected_time) in enumerate(zip(samples, sample_times)):
            if not isinstance(sample, dict):
                raise ProtocolError(f"reference-sample-not-object:{case_id}:{sample_index}")
            actual_time = finite_number(
                sample.get("time_s"), f"{case_id}.sample[{sample_index}].time")
            if actual_time != expected_time:
                raise ProtocolError(f"reference-sample-time-mismatch:{case_id}:{sample_index}")
            composition = sample.get("major_species_mass_fractions")
            if not isinstance(composition, dict) or set(composition) != set(SPECIES):
                raise ProtocolError(f"reference-sample-composition-invalid:{case_id}:{sample_index}")
            for name in SPECIES:
                finite_number(
                    composition[name],
                    f"{case_id}.sample[{sample_index}].{name}",
                )
            finite_number(
                sample.get("temperature_k"),
                f"{case_id}.sample[{sample_index}].temperature",
            )
            finite_number(
                sample.get("heat_release_rate_w_m3"),
                f"{case_id}.sample[{sample_index}].qdot",
            )
    return sample_times, cases


def run_rfam_case(
    runner: Path,
    case: dict[str, Any],
    sample_times: list[float],
) -> tuple[dict[str, str], list[dict[str, Any]]]:
    initial = case["initial_major_species_mass_fractions"]
    command = [
        str(runner),
        case["id"],
        format(float(case["temperature_k"]), ".17g"),
        format(float(case["pressure_pa"]), ".17g"),
        str(len(sample_times)),
        "0",
        *[format(float(initial[name]), ".17g") for name in SPECIES],
        *[format(value, ".17g") for value in sample_times],
    ]
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    if result.returncode != 0:
        raise ProtocolError(
            f"runner-rejected:{case['id']}:rc={result.returncode}:"
            f"{result.stdout.strip()}:{result.stderr.strip()}"
        )
    return parse_runner_output(
        result.stdout, result.stderr, case["id"], len(sample_times))


def compare_case(
    case: dict[str, Any],
    rfam_headers: dict[str, str],
    rfam_samples: list[dict[str, Any]],
) -> dict[str, Any]:
    detailed_samples = case["samples"]
    detailed_idt = finite_number(case["ignition_delay_s"], "detailed-idt")
    if rfam_headers.get("ignition_detected") != "1":
        raise ProtocolError(f"rfam-ignition-not-detected:{case['id']}")
    rfam_idt = finite_number(rfam_headers.get("ignition_delay_s"), "rfam-idt")
    if rfam_idt <= 0.0:
        raise ProtocolError(f"rfam-ignition-delay-invalid:{case['id']}")

    maximum_temperature_error = 0.0
    maximum_species_l1 = 0.0
    for detailed, rfam in zip(detailed_samples, rfam_samples):
        maximum_temperature_error = max(
            maximum_temperature_error,
            abs(
                finite_number(rfam["temperature_k"], "rfam-temperature")
                - finite_number(detailed["temperature_k"], "detailed-temperature")
            ),
        )
        detailed_y = detailed["major_species_mass_fractions"]
        species_l1 = sum(
            abs(
                finite_number(rfam["mass_fractions"][index], "rfam-Y")
                - finite_number(detailed_y[name], "detailed-Y")
            )
            for index, name in enumerate(SPECIES)
        )
        maximum_species_l1 = max(maximum_species_l1, species_l1)

    detailed_qdot0 = finite_number(
        detailed_samples[0]["heat_release_rate_w_m3"], "detailed-qdot0")
    rfam_qdot0 = finite_number(rfam_samples[0]["qdot_w_m3"], "rfam-qdot0")
    final_temperature_delta = finite_number(
        rfam_samples[-1]["temperature_k"], "rfam-final-temperature")
    final_temperature_delta -= finite_number(
        detailed_samples[-1]["temperature_k"], "detailed-final-temperature")

    return {
        "case_id": case["id"],
        "partition": case.get("partition"),
        "initial_temperature_k": stable(float(case["temperature_k"])),
        "pressure_pa": stable(float(case["pressure_pa"])),
        "detailed_ignition_delay_s": stable(detailed_idt),
        "rfam_ignition_delay_s": stable(rfam_idt),
        "rfam_to_detailed_ignition_delay_ratio": stable(rfam_idt / detailed_idt),
        "detailed_initial_qdot_w_m3": stable(detailed_qdot0),
        "rfam_initial_qdot_w_m3": stable(rfam_qdot0),
        "initial_qdot_absolute_ratio": stable(
            abs(rfam_qdot0) / max(abs(detailed_qdot0), 1.0)
        ),
        "initial_qdot_sign_matches": math.copysign(1.0, rfam_qdot0)
        == math.copysign(1.0, detailed_qdot0),
        "final_temperature_delta_k": stable(final_temperature_delta),
        "maximum_temperature_absolute_error_k": stable(maximum_temperature_error),
        "maximum_major_species_l1": stable(maximum_species_l1),
    }


def generate(reference_path: Path, runner: Path) -> dict[str, Any]:
    reference = load_json(reference_path, "detailed-reference")
    sample_times, cases = validate_reference(reference)
    if not runner.is_file():
        raise ProtocolError("rfam-runner-missing")

    comparisons = []
    for case in cases:
        headers, samples = run_rfam_case(runner, case, sample_times)
        comparisons.append(compare_case(case, headers, samples))

    ignition_ratios = [
        item["rfam_to_detailed_ignition_delay_ratio"] for item in comparisons
    ]
    temperature_deltas = [item["final_temperature_delta_k"] for item in comparisons]
    qdot_ratios = [item["initial_qdot_absolute_ratio"] for item in comparisons]
    return {
        "schema": "rfam.r2-detailed-candidate-screen/v1",
        "status": "internal-screening-only-no-fit-no-physical-validation",
        "source": {
            "detailed_reference_path": str(reference_path),
            "detailed_reference_sha256": sha256(reference_path),
            "detailed_mechanism_sha256": reference["mechanism_sha256"],
            "detailed_generator_sha256": reference["generator_sha256"],
            "cantera_version": reference.get("cantera_version"),
            "rfam_runner_path": str(runner),
            "rfam_runner_sha256": sha256(runner),
        },
        "scope": {
            "reactor": "ideal-gas-constant-pressure-high-temperature-screen",
            "sample_times_s": [stable(value) for value in sample_times],
            "species_order": SPECIES,
            "rights_clear_for_fit": False,
            "fit_executed": False,
            "acceptance_thresholds_defined": False,
        },
        "comparisons": comparisons,
        "aggregate": {
            "minimum_rfam_to_detailed_ignition_delay_ratio": stable(min(ignition_ratios)),
            "maximum_rfam_to_detailed_ignition_delay_ratio": stable(max(ignition_ratios)),
            "minimum_final_temperature_delta_k": stable(min(temperature_deltas)),
            "maximum_final_temperature_delta_k": stable(max(temperature_deltas)),
            "minimum_initial_qdot_absolute_ratio": stable(min(qdot_ratios)),
            "maximum_initial_qdot_absolute_ratio": stable(max(qdot_ratios)),
            "all_initial_qdot_signs_match": all(
                item["initial_qdot_sign_matches"] for item in comparisons
            ),
        },
        "approvals": {
            "mechanism_selected_for_fit": False,
            "experimental_data_selected": False,
            "r1_r2_fitting_approved": False,
            "physical_validation_approved": False,
            "production_ready": False,
        },
        "interpretation": [
            "The current candidate must not be treated as an R2-only calibration problem.",
            "R1/R2 rates and six-species endpoint thermochemistry require separate scoring.",
            "A rights-cleared detailed mechanism and disjoint experimental holdout remain required.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--detailed-reference", type=Path, required=True)
    parser.add_argument("--rfam-runner", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        document = generate(
            args.detailed_reference.resolve(), args.rfam_runner.resolve())
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(document, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
            + "\n",
            encoding="utf-8",
        )
    except (OSError, ProtocolError, ValueError) as error:
        print(f"R2_DETAILED_CANDIDATE_SCREEN_FAIL reason={error}", file=sys.stderr)
        return 2
    print(
        "R2_DETAILED_CANDIDATE_SCREEN_COMPLETE "
        "scope=internal-no-fit-no-physical-validation"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
