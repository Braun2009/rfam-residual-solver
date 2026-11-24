#!/usr/bin/env python3
"""Run the RFAM TwoStepCO 0-D implementation protocol on six frozen cases.

The Cantera document supplies only the predeclared initial conditions, case
partitions, and output times.  This tool deliberately does not score the RFAM
trajectory against the local GRI30 trajectory: a detailed mechanism, an
experimental training/holdout set, and a common thermodynamic basis are still
unselected.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any


SPECIES = ["C3H8", "O2", "CO", "CO2", "H2O", "N2"]
MOLECULAR_WEIGHTS = [
    0.04409562,
    0.0319988,
    0.0280101,
    0.0440095,
    0.01801528,
    0.0280134,
]
ELEMENT_ATOMS = [
    [3, 8, 0, 0],
    [0, 0, 2, 0],
    [1, 0, 1, 0],
    [1, 0, 2, 0],
    [0, 2, 1, 0],
    [0, 0, 0, 2],
]
GAS_CONSTANT = 8.314462618
HEADER_KEYS = {
    "schema",
    "status",
    "case_id",
    "reactor",
    "chemistry",
    "energy_basis",
    "max_internal_step_s",
    "species_order",
    "sample_count",
    "ignition_definition",
    "ignition_detected",
    "ignition_delay_s",
    "positivity_strategy",
    "depletion_event_count",
    "interval_bisection_count",
    "roundoff_clamp_count",
    "accepted_physical_time_s",
    "implicit_time_extension_s",
    "fit_approved",
    "physical_validation_approved",
    "production_ready",
}
SAMPLE_KEYS = {
    "sample",
    "time_s",
    "temperature_k",
    "pressure_pa",
    "density_kg_m3",
    "sensible_enthalpy_j_kg",
    "total_enthalpy_j_kg",
    "qdot_w_m3",
    "mass_fractions",
}


class ProtocolError(ValueError):
    pass


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise ProtocolError(f"artifact-open-failed:{path}:{error}") from error


def reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ProtocolError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        document = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=reject_duplicate_pairs,
        )
    except OSError as error:
        raise ProtocolError(f"{label}-open-failed:{error}") from error
    except json.JSONDecodeError as error:
        raise ProtocolError(f"{label}-json-invalid:{error.msg}") from error
    if not isinstance(document, dict):
        raise ProtocolError(f"{label}-not-object")
    return document


def finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ProtocolError(f"invalid-number:{label}")
    try:
        result = float(value)
    except ValueError as error:
        raise ProtocolError(f"invalid-number:{label}") from error
    if not math.isfinite(result):
        raise ProtocolError(f"nonfinite-number:{label}")
    return result


def integer_text(value: str, label: str) -> int:
    if not value or any(character < "0" or character > "9" for character in value):
        raise ProtocolError(f"invalid-integer:{label}")
    return int(value)


def close(actual: float, expected: float, relative: float, absolute: float) -> bool:
    return abs(actual - expected) <= absolute + relative * abs(expected)


def require_close(
    actual: float,
    expected: float,
    relative: float,
    absolute: float,
    label: str,
) -> None:
    if not close(actual, expected, relative, absolute):
        raise ProtocolError(
            f"value-mismatch:{label}:actual={actual:.17g}:expected={expected:.17g}"
        )


def validate_upstream(
    root: Path,
    contract_path: Path,
    reference_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    validator = root / "tools/validate_r2_0d_calibration.py"
    result = subprocess.run(
        [
            sys.executable,
            str(validator),
            "--root",
            str(root),
            "--contract",
            str(contract_path),
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    if (
        result.returncode != 0
        or "R2_0D_CALIBRATION_CONTRACT_BLOCKED" not in result.stdout
        or "R2_0D_CALIBRATION_CONTRACT_FAIL" in result.stdout + result.stderr
    ):
        raise ProtocolError(
            "upstream-contract-rejected:"
            f"rc={result.returncode}:{result.stdout.strip()}:{result.stderr.strip()}"
        )

    contract = load_json(contract_path, "contract")
    reference = load_json(reference_path, "reference")
    artifacts = contract.get("artifacts")
    if not isinstance(artifacts, list):
        raise ProtocolError("contract-artifacts-invalid")
    reference_artifacts = [
        artifact
        for artifact in artifacts
        if isinstance(artifact, dict) and artifact.get("id") == "reference_output"
    ]
    if len(reference_artifacts) != 1:
        raise ProtocolError("contract-reference-artifact-invalid")
    expected_reference = (root / reference_artifacts[0]["path"]).resolve()
    if reference_path.resolve() != expected_reference:
        raise ProtocolError("reference-path-not-bound-by-contract")
    if reference_artifacts[0].get("sha256") != sha256(reference_path):
        raise ProtocolError("reference-sha256-not-bound-by-contract")
    if contract.get("status") != "blocked":
        raise ProtocolError("upstream-contract-not-blocked")
    approvals = contract.get("approvals")
    if not isinstance(approvals, dict) or not approvals or any(
        value is not False for value in approvals.values()
    ):
        raise ProtocolError("upstream-approval-enabled")
    protocol = contract.get("protocol")
    if not isinstance(protocol, dict):
        raise ProtocolError("upstream-protocol-invalid")
    if protocol.get("rfam_0d_runner_status") != "not-implemented":
        raise ProtocolError("upstream-runner-baseline-status-mismatch")
    if protocol.get("common_basis_status") != (
        "not-established-for-rfam-vs-cantera-trajectory-comparison"
    ):
        raise ProtocolError("upstream-common-basis-status-mismatch")

    case_ids = protocol.get("protocol_training_case_ids", []) + protocol.get(
        "protocol_holdout_case_ids", []
    )
    cases = reference.get("cases")
    if not isinstance(cases, list) or [case.get("id") for case in cases] != case_ids:
        raise ProtocolError("reference-case-order-mismatch")
    if reference.get("major_species") != SPECIES:
        raise ProtocolError("reference-species-order-mismatch")
    sample_times = protocol.get("reactor", {}).get("sample_times_s")
    if sample_times != reference.get("sample_times_s"):
        raise ProtocolError("reference-sample-times-mismatch")
    return contract, reference


def parse_sample(line: str, expected_index: int) -> dict[str, Any]:
    fields: dict[str, str] = {}
    for token in line.split("|"):
        if "=" not in token:
            raise ProtocolError(f"runner-sample-token-invalid:{expected_index}")
        key, value = token.split("=", 1)
        if not key or key in fields:
            raise ProtocolError(f"runner-sample-key-invalid:{expected_index}:{key}")
        fields[key] = value
    if set(fields) != SAMPLE_KEYS:
        raise ProtocolError(f"runner-sample-shape-mismatch:{expected_index}")
    if integer_text(fields["sample"], "sample-index") != expected_index:
        raise ProtocolError(f"runner-sample-index-mismatch:{expected_index}")
    composition_text = fields.pop("mass_fractions")
    composition = [
        finite_number(item, f"sample[{expected_index}].Y")
        for item in composition_text.split(",")
    ]
    if len(composition) != len(SPECIES):
        raise ProtocolError(f"runner-composition-size-mismatch:{expected_index}")
    result: dict[str, Any] = {
        key: finite_number(value, f"sample[{expected_index}].{key}")
        for key, value in fields.items()
        if key != "sample"
    }
    result["mass_fractions"] = composition
    return result


def parse_runner_output(
    stdout: str,
    stderr: str,
    case_id: str,
    sample_count: int,
) -> tuple[dict[str, str], list[dict[str, Any]]]:
    if stderr.strip():
        raise ProtocolError(f"runner-stderr-not-empty:{case_id}:{stderr.strip()}")
    lines = stdout.splitlines()
    if not lines or lines[-1] != "R2_0D_RFAM_RUNNER_OK":
        raise ProtocolError(f"runner-success-marker-missing:{case_id}")
    if lines.count("R2_0D_RFAM_RUNNER_OK") != 1:
        raise ProtocolError(f"runner-success-marker-count:{case_id}")
    headers: dict[str, str] = {}
    sample_lines: list[str] = []
    for line in lines[:-1]:
        if line.startswith("sample="):
            sample_lines.append(line)
            continue
        if "=" not in line:
            raise ProtocolError(f"runner-header-line-invalid:{case_id}")
        key, value = line.split("=", 1)
        if not key or key in headers:
            raise ProtocolError(f"runner-header-key-invalid:{case_id}:{key}")
        headers[key] = value
    if set(headers) != HEADER_KEYS:
        raise ProtocolError(f"runner-header-shape-mismatch:{case_id}")
    expected_headers = {
        "schema": "rfam.r2-0d-rfam-trajectory/v1",
        "status": "implementation-protocol-only-not-fit-eligible",
        "case_id": case_id,
        "reactor": "ideal-gas-constant-pressure-sensible-enthalpy",
        "chemistry": "rfam-twostepco-candidate",
        "energy_basis": "sensible-enthalpy-plus-explicit-formation-heat",
        "max_internal_step_s": "0",
        "species_order": ",".join(SPECIES),
        "ignition_definition": "first-root-T-minus-T0-equals-400K",
        "positivity_strategy": (
            "fuel-depletion-root-plus-rejected-interval-bisection"
        ),
        "implicit_time_extension_s": "0",
        "fit_approved": "0",
        "physical_validation_approved": "0",
        "production_ready": "0",
    }
    for key, value in expected_headers.items():
        if headers[key] != value:
            raise ProtocolError(f"runner-header-value-mismatch:{case_id}:{key}")
    if integer_text(headers["sample_count"], "sample-count") != sample_count:
        raise ProtocolError(f"runner-sample-count-mismatch:{case_id}")
    if len(sample_lines) != sample_count:
        raise ProtocolError(f"runner-sample-line-count-mismatch:{case_id}")
    ignition_detected = integer_text(
        headers["ignition_detected"], "ignition-detected"
    )
    if ignition_detected not in (0, 1):
        raise ProtocolError(f"runner-ignition-flag-invalid:{case_id}")
    ignition_delay = finite_number(headers["ignition_delay_s"], "ignition-delay")
    if (ignition_detected == 1) != (ignition_delay > 0.0):
        raise ProtocolError(f"runner-ignition-delay-invalid:{case_id}")
    for count_key in (
        "depletion_event_count",
        "interval_bisection_count",
        "roundoff_clamp_count",
    ):
        count = integer_text(headers[count_key], count_key)
        if count > 1000000:
            raise ProtocolError(f"runner-diagnostic-count-invalid:{case_id}:{count_key}")
    finite_number(headers["accepted_physical_time_s"], "accepted-physical-time")
    samples = [
        parse_sample(line, index) for index, line in enumerate(sample_lines)
    ]
    return headers, samples


def elemental_inventory(composition: list[float]) -> list[float]:
    return [
        sum(
            composition[species]
            / MOLECULAR_WEIGHTS[species]
            * ELEMENT_ATOMS[species][element]
            for species in range(len(SPECIES))
        )
        for element in range(4)
    ]


def validate_case(
    runner: Path,
    case: dict[str, Any],
    sample_times: list[float],
) -> tuple[dict[str, Any], int]:
    case_id = case["id"]
    initial = case.get("initial_major_species_mass_fractions")
    if not isinstance(initial, dict) or set(initial) != set(SPECIES):
        raise ProtocolError(f"reference-initial-composition-invalid:{case_id}")
    initial_composition = [
        finite_number(initial[species], f"{case_id}.initial.{species}")
        for species in SPECIES
    ]
    command = [
        str(runner),
        case_id,
        format(finite_number(case["temperature_k"], "initial-temperature"), ".17g"),
        format(finite_number(case["pressure_pa"], "initial-pressure"), ".17g"),
        str(len(sample_times)),
        "0",
        *[format(value, ".17g") for value in initial_composition],
        *[format(finite_number(value, "sample-time"), ".17g") for value in sample_times],
    ]
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    if result.returncode != 0:
        raise ProtocolError(
            f"runner-rejected:{case_id}:rc={result.returncode}:"
            f"{result.stdout.strip()}:{result.stderr.strip()}"
        )
    headers, samples = parse_runner_output(
        result.stdout, result.stderr, case_id, len(sample_times)
    )

    checks = 0
    initial_elements = elemental_inventory(initial_composition)
    initial_total_enthalpy = samples[0]["total_enthalpy_j_kg"]
    pressure = finite_number(case["pressure_pa"], "case-pressure")
    for index, (sample, expected_time) in enumerate(zip(samples, sample_times)):
        require_close(sample["time_s"], expected_time, 0.0, 2.0e-18, f"{case_id}.time")
        require_close(sample["pressure_pa"], pressure, 0.0, 1.0e-9, f"{case_id}.pressure")
        if not minimum_state_envelope(sample):
            raise ProtocolError(f"runner-state-out-of-envelope:{case_id}:{index}")
        composition = sample["mass_fractions"]
        if any(value < 0.0 for value in composition):
            raise ProtocolError(f"runner-negative-composition:{case_id}:{index}")
        require_close(sum(composition), 1.0, 0.0, 2.0e-9, f"{case_id}.sumY")
        inverse_molecular_weight = sum(
            value / molecular_weight
            for value, molecular_weight in zip(composition, MOLECULAR_WEIGHTS)
        )
        expected_density = pressure / (
            GAS_CONSTANT * inverse_molecular_weight * sample["temperature_k"]
        )
        require_close(
            sample["density_kg_m3"], expected_density, 2.0e-12, 1.0e-12,
            f"{case_id}.ideal-density",
        )
        current_elements = elemental_inventory(composition)
        for element, (actual, expected) in enumerate(
            zip(current_elements, initial_elements)
        ):
            require_close(actual, expected, 2.0e-10, 2.0e-12, f"{case_id}.element[{element}]")
        require_close(
            sample["total_enthalpy_j_kg"],
            initial_total_enthalpy,
            2.0e-10,
            1.0e-5,
            f"{case_id}.total-enthalpy",
        )
        checks += 9

    require_close(
        finite_number(
            headers["accepted_physical_time_s"], "accepted-physical-time"
        ),
        sample_times[-1],
        0.0,
        2.0e-18,
        f"{case_id}.accepted-physical-time",
    )
    require_close(
        finite_number(
            headers["implicit_time_extension_s"], "implicit-time-extension"
        ),
        0.0,
        0.0,
        0.0,
        f"{case_id}.implicit-time-extension",
    )
    checks += 2

    require_close(
        samples[0]["temperature_k"],
        finite_number(case["temperature_k"], "initial-temperature"),
        0.0,
        2.0e-10,
        f"{case_id}.initial-temperature",
    )
    for species, (actual, expected) in enumerate(
        zip(samples[0]["mass_fractions"], initial_composition)
    ):
        require_close(actual, expected, 0.0, 2.0e-12, f"{case_id}.initial-Y[{species}]")
    reference_density = finite_number(
        case["samples"][0]["density_kg_m3"], "reference-initial-density"
    )
    initial_density_gap = abs(samples[0]["density_kg_m3"] - reference_density) / reference_density
    # This is only an input-mapping sanity bound.  It is intentionally not a
    # fit metric; RFAM and Cantera molecular-weight authorities are not yet a
    # frozen common comparison basis.
    if initial_density_gap > 1.0e-3:
        raise ProtocolError(f"initial-density-mapping-gap:{case_id}")
    checks += 8

    return {
        "id": case_id,
        "partition": case["partition"],
        "runner_headers": headers,
        "initial_density_relative_mapping_gap": initial_density_gap,
        "samples": samples,
    }, checks


def minimum_state_envelope(sample: dict[str, Any]) -> bool:
    return (
        200.0 <= sample["temperature_k"] <= 8000.0
        and sample["pressure_pa"] > 0.0
        and sample["density_kg_m3"] > 0.0
        and math.isfinite(sample["sensible_enthalpy_j_kg"])
        and math.isfinite(sample["total_enthalpy_j_kg"])
        and math.isfinite(sample["qdot_w_m3"])
    )


def run_protocol(args: argparse.Namespace) -> tuple[dict[str, Any], int]:
    root = args.root.resolve()
    contract_path = args.contract.resolve()
    reference_path = args.reference.resolve()
    runner = args.runner.resolve()
    if not runner.is_file():
        raise ProtocolError(f"runner-open-failed:{runner}")
    contract, reference = validate_upstream(root, contract_path, reference_path)
    sample_times = [
        finite_number(value, "reference-sample-time")
        for value in reference["sample_times_s"]
    ]
    output_cases = []
    checks = 0
    for case in reference["cases"]:
        output, case_checks = validate_case(runner, case, sample_times)
        output_cases.append(output)
        checks += case_checks
    output_document = {
        "schema": "rfam.r2-0d-rfam-protocol-output/v1",
        "status": "implementation-protocol-only-no-fit-or-physical-validation",
        "upstream": {
            "contract_path": contract_path.relative_to(root).as_posix(),
            "contract_sha256": sha256(contract_path),
            "reference_path": reference_path.relative_to(root).as_posix(),
            "reference_sha256": sha256(reference_path),
            "common_basis_status": contract["protocol"]["common_basis_status"],
        },
        "runner": {
            "path": str(runner),
            "sha256": sha256(runner),
            "reactor": "ideal-gas-constant-pressure-sensible-enthalpy",
            "chemistry": "rfam-twostepco-candidate",
            "case_count": len(output_cases),
            "sample_count_per_case": len(sample_times),
        },
        "comparison": {
            "cantera_numeric_fit_score_generated": False,
            "ignition_delay_fit_metric_enabled": False,
            "trajectory_fit_metrics_enabled": False,
            "heat_release_fit_metric_enabled": False,
            "reason": "detailed mechanism, experiment, tolerances and common basis remain unselected",
        },
        "approvals": {
            "fit_approved": False,
            "physical_validation_approved": False,
            "production_ready": False,
            "reacting_hot_start_approved": False,
        },
        "cases": output_cases,
    }
    return output_document, checks


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    try:
        document, checks = run_protocol(args)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps(
                    document,
                    sort_keys=True,
                    separators=(",", ":"),
                    ensure_ascii=True,
                )
                + "\n",
                encoding="utf-8",
            )
    except (ProtocolError, OSError, ValueError) as error:
        print(f"R2_0D_RFAM_PROTOCOL_FAIL {error}", file=sys.stderr)
        return 2
    except Exception as error:  # fail closed without exposing a traceback
        print(
            "R2_0D_RFAM_PROTOCOL_FAIL unexpected:"
            + type(error).__name__,
            file=sys.stderr,
        )
        return 2
    print(
        "R2_0D_RFAM_PROTOCOL_BLOCKED "
        f"cases={len(document['cases'])} "
        f"samples={sum(len(case['samples']) for case in document['cases'])} "
        f"checks={checks}"
    )
    print("fit_approved=0 physical_validation_approved=0 production_ready=0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
