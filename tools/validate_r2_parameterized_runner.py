#!/usr/bin/env python3
"""Validate the isolated parameterized RFAM R1/R2 trajectory runner."""

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
DEFAULT_PARAMETERS = [1.0e10, 125520.0, 0.1, 1.65, 2.0e12, 167360.0]
PARAMETER_NAMES = [
    "R1_A_SI",
    "R1_Ea_J_per_mol",
    "R1_order_C3H8",
    "R1_order_O2",
    "R2_A_SI",
    "R2_Ea_J_per_mol",
]


class ValidationError(RuntimeError):
    pass


def object_without_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValidationError("duplicate-json-key:" + key)
        result[key] = value
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_hash(path: Path, expected: str, label: str) -> None:
    if not path.is_file():
        raise ValidationError("missing-artifact:" + label)
    actual = sha256(path)
    if actual != expected:
        raise ValidationError(
            "sha256-mismatch:" + label + ":expected=" + expected + ":actual=" + actual
        )


def finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValidationError("non-numeric:" + label)
    result = float(value)
    if not math.isfinite(result):
        raise ValidationError("nonfinite:" + label)
    return result


def load_reference(path: Path) -> tuple[list[float], list[dict[str, Any]]]:
    try:
        document = json.loads(
            path.read_text(encoding="ascii"),
            object_pairs_hook=object_without_duplicates,
            parse_constant=lambda token: (_ for _ in ()).throw(
                ValidationError("nonfinite-json-token:" + token)
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValidationError("reference-read-failed:" + str(error)) from error
    if not isinstance(document, dict):
        raise ValidationError("reference-root-not-object")
    if document.get("schema") != "rfam.r2-0d-cantera-reference/v1":
        raise ValidationError("reference-schema-mismatch")
    if document.get("major_species") != SPECIES:
        raise ValidationError("reference-species-mismatch")
    times_raw = document.get("sample_times_s")
    cases = document.get("cases")
    if not isinstance(times_raw, list) or not isinstance(cases, list) or not cases:
        raise ValidationError("reference-shape-invalid")
    times = [finite(value, "sample-time") for value in times_raw]
    if len(times) < 2 or times[0] != 0.0 or any(
        right <= left for left, right in zip(times, times[1:])
    ):
        raise ValidationError("sample-times-invalid")
    case_ids: set[str] = set()
    for case in cases:
        if not isinstance(case, dict):
            raise ValidationError("reference-case-not-object")
        case_id = case.get("id")
        initial = case.get("initial_major_species_mass_fractions")
        if (
            not isinstance(case_id, str)
            or not case_id
            or case_id in case_ids
            or not isinstance(initial, dict)
            or set(initial) != set(SPECIES)
        ):
            raise ValidationError("reference-case-identity-invalid")
        case_ids.add(case_id)
        finite(case.get("temperature_k"), case_id + ".temperature")
        finite(case.get("pressure_pa"), case_id + ".pressure")
        values = [finite(initial[name], case_id + "." + name) for name in SPECIES]
        if any(value < 0.0 for value in values) or abs(math.fsum(values) - 1.0) > 2e-9:
            raise ValidationError("reference-case-composition-invalid:" + case_id)
    return times, cases


def trajectory_arguments(case: dict[str, Any], times: list[float]) -> list[str]:
    initial = case["initial_major_species_mass_fractions"]
    return [
        case["id"],
        format(float(case["temperature_k"]), ".17g"),
        format(float(case["pressure_pa"]), ".17g"),
        str(len(times)),
        "0",
        *[format(float(initial[name]), ".17g") for name in SPECIES],
        *[format(value, ".17g") for value in times],
    ]


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=False, text=True, capture_output=True)


def parameter_prefix(parameters: list[float]) -> str:
    lines = ["fit_parameter_schema=rfam.r2-joint-fit-parameters/v1"]
    lines.extend(
        "fit_" + name + "=" + format(value, ".17g")
        for name, value in zip(PARAMETER_NAMES, parameters)
    )
    return "\n".join(lines) + "\n"


def validate(
    reference: Path,
    fixed_runner: Path,
    fit_runner: Path,
    expected_reference_sha256: str,
    expected_fixed_sha256: str,
    expected_fit_sha256: str,
) -> int:
    checks = 0
    for value, label in (
        (expected_reference_sha256, "reference"),
        (expected_fixed_sha256, "fixed-runner"),
        (expected_fit_sha256, "fit-runner"),
    ):
        if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
            raise ValidationError("invalid-expected-sha256:" + label)
    require_hash(reference, expected_reference_sha256, "reference")
    require_hash(fixed_runner, expected_fixed_sha256, "fixed-runner")
    require_hash(fit_runner, expected_fit_sha256, "fit-runner")
    checks += 3

    times, cases = load_reference(reference)
    prefix = parameter_prefix(DEFAULT_PARAMETERS)
    for case in cases:
        arguments = trajectory_arguments(case, times)
        fixed = run([str(fixed_runner), *arguments])
        fitted = run(
            [
                str(fit_runner),
                *[format(value, ".17g") for value in DEFAULT_PARAMETERS],
                *arguments,
            ]
        )
        if fixed.returncode != 0 or fitted.returncode != 0:
            raise ValidationError("default-runner-failed:" + case["id"])
        if fixed.stderr != fitted.stderr:
            raise ValidationError("default-stderr-mismatch:" + case["id"])
        if not fitted.stdout.startswith(prefix):
            raise ValidationError("fit-parameter-prefix-mismatch:" + case["id"])
        if fitted.stdout[len(prefix) :] != fixed.stdout:
            raise ValidationError("default-trajectory-byte-mismatch:" + case["id"])
        checks += 3

    perturbed = list(DEFAULT_PARAMETERS)
    perturbed[0] *= 0.1
    first_arguments = trajectory_arguments(cases[0], times)
    default_run = run(
        [
            str(fit_runner),
            *[format(value, ".17g") for value in DEFAULT_PARAMETERS],
            *first_arguments,
        ]
    )
    perturbed_run = run(
        [
            str(fit_runner),
            *[format(value, ".17g") for value in perturbed],
            *first_arguments,
        ]
    )
    if default_run.returncode != 0 or perturbed_run.returncode != 0:
        raise ValidationError("perturbed-runner-failed")
    if not perturbed_run.stdout.startswith(parameter_prefix(perturbed)):
        raise ValidationError("perturbed-parameter-prefix-mismatch")
    default_tail = default_run.stdout[len(prefix) :]
    perturbed_prefix = parameter_prefix(perturbed)
    if perturbed_run.stdout[len(perturbed_prefix) :] == default_tail:
        raise ValidationError("parameter-perturbation-has-no-effect")
    checks += 3

    invalid_vectors = [
        ["true", "125520", "0.1", "1.65", "2e12", "167360"],
        ["99999", "125520", "0.1", "1.65", "2e12", "167360"],
        ["1e10", "-1", "0.1", "1.65", "2e12", "167360"],
        ["1e10", "125520", "2.1", "1.65", "2e12", "167360"],
        ["1e10", "125520", "0.1", "2.6", "2e12", "167360"],
        ["1e10", "125520", "0.1", "1.65", "1e21", "167360"],
        ["1e10", "125520", "0.1", "1.65", "2e12", "nan"],
    ]
    for index, parameters in enumerate(invalid_vectors):
        result = run([str(fit_runner), *parameters, *first_arguments])
        if result.returncode != 2 or "R2_FIT_RUNNER_ERROR" not in result.stderr:
            raise ValidationError("invalid-vector-not-rejected:" + str(index))
        if "schema=rfam.r2-0d-rfam-trajectory/v1" in result.stdout:
            raise ValidationError("invalid-vector-emitted-trajectory:" + str(index))
        checks += 2

    print(
        "R2_PARAMETERIZED_RUNNER_PASS cases="
        + str(len(cases))
        + " checks="
        + str(checks)
        + " default_byte_equivalent=1 physical_fit_executed=0"
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--fixed-runner", type=Path, required=True)
    parser.add_argument("--fit-runner", type=Path, required=True)
    parser.add_argument("--expected-reference-sha256", required=True)
    parser.add_argument("--expected-fixed-sha256", required=True)
    parser.add_argument("--expected-fit-sha256", required=True)
    arguments = parser.parse_args()
    try:
        return validate(
            arguments.reference.resolve(),
            arguments.fixed_runner.resolve(),
            arguments.fit_runner.resolve(),
            arguments.expected_reference_sha256,
            arguments.expected_fixed_sha256,
            arguments.expected_fit_sha256,
        )
    except (ValidationError, OSError) as error:
        print("R2_PARAMETERIZED_RUNNER_FAIL " + str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
