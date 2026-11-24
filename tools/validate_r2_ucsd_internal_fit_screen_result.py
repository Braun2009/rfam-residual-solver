#!/usr/bin/env python3
"""Validate the frozen UCSD internal fit-screen package without executing it."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence


EXPECTED_RESULT_SHA256 = "b95278870930eacad51ee6870223627f6e715090c0e1e807c349514f5e8bf353"
EXPECTED_RESULT_PATHS = {
    "F:/tmp/a408-r2-ucsd-internal-fit-screen-result-v4.json",
    "/mnt/f/tmp/a408-r2-ucsd-internal-fit-screen-result-v4.json",
}
EXPECTED_EXTERNAL = {
    "training_reference": "bcc2a6251384734f44d5b31cbeb10ebd5a77ae265eacdb0e1aaba82040c6d71d",
    "holdout_reference": "856f86914e61607c4a5f3113b72fbb79d3a68562162a2b0bf43a29b4c3865def",
    "parameterized_runner": "d080d0ed5a17fb029fc5e01fbad5dbe6f728f88f0686f0292606a8e81344f9ad",
}
EXPECTED_REPOSITORY = {
    "tools/run_r2_ucsd_internal_fit_screen.py": "2b230c4fa4434ad9506a06019617cca2dec18a7737093e5589788cc0ce2ff407",
    "tests/test_r2_ucsd_internal_fit_screen.py": "70790e90df7a2119de7521d10ec17df43b47355144cf342e4c080298e11debe9",
    "tests/reference/r2_ucsd_internal_fit_screen_contract_v1.json": "8cdf5bb9b11fccbaf4aafef07840bc0c5cbe22a3b99d9ff685eef5e62d703796",
    "tools/r2_0d_rfam_fit_runner.cpp": "ee0accda3f6d414f0757ec3723c24f93a0da7560714ab866682ea9fecf2c2c73",
    "tools/validate_r2_parameterized_runner.py": "0a1c89f10f265ccf23c753eaeefec32fcb4c86864ca09504cd8fbd3672fe0271",
    "tools/split_r2_ucsd_internal_reference.py": "0372b833453c8d8c479238a426f17bad5564ddaa18765a071d13b09baa421985",
    "tools/run_r2_0d_rfam_protocol.py": "5a4642c37ce4b07ca97db82a6af218d39e039b8a73e417d20efb6094d43f58a8",
    "tools/run_r2_joint_fit_synthetic.py": "0669e6d7c2d81fc42a71dfd96688fdb13d1c821dbe953c8a6f1f9bbefd40f99d",
    "tools/screen_r2_detailed_candidate.py": "6e88e92ebe84db8ae378f95399de4764ca55f122ef9156ffb5c42e48874a14e7",
}
EXPECTED_SOURCE = {
    "contract_path": "tests/reference/r2_ucsd_internal_fit_screen_contract_v1.json",
    "contract_sha256": EXPECTED_REPOSITORY[
        "tests/reference/r2_ucsd_internal_fit_screen_contract_v1.json"
    ],
    "execution": {
        "platform_machine": "x86_64",
        "platform_system": "Linux",
        "python_version": "3.8.10",
        "runner_timeout_seconds": 30,
    },
    "training_reference": {
        "logical_path": "internal://a408-r2-ucsd/training-v1.json",
        "sha256": EXPECTED_EXTERNAL["training_reference"],
    },
    "holdout_reference": {
        "logical_path": "internal://a408-r2-ucsd/holdout-v1.json",
        "sha256": EXPECTED_EXTERNAL["holdout_reference"],
    },
    "parameterized_runner": {
        "logical_path": "build://rfam-r2-0d-fit-runner/linux-x86_64-v1",
        "sha256": EXPECTED_EXTERNAL["parameterized_runner"],
    },
    "repository_artifacts": [
        {"path": path, "sha256": EXPECTED_REPOSITORY[path]}
        for path in (
            "tools/r2_0d_rfam_fit_runner.cpp",
            "tools/validate_r2_parameterized_runner.py",
            "tools/split_r2_ucsd_internal_reference.py",
            "tools/run_r2_0d_rfam_protocol.py",
            "tools/run_r2_joint_fit_synthetic.py",
            "tools/screen_r2_detailed_candidate.py",
        )
    ],
}
EXPECTED_CLAIMS = {
    "experimental_fit_executed": False,
    "optimizer_executed_against_internal_detailed_candidate": True,
    "parameter_vector_publishable": False,
    "parameterized_runner_default_byte_equivalence_is_separate_required_evidence": True,
    "physical_validation_claimed": False,
    "production_ready": False,
}
EXPECTED_LIMITS = [
    "optimizer output is an internal model-structure diagnostic only",
    "protocol training and holdout labels are not experimental validation partitions",
    "a low objective cannot cure missing radical species or grant mechanism rights",
    "holdout data cannot be read by the objective during parameter selection",
    "no parameter may enter runtime defaults or a production manifest",
]
TRAINING_CASES = ["train_phi0p8_diluted", "train_phi1p0", "train_phi1p2"]
HOLDOUT_CASES = ["holdout_phi0p6_diluted", "holdout_phi1p1", "holdout_phi1p4"]
DEFAULT_PARAMETERS = [10000000000.0, 125520.0, 0.1, 1.65, 2000000000000.0, 167360.0]
PARAMETER_POLICIES = [
    ("log10", 100000.0, 1.0e20),
    ("linear", 0.0, 300000.0),
    ("linear", 0.0, 2.0),
    ("linear", 0.0, 2.5),
    ("log10", 100000.0, 1.0e20),
    ("linear", 0.0, 300000.0),
]


class ValidationError(RuntimeError):
    pass


def strict_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValidationError("duplicate-json-key:" + key)
        result[key] = value
    return result


def reject_constant(token: str) -> None:
    raise ValidationError("nonfinite-json-token:" + token)


def load_json(path: Path) -> tuple[dict[str, Any], bytes]:
    try:
        raw = path.read_bytes()
        document = json.loads(
            raw.decode("ascii"),
            object_pairs_hook=strict_pairs,
            parse_constant=reject_constant,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValidationError("result-read-failed:" + str(error)) from error
    if not isinstance(document, dict):
        raise ValidationError("result-root-not-object")
    return document, raw


def digest_file(path: Path) -> str:
    result = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                result.update(block)
    except OSError as error:
        raise ValidationError("artifact-read-failed:" + str(path)) from error
    return result.hexdigest()


def exact_keys(value: Any, expected: Sequence[str], label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValidationError("not-object:" + label)
    if set(value) != set(expected):
        raise ValidationError("key-schema-mismatch:" + label)
    return value


def exact_tree(actual: Any, expected: Any, label: str) -> None:
    if type(actual) is not type(expected):
        raise ValidationError("type-mismatch:" + label)
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise ValidationError("key-schema-mismatch:" + label)
        for key in expected:
            exact_tree(actual[key], expected[key], label + "." + key)
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValidationError("list-length-mismatch:" + label)
        for index, (left, right) in enumerate(zip(actual, expected)):
            exact_tree(left, right, label + "[" + str(index) + "]")
    elif actual != expected:
        raise ValidationError("value-mismatch:" + label)


def strict_int(value: Any, label: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValidationError("type-mismatch:" + label)
    if value < minimum:
        raise ValidationError("range-mismatch:" + label)
    return value


def strict_float(value: Any, label: str, minimum: float = 0.0) -> float:
    if not isinstance(value, float):
        raise ValidationError("type-mismatch:" + label)
    if not math.isfinite(value) or value < minimum:
        raise ValidationError("range-mismatch:" + label)
    return value


def require_sha(value: Any, label: str) -> None:
    if not isinstance(value, str) or len(value) != 64 or any(
        character not in "0123456789abcdef" for character in value
    ):
        raise ValidationError("sha256-invalid:" + label)


def case_objectives(value: Any, case_ids: Sequence[str], label: str) -> dict[str, Any]:
    result = exact_keys(value, case_ids, label)
    for case_id in case_ids:
        strict_float(result[case_id], label + "." + case_id)
    return result


def close_sum(total: float, components: Mapping[str, Any], label: str) -> None:
    expected = math.fsum(float(value) for value in components.values())
    tolerance = max(1.0e-9, 5.0e-12 * max(abs(total), abs(expected), 1.0))
    if abs(total - expected) > tolerance:
        raise ValidationError("objective-sum-mismatch:" + label)


def canonical_bytes(document: dict[str, Any]) -> bytes:
    return (
        json.dumps(
            document,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("ascii")


def validate_parameter_mapping(
    coordinates: Any,
    parameters: Any,
    label: str,
) -> tuple[float, ...]:
    if not isinstance(coordinates, list) or len(coordinates) != 6:
        raise ValidationError(label + "-coordinate-shape-mismatch")
    if not isinstance(parameters, list) or len(parameters) != 6:
        raise ValidationError(label + "-parameter-shape-mismatch")
    transformed_key = []
    for index, (coordinate_value, actual_value, policy) in enumerate(
        zip(coordinates, parameters, PARAMETER_POLICIES)
    ):
        coordinate = strict_float(
            coordinate_value, label + ".coordinate[" + str(index) + "]"
        )
        actual = strict_float(
            actual_value, label + ".parameter[" + str(index) + "]"
        )
        if coordinate > 1.0:
            raise ValidationError(label + "-coordinate-range-mismatch")
        transform, lower, upper = policy
        transformed = (
            math.log10(lower)
            + coordinate * (math.log10(upper) - math.log10(lower))
            if transform == "log10"
            else lower + coordinate * (upper - lower)
        )
        expected = 10.0 ** transformed if transform == "log10" else transformed
        tolerance = max(1.0e-12, 5.0e-12 * max(abs(actual), abs(expected), 1.0))
        if abs(actual - expected) > tolerance:
            raise ValidationError(label + "-parameter-transform-mismatch:" + str(index))
        transformed_key.append(transformed)
    return tuple(transformed_key)


def validate_document(document: dict[str, Any]) -> int:
    checks = 0
    exact_keys(
        document,
        [
            "baseline", "claims", "limits", "runner_accounting", "schema",
            "selection", "source", "starts", "status",
        ],
        "root",
    )
    if document["schema"] != "rfam.r2-ucsd-internal-fit-screen-result/v1":
        raise ValidationError("schema-mismatch")
    if document["status"] != "blocked-internal-route-screen-only":
        raise ValidationError("status-mismatch")
    checks += 3

    exact_tree(document["source"], EXPECTED_SOURCE, "source")
    exact_tree(document["claims"], EXPECTED_CLAIMS, "claims")
    exact_tree(document["limits"], EXPECTED_LIMITS, "limits")
    checks += 3

    baseline = exact_keys(
        document["baseline"],
        ["parameter_vector", "training_case_objectives", "training_objective"],
        "baseline",
    )
    exact_tree(baseline["parameter_vector"], DEFAULT_PARAMETERS, "baseline.parameter_vector")
    baseline_cases = case_objectives(
        baseline["training_case_objectives"], TRAINING_CASES, "baseline.training_case_objectives"
    )
    baseline_objective = strict_float(
        baseline["training_objective"], "baseline.training_objective"
    )
    close_sum(baseline_objective, baseline_cases, "baseline")
    checks += 4

    selection = exact_keys(
        document["selection"],
        [
            "holdout_case_objectives", "holdout_objective",
            "selected_normalized_coordinates", "selected_parameter_vector",
            "selected_start_index", "selected_start_kind", "start_count",
            "training_case_objectives", "training_objective", "valid_start_count",
        ],
        "selection",
    )
    training_cases = case_objectives(
        selection["training_case_objectives"], TRAINING_CASES,
        "selection.training_case_objectives",
    )
    holdout_cases = case_objectives(
        selection["holdout_case_objectives"], HOLDOUT_CASES,
        "selection.holdout_case_objectives",
    )
    training_objective = strict_float(
        selection["training_objective"], "selection.training_objective"
    )
    holdout_objective = strict_float(
        selection["holdout_objective"], "selection.holdout_objective"
    )
    close_sum(training_objective, training_cases, "selection.training")
    close_sum(holdout_objective, holdout_cases, "selection.holdout")
    coordinates = selection["selected_normalized_coordinates"]
    parameters = selection["selected_parameter_vector"]
    validate_parameter_mapping(coordinates, parameters, "selected")
    selected_index = strict_int(selection["selected_start_index"], "selection.selected_start_index")
    if not isinstance(selection["selected_start_kind"], str):
        raise ValidationError("type-mismatch:selection.selected_start_kind")
    if strict_int(selection["start_count"], "selection.start_count", 1) != 5:
        raise ValidationError("start-count-mismatch")
    valid_start_count = strict_int(
        selection["valid_start_count"], "selection.valid_start_count", 1
    )
    checks += 12

    starts = document["starts"]
    if not isinstance(starts, list) or len(starts) != 5:
        raise ValidationError("starts-shape-mismatch")
    expected_kinds = ["midpoint", "halton-1", "halton-2", "halton-3", "halton-4"]
    evaluation_total = 0
    valid_count = 0
    transformed_keys: list[tuple[float, ...]] = []
    for index, (start, expected_kind) in enumerate(zip(starts, expected_kinds)):
        entry = exact_keys(
            start,
            [
                "accepted_coordinates", "cache_hits", "completed_sweeps",
                "contractions", "evaluation_count", "index", "kind", "objective",
                "selected_normalized_coordinates", "selected_parameter_vector",
                "termination", "trace_sha256", "valid",
            ],
            "starts[" + str(index) + "]",
        )
        if strict_int(entry["index"], "start.index") != index or entry["kind"] != expected_kind:
            raise ValidationError("start-identity-mismatch:" + str(index))
        if entry["valid"] is not True:
            raise ValidationError("start-validity-mismatch:" + str(index))
        valid_count += 1
        objective = strict_float(entry["objective"], "start.objective")
        transformed_keys.append(validate_parameter_mapping(
            entry["selected_normalized_coordinates"],
            entry["selected_parameter_vector"],
            "start[" + str(index) + "]",
        ))
        evaluation_count = strict_int(entry["evaluation_count"], "start.evaluation_count", 1)
        if evaluation_count > 800:
            raise ValidationError("start-evaluation-budget-mismatch:" + str(index))
        evaluation_total += evaluation_count
        strict_int(entry["cache_hits"], "start.cache_hits")
        strict_int(entry["accepted_coordinates"], "start.accepted_coordinates")
        strict_int(entry["completed_sweeps"], "start.completed_sweeps")
        if strict_int(entry["contractions"], "start.contractions") != 4:
            raise ValidationError("start-contraction-mismatch:" + str(index))
        if entry["termination"] != "minimum-step-reached":
            raise ValidationError("start-termination-mismatch:" + str(index))
        require_sha(entry["trace_sha256"], "start.trace_sha256")
        if index == selected_index:
            if entry["kind"] != selection["selected_start_kind"]:
                raise ValidationError("selected-start-kind-mismatch")
            if objective != training_objective:
                raise ValidationError("selected-start-objective-mismatch")
            if (
                entry["selected_normalized_coordinates"] != coordinates
                or entry["selected_parameter_vector"] != parameters
            ):
                raise ValidationError("selected-start-parameter-mismatch")
        checks += 12
    if valid_count != valid_start_count:
        raise ValidationError("valid-start-count-mismatch")
    best_index = min(
        range(len(starts)),
        key=lambda index: (starts[index]["objective"], transformed_keys[index]),
    )
    if best_index != selected_index:
        raise ValidationError("global-selection-mismatch")
    checks += 2

    accounting = exact_keys(
        document["runner_accounting"],
        [
            "holdout_access_attempts_during_selection", "holdout_case_process_attempts",
            "holdout_case_process_successes", "holdout_objective_calls_after_selection",
            "invalid_candidate_reasons", "training_case_process_attempts",
            "training_case_process_successes", "training_objective_calls",
        ],
        "runner_accounting",
    )
    values = {
        key: strict_int(accounting[key], "runner_accounting." + key)
        for key in accounting
        if key != "invalid_candidate_reasons"
    }
    invalid = exact_keys(
        accounting["invalid_candidate_reasons"],
        ["ignition-not-detected:train_phi0p8_diluted"],
        "runner_accounting.invalid_candidate_reasons",
    )
    invalid_total = strict_int(
        invalid["ignition-not-detected:train_phi0p8_diluted"],
        "runner_accounting.invalid-count",
    )
    expected_training_calls = 1 + evaluation_total
    if values["training_objective_calls"] != expected_training_calls:
        raise ValidationError("training-objective-accounting-mismatch")
    if values["training_case_process_attempts"] != 3 * expected_training_calls:
        raise ValidationError("training-attempt-accounting-mismatch")
    if values["training_case_process_successes"] != values["training_case_process_attempts"]:
        raise ValidationError("training-success-accounting-mismatch")
    if (
        values["holdout_access_attempts_during_selection"] != 0
        or values["holdout_objective_calls_after_selection"] != 1
        or values["holdout_case_process_attempts"] != 3
        or values["holdout_case_process_successes"] != 3
    ):
        raise ValidationError("holdout-accounting-mismatch")
    if invalid_total > expected_training_calls:
        raise ValidationError("invalid-candidate-accounting-mismatch")
    checks += 7
    return checks


def verify_repository(root: Path) -> int:
    if not root.is_dir():
        raise ValidationError("repository-root-missing")
    for relative, expected_sha256 in EXPECTED_REPOSITORY.items():
        path = root / relative
        if not path.is_file() or digest_file(path) != expected_sha256:
            raise ValidationError("repository-artifact-mismatch:" + relative)
    return len(EXPECTED_REPOSITORY)


def verify_external(
    training_reference: Path,
    holdout_reference: Path,
    runner: Path,
) -> int:
    for label, path in (
        ("training_reference", training_reference),
        ("holdout_reference", holdout_reference),
        ("parameterized_runner", runner),
    ):
        if not path.is_file() or digest_file(path) != EXPECTED_EXTERNAL[label]:
            raise ValidationError("external-artifact-mismatch:" + label)
    return 3


def validate_package(
    root: Path,
    result_path: Path,
    training_reference: Path,
    holdout_reference: Path,
    runner: Path,
) -> int:
    document, raw = load_json(result_path)
    checks = validate_document(document)
    if canonical_bytes(document) != raw:
        raise ValidationError("result-not-canonical")
    if hashlib.sha256(raw).hexdigest() != EXPECTED_RESULT_SHA256:
        raise ValidationError("result-sha256-mismatch")
    if result_path.resolve().as_posix() not in EXPECTED_RESULT_PATHS:
        raise ValidationError("result-path-mismatch")
    checks += 3
    checks += verify_repository(root)
    checks += verify_external(training_reference, holdout_reference, runner)
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--training-reference", type=Path, required=True)
    parser.add_argument("--holdout-reference", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        checks = validate_package(
            arguments.root.resolve(),
            arguments.result.resolve(),
            arguments.training_reference.resolve(),
            arguments.holdout_reference.resolve(),
            arguments.runner.resolve(),
        )
    except Exception as error:
        print("R2_UCSD_INTERNAL_FIT_SCREEN_RESULT_FAIL " + str(error), file=sys.stderr)
        return 2
    print(
        "R2_UCSD_INTERNAL_FIT_SCREEN_RESULT_BLOCKED checks=" + str(checks)
        + " physical_validation=0 production_ready=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
