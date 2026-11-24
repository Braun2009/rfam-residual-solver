#!/usr/bin/env python3
"""Fail-closed package gate for deterministic synthetic R2 fit execution."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Sequence, Tuple


BLOCKED_MARKER = "R2_JOINT_FIT_SYNTHETIC_PACKAGE_BLOCKED"
FAIL_MARKER = "R2_JOINT_FIT_SYNTHETIC_PACKAGE_FAIL"

AUTHORITIES = {
    "execution_contract": {
        "path": "tests/reference/r2_joint_fit_synthetic_execution_contract_v1.json",
        "sha256": "4c3db0f2711a00d93d8bdff4b08a72df41107db1a930cd072cc10fa7d4b0e268",
        "byte_count": 6028,
    },
    "implementation": {
        "path": "tools/run_r2_joint_fit_synthetic.py",
        "sha256": "0669e6d7c2d81fc42a71dfd96688fdb13d1c821dbe953c8a6f1f9bbefd40f99d",
        "byte_count": 71874,
    },
    "focused_test": {
        "path": "tests/test_r2_joint_fit_synthetic_execution.py",
        "sha256": "a746890a0ccde9d28b0d7396835c71349c04b9ca545d368bd84544898cf39736",
        "byte_count": 28122,
    },
    "documentation": {
        "path": "docs/r2_joint_fit_synthetic_execution.md",
        "sha256": "fb673f2e0b7559102d07164ac4a7191d317af1d33d16e6f7adff598a97d94e7f",
        "byte_count": 5600,
    },
    "deterministic_result": {
        "path": "tests/reference/r2_joint_fit_synthetic_result_v1.json",
        "sha256": "582a00acc9779f0de4d1be2ecbe9efb91b5720fd464cb84e8d7c743d5aaf70ca",
        "byte_count": 29320,
    },
}

CANONICAL_RESULT_SHA256 = (
    "582a00acc9779f0de4d1be2ecbe9efb91b5720fd464cb84e8d7c743d5aaf70ca"
)

RESULT_FROZEN_V3_AUTHORITIES = {
    "objective_contract": {
        "path": "tests/reference/r2_joint_fit_objective_contract_v1.json",
        "sha256": "4e15ed2c7cca709b5bd54f04faa6cd2db50347e31dacc9676213b3fe375679e7",
        "byte_count": 15830,
    },
    "synthetic_fixture": {
        "path": "tests/reference/r2_joint_fit_synthetic_fixture_v1.json",
        "sha256": "6cef95b7fd3d56eee4a009214a5205d539c69681712e3088be6b7418a1b2199c",
        "byte_count": 6652,
    },
    "objective_validator": {
        "path": "tools/validate_r2_joint_fit_objective.py",
        "sha256": "7e04bc983e3f3312b0fed1afb3433ca83fcdb4d93e88fd0c73fa526d1ce5bb50",
        "byte_count": 35024,
    },
    "objective_test": {
        "path": "tests/test_r2_joint_fit_objective_contract.py",
        "sha256": "ec0005ccba7715dd00d6eddeaa9b5297ef9a2c3bea9a446c87cfff37fab5c3a8",
        "byte_count": 13622,
    },
    "objective_documentation": {
        "path": "docs/r2_joint_fit_objective_contract.md",
        "sha256": "cece29648591d19f498608dee9bf0f3510d82199f710ec506b56a4cf77197538",
        "byte_count": 6529,
    },
}

RESULT_CLAIMS = {
    "synthetic_optimizer_executed": True,
    "physical_optimizer_executed": False,
    "physical_fit_executed": False,
    "parameter_vector_publishable": False,
    "synthetic_identifiability_is_physical_evidence": False,
    "physical_uncertainty_claimed": False,
    "runner_equivalent_to_rfam_or_cantera": False,
    "hem_chemistry_enabled": False,
    "fluent_execution_approved": False,
    "reacting_hot_start_approved": False,
    "production_ready": False,
}

RESULT_SCOPE_LIMITS = [
    "synthetic affine fixture only; no RFAM or Cantera kinetics executed",
    "no detailed mechanism or experimental training or holdout data selected",
    "covariance is a deterministic software diagnostic, not physical uncertainty",
    "parameter and predictive intervals remain undefined and unpublished",
    "no HEM Fluent hot-start or production capability change",
]

RESULT_PARAMETER_ORDER = [
    "R1_A_SI",
    "R1_Ea_J_per_mol",
    "R1_order_C3H8",
    "R1_order_O2",
    "R2_A_SI",
    "R2_Ea_J_per_mol",
]

RESULT_REQUEST_FIELDS = [
    "ordered-parameter-vector",
    "initial-temperature-k",
    "initial-pressure-pa",
    "equivalence-ratio",
    "diluent-mass-fraction",
    "ordered-initial-mass-fraction-vector",
    "ordered-sample-times-s",
    "energy-basis-id",
]

RESULT_UNCERTAINTY_OUTPUTS = {
    "status": "not-publishable-synthetic-covariance-diagnostic-only",
    "marginal_and_joint_parameter_intervals": (
        "not-computed-no-frozen-coverage-distribution-rule"
    ),
    "training_and_holdout_predictive_intervals": (
        "not-computed-no-frozen-observation-uncertainty-model"
    ),
    "domain_of_applicability": "not-established-by-synthetic-fixture",
    "physical_evidence": False,
}

RESULT_ROOT_KEYS = {
    "schema", "status", "watermark", "authorities", "runner", "optimizer",
    "holdout", "identifiability", "claims", "scope_limits",
}
RESULT_SECTION_KEYS = {
    "runner": {
        "adapter_id", "affine_formula", "parameter_order", "required_case_request_fields",
        "initial_state_authority", "case_requests_per_training_call",
        "case_requests_per_holdout_call", "anchor_normalized_coordinates",
        "selection_training_calls", "holdout_calls_during_selection",
        "holdout_calls_after_selection", "identifiability_training_calls", "total_calls",
    },
    "optimizer": {
        "algorithm", "start_count", "selected_start_index", "selected_start_kind",
        "selected_normalized_coordinates", "selected_transformed_parameters",
        "selected_parameter_vector", "training_objective", "training_case_objectives",
        "global_selection", "per_start",
    },
    "holdout": {
        "evaluation_count", "objective", "case_objectives", "used_during_selection",
    },
    "identifiability": {
        "status", "coordinate_basis", "row_count", "parameter_count", "difference_modes",
        "runner_calls", "sensitivity_matrix", "gram_matrix", "jacobi_sweeps",
        "singular_values", "rank", "required_rank", "condition_number",
        "maximum_condition_number", "residual_sum_squares", "residual_variance",
        "covariance_matrix", "correlation_matrix",
        "maximum_absolute_offdiagonal_correlation", "rank_condition_correlation_pass",
        "uncertainty_outputs",
    },
}
RESULT_PER_START_KEYS = {
    "start_index", "start_kind", "start_normalized_coordinates",
    "selected_normalized_coordinates", "objective", "valid", "runner_evaluations",
    "cache_hits", "out_of_bound_skips", "completed_sweeps", "step_contractions",
    "accepted_coordinates", "termination", "trace_event_count", "trace_sha256",
}


class PackageError(RuntimeError):
    pass


def _strict_pairs(pairs: Sequence[Tuple[str, Any]]) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise PackageError("duplicate-json-key:" + key)
        result[key] = value
    return result


def _reject_constant(value: str) -> None:
    raise PackageError("nonfinite-json-number:" + value)


def _load_json(raw: bytes, label: str) -> Dict[str, Any]:
    try:
        value = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=_strict_pairs,
            parse_constant=_reject_constant,
        )
    except PackageError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError, ValueError) as exc:
        raise PackageError("invalid-json:" + label) from exc
    if not isinstance(value, dict):
        raise PackageError("json-root-not-object:" + label)
    return value


def _sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _resolve(root: Path, relative: str, label: str) -> Path:
    if not relative or "\\" in relative:
        raise PackageError("invalid-relative-path:" + label)
    path = (root / relative).resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise PackageError("path-escapes-root:" + label) from exc
    return path


def _read_authorities(root: Path) -> Dict[str, Tuple[Path, bytes]]:
    payloads: Dict[str, Tuple[Path, bytes]] = {}
    for name, authority in AUTHORITIES.items():
        path = _resolve(root, authority["path"], name + ".path")
        try:
            raw = path.read_bytes()
        except OSError as exc:
            raise PackageError("artifact-read-failed:" + name) from exc
        if len(raw) != authority["byte_count"]:
            raise PackageError("artifact-byte-count-mismatch:" + name)
        if _sha256(raw) != authority["sha256"]:
            raise PackageError("artifact-sha256-mismatch:" + name)
        payloads[name] = (path, raw)
    return payloads


def _all_finite(value: Any) -> bool:
    if isinstance(value, float):
        return math.isfinite(value)
    if isinstance(value, list):
        return all(_all_finite(item) for item in value)
    if isinstance(value, dict):
        return all(_all_finite(item) for item in value.values())
    return True


def _exact(actual: Any, expected: Any, label: str) -> None:
    if type(actual) is not type(expected):
        raise PackageError("result-contract-mismatch:" + label)
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise PackageError("result-contract-mismatch:" + label)
        for key in expected:
            _exact(actual[key], expected[key], label + "." + key)
        return
    if isinstance(expected, (list, tuple)):
        if len(actual) != len(expected):
            raise PackageError("result-contract-mismatch:" + label)
        for index, (actual_item, expected_item) in enumerate(zip(actual, expected)):
            _exact(actual_item, expected_item, label + "[" + str(index) + "]")
        return
    if actual != expected:
        raise PackageError("result-contract-mismatch:" + label)


def _exact_keys(value: Any, expected: set, label: str) -> Dict[str, Any]:
    if not isinstance(value, dict):
        raise PackageError("result-object-required:" + label)
    if set(value) != expected:
        raise PackageError("result-keys-mismatch:" + label)
    return value


def _exact_int(value: Any, expected: int, label: str) -> None:
    if type(value) is not int or value != expected:
        raise PackageError("result-integer-mismatch:" + label)


def _finite_float(value: Any, label: str) -> float:
    if type(value) is not float or not math.isfinite(value):
        raise PackageError("result-finite-float-required:" + label)
    return value


def _finite_vector(value: Any, size: int, label: str) -> None:
    if not isinstance(value, list) or len(value) != size:
        raise PackageError("result-vector-shape-mismatch:" + label)
    for index, item in enumerate(value):
        _finite_float(item, label + "[" + str(index) + "]")


def _finite_matrix(value: Any, rows: int, columns: int, label: str) -> None:
    if not isinstance(value, list) or len(value) != rows:
        raise PackageError("result-matrix-shape-mismatch:" + label)
    for index, row in enumerate(value):
        _finite_vector(row, columns, label + "[" + str(index) + "]")


def _verify_result(document: Dict[str, Any]) -> None:
    _exact_keys(document, RESULT_ROOT_KEYS, "root")
    _exact(
        document["schema"],
        "rfam.r2-joint-fit-synthetic-execution-result/v1",
        "schema",
    )
    _exact(
        document["status"],
        "blocked-synthetic-optimizer-execution-only",
        "status",
    )
    _exact(document["watermark"], "SYNTHETIC_NONPHYSICAL_TEST_ONLY", "watermark")
    if not _all_finite(document):
        raise PackageError("result-nonfinite")
    authorities = _exact_keys(
        document["authorities"], {"execution_contract", "frozen_v3"}, "authorities"
    )
    _exact(authorities["execution_contract"], AUTHORITIES["execution_contract"], "authorities.execution_contract")
    _exact(authorities["frozen_v3"], RESULT_FROZEN_V3_AUTHORITIES, "authorities.frozen_v3")

    runner = _exact_keys(document["runner"], RESULT_SECTION_KEYS["runner"], "runner")
    optimizer = _exact_keys(
        document["optimizer"], RESULT_SECTION_KEYS["optimizer"], "optimizer"
    )
    holdout = _exact_keys(document["holdout"], RESULT_SECTION_KEYS["holdout"], "holdout")
    identifiability = _exact_keys(
        document["identifiability"],
        RESULT_SECTION_KEYS["identifiability"],
        "identifiability",
    )
    _exact(document["claims"], RESULT_CLAIMS, "claims")
    _exact(document["scope_limits"], RESULT_SCOPE_LIMITS, "scope_limits")

    _exact(runner["adapter_id"], "fixture-anchored-affine-response-v1", "runner.adapter_id")
    _exact(
        runner["affine_formula"],
        "anchor_prediction+slope*(u-anchor_u)",
        "runner.affine_formula",
    )
    _exact(runner["parameter_order"], RESULT_PARAMETER_ORDER, "runner.parameter_order")
    _exact(
        runner["required_case_request_fields"],
        RESULT_REQUEST_FIELDS,
        "runner.required_case_request_fields",
    )
    _exact(
        runner["initial_state_authority"],
        "ordered-initial-mass-fraction-vector-no-reconstruction",
        "runner.initial_state_authority",
    )
    for key, expected in (
        ("case_requests_per_training_call", 2),
        ("case_requests_per_holdout_call", 2),
        ("selection_training_calls", 3284),
        ("holdout_calls_during_selection", 0),
        ("holdout_calls_after_selection", 1),
        ("identifiability_training_calls", 13),
        ("total_calls", 3298),
    ):
        _exact_int(runner[key], expected, "runner." + key)
    _finite_vector(runner["anchor_normalized_coordinates"], 6, "runner.anchor")

    _exact(optimizer["algorithm"], "bounded-coordinate-pattern-search-v1", "optimizer.algorithm")
    _exact_int(optimizer["start_count"], 17, "optimizer.start_count")
    _exact_int(optimizer["selected_start_index"], 12, "optimizer.selected_start_index")
    _exact(optimizer["selected_start_kind"], "halton-12", "optimizer.selected_start_kind")
    _finite_vector(optimizer["selected_normalized_coordinates"], 6, "optimizer.selected_u")
    _finite_vector(optimizer["selected_transformed_parameters"], 6, "optimizer.selected_transformed")
    _finite_vector(optimizer["selected_parameter_vector"], 6, "optimizer.selected_parameters")
    _finite_float(optimizer["training_objective"], "optimizer.training_objective")
    training_cases = _exact_keys(
        optimizer["training_case_objectives"],
        {"SYNTHETIC_TRAIN_001", "SYNTHETIC_TRAIN_002"},
        "optimizer.training_case_objectives",
    )
    for key, value in training_cases.items():
        _finite_float(value, "optimizer.training_case_objectives." + key)
    _exact(
        optimizer["global_selection"],
        "after-all-starts-choose-lowest-valid-objective-then-lexicographically-smallest-transformed-vector",
        "optimizer.global_selection",
    )
    starts = optimizer["per_start"]
    if not isinstance(starts, list) or len(starts) != 17:
        raise PackageError("result-per-start-count-mismatch")
    for index, value in enumerate(starts):
        start = _exact_keys(value, RESULT_PER_START_KEYS, "optimizer.per_start")
        _exact_int(start["start_index"], index, "optimizer.per_start.start_index")
        _exact(
            start["start_kind"],
            "midpoint" if index == 0 else "halton-" + str(index),
            "optimizer.per_start.start_kind",
        )
        _finite_vector(start["start_normalized_coordinates"], 6, "optimizer.per_start.start_u")
        _finite_vector(start["selected_normalized_coordinates"], 6, "optimizer.per_start.selected_u")
        _finite_float(start["objective"], "optimizer.per_start.objective")
        _exact(start["valid"], True, "optimizer.per_start.valid")
        for key in (
            "runner_evaluations", "cache_hits", "out_of_bound_skips", "completed_sweeps",
            "step_contractions", "accepted_coordinates", "trace_event_count",
        ):
            if type(start[key]) is not int or start[key] < 0:
                raise PackageError("result-nonnegative-integer-required:optimizer.per_start." + key)
        if not isinstance(start["termination"], str) or not start["termination"]:
            raise PackageError("result-string-required:optimizer.per_start.termination")
        trace_sha = start["trace_sha256"]
        if not isinstance(trace_sha, str) or len(trace_sha) != 64 or any(
            character not in "0123456789abcdef" for character in trace_sha
        ):
            raise PackageError("result-sha256-required:optimizer.per_start.trace_sha256")

    _exact_int(holdout["evaluation_count"], 1, "holdout.evaluation_count")
    _finite_float(holdout["objective"], "holdout.objective")
    holdout_cases = _exact_keys(
        holdout["case_objectives"],
        {"SYNTHETIC_HOLDOUT_001", "SYNTHETIC_HOLDOUT_002"},
        "holdout.case_objectives",
    )
    for key, value in holdout_cases.items():
        _finite_float(value, "holdout.case_objectives." + key)
    _exact(holdout["used_during_selection"], False, "holdout.used_during_selection")

    _exact(identifiability["status"], "synthetic-full-rank", "identifiability.status")
    _exact(
        identifiability["coordinate_basis"],
        "normalized-transformed-parameter-u",
        "identifiability.coordinate_basis",
    )
    for key, expected in (
        ("row_count", 58), ("parameter_count", 6), ("runner_calls", 13),
        ("rank", 6), ("required_rank", 6),
    ):
        _exact_int(identifiability[key], expected, "identifiability." + key)
    _exact(
        identifiability["difference_modes"],
        ["central"] * 6,
        "identifiability.difference_modes",
    )
    _finite_matrix(identifiability["sensitivity_matrix"], 58, 6, "identifiability.sensitivity")
    _finite_matrix(identifiability["gram_matrix"], 6, 6, "identifiability.gram")
    if type(identifiability["jacobi_sweeps"]) is not int or identifiability["jacobi_sweeps"] < 0:
        raise PackageError("result-nonnegative-integer-required:identifiability.jacobi_sweeps")
    _finite_vector(identifiability["singular_values"], 6, "identifiability.singular_values")
    for key in (
        "condition_number", "maximum_condition_number", "residual_sum_squares",
        "residual_variance", "maximum_absolute_offdiagonal_correlation",
    ):
        _finite_float(identifiability[key], "identifiability." + key)
    _finite_matrix(identifiability["covariance_matrix"], 6, 6, "identifiability.covariance")
    _finite_matrix(identifiability["correlation_matrix"], 6, 6, "identifiability.correlation")
    _exact(
        identifiability["rank_condition_correlation_pass"],
        True,
        "identifiability.rank_condition_correlation_pass",
    )
    _exact(
        identifiability["uncertainty_outputs"],
        RESULT_UNCERTAINTY_OUTPUTS,
        "identifiability.uncertainty_outputs",
    )

    try:
        canonical = (
            json.dumps(document, indent=2, ensure_ascii=True, allow_nan=False) + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise PackageError("result-canonical-serialization-failed") from exc
    if _sha256(canonical) != CANONICAL_RESULT_SHA256:
        raise PackageError("result-canonical-identity-mismatch")


def validate(root: Path) -> Tuple[int, str]:
    root = root.resolve()
    if not root.is_dir():
        raise PackageError("root-not-directory")

    # All package authorities are authenticated before the implementation is run.
    payloads = _read_authorities(root)
    expected_raw = payloads["deterministic_result"][1]
    expected_document = _load_json(expected_raw, "deterministic-result")
    _verify_result(expected_document)

    implementation = payloads["implementation"][0]
    with tempfile.TemporaryDirectory(prefix="rfam-r2-synthetic-package-") as directory:
        generated = Path(directory) / "result.json"
        environment = os.environ.copy()
        environment["PYTHONHASHSEED"] = "0"
        try:
            completed = subprocess.run(
                [
                    sys.executable,
                    str(implementation),
                    "--root",
                    str(root),
                    "--output",
                    str(generated),
                ],
                check=False,
                text=True,
                capture_output=True,
                timeout=120,
                env=environment,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise PackageError("implementation-execution-failed") from exc
        combined = completed.stdout + completed.stderr
        if completed.returncode != 0:
            raise PackageError("implementation-returned-nonzero")
        if "R2_JOINT_FIT_SYNTHETIC_EXECUTION_BLOCKED" not in completed.stdout:
            raise PackageError("implementation-blocked-marker-missing")
        if "R2_JOINT_FIT_SYNTHETIC_EXECUTION_FAIL" in combined:
            raise PackageError("implementation-failure-marker-present")
        try:
            actual_raw = generated.read_bytes()
        except OSError as exc:
            raise PackageError("implementation-result-missing") from exc
    if actual_raw != expected_raw:
        raise PackageError("deterministic-result-byte-mismatch")
    actual_document = _load_json(actual_raw, "regenerated-result")
    _verify_result(actual_document)
    return len(AUTHORITIES) * 3 + 18, _sha256(actual_raw)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    try:
        checks, result_hash = validate(args.root)
    except PackageError as exc:
        print(FAIL_MARKER + " " + str(exc), file=sys.stderr)
        return 2
    except Exception as exc:  # Every malformed package fails without traceback.
        print(FAIL_MARKER + " unexpected:" + type(exc).__name__, file=sys.stderr)
        return 2
    print(
        BLOCKED_MARKER
        + " synthetic_only=1 authorities=" + str(len(AUTHORITIES))
        + " starts=17 selection_calls=3284 holdout_calls=1 identifiability_rank=6"
        + " physical_fit_executed=0 production_ready=0"
        + " result_sha256=" + result_hash
        + " checks=" + str(checks)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
