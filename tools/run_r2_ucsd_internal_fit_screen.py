#!/usr/bin/env python3
"""Run a blocked, internal-only R1/R2 route screen against UCSD trajectories."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
import os
import platform
import subprocess
import sys
import tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Sequence

class ScreenError(RuntimeError):
    pass


EXPECTED_PARAMETER_NAMES = [
    "R1_A_SI",
    "R1_Ea_J_per_mol",
    "R1_order_C3H8",
    "R1_order_O2",
    "R2_A_SI",
    "R2_Ea_J_per_mol",
]
DEFAULT_PARAMETERS = [1.0e10, 125520.0, 0.1, 1.65, 2.0e12, 167360.0]
EXPECTED_CONTRACT_SHA256 = (
    "8cdf5bb9b11fccbaf4aafef07840bc0c5cbe22a3b99d9ff685eef5e62d703796"
)
EXPECTED_EXECUTION = {
    "platform_system": "Linux",
    "platform_machine": "x86_64",
    "python_version": "3.8.10",
    "runner_timeout_seconds": 30,
}
EXPECTED_ARTIFACTS = {
    "source_reference_sha256": (
        "ffcf5fb0c88b1cf2a4be465bfe2b53b9b306012154fc74c02b45a2b5a1c3984f"
    ),
    "training_reference": {
        "logical_path": "internal://a408-r2-ucsd/training-v1.json",
        "sha256": "bcc2a6251384734f44d5b31cbeb10ebd5a77ae265eacdb0e1aaba82040c6d71d",
    },
    "holdout_reference": {
        "logical_path": "internal://a408-r2-ucsd/holdout-v1.json",
        "sha256": "856f86914e61607c4a5f3113b72fbb79d3a68562162a2b0bf43a29b4c3865def",
    },
    "parameterized_runner": {
        "logical_path": "build://rfam-r2-0d-fit-runner/linux-x86_64-v1",
        "sha256": "d080d0ed5a17fb029fc5e01fbad5dbe6f728f88f0686f0292606a8e81344f9ad",
    },
    "repository_artifacts": [
        {
            "path": "tools/r2_0d_rfam_fit_runner.cpp",
            "sha256": "ee0accda3f6d414f0757ec3723c24f93a0da7560714ab866682ea9fecf2c2c73",
        },
        {
            "path": "tools/validate_r2_parameterized_runner.py",
            "sha256": "0a1c89f10f265ccf23c753eaeefec32fcb4c86864ca09504cd8fbd3672fe0271",
        },
        {
            "path": "tools/split_r2_ucsd_internal_reference.py",
            "sha256": "0372b833453c8d8c479238a426f17bad5564ddaa18765a071d13b09baa421985",
        },
        {
            "path": "tools/run_r2_0d_rfam_protocol.py",
            "sha256": "5a4642c37ce4b07ca97db82a6af218d39e039b8a73e417d20efb6094d43f58a8",
        },
        {
            "path": "tools/run_r2_joint_fit_synthetic.py",
            "sha256": "0669e6d7c2d81fc42a71dfd96688fdb13d1c821dbe953c8a6f1f9bbefd40f99d",
        },
        {
            "path": "tools/screen_r2_detailed_candidate.py",
            "sha256": "6e88e92ebe84db8ae378f95399de4764ca55f122ef9156ffb5c42e48874a14e7",
        },
    ],
}
EXPECTED_PARAMETERS = [
    ("R1_A_SI", "log10", 100000.0, 1.0e20),
    ("R1_Ea_J_per_mol", "linear", 0.0, 300000.0),
    ("R1_order_C3H8", "linear", 0.0, 2.0),
    ("R1_order_O2", "linear", 0.0, 2.5),
    ("R2_A_SI", "log10", 100000.0, 1.0e20),
    ("R2_Ea_J_per_mol", "linear", 0.0, 300000.0),
]
EXPECTED_OBJECTIVE_POLICY = {
    "observable_order": [
        "ignition_delay", "qdot", "temperature", "C3H8", "O2", "CO", "CO2", "H2O"
    ],
    "ignition_delay_definition": (
        "first-time-temperature-minus-initial-temperature-reaches-400-K"
    ),
    "ignition_delay_residual": "natural-log-ratio-divided-by-sigma-log",
    "history_residual": "prediction-minus-reference-divided-by-fixed-scale",
    "case_aggregation": "arithmetic-mean-in-frozen-case-order",
    "reduction": "python-math-fsum-binary64-in-case-observable-sample-order",
    "nonignition_policy": "invalid-candidate-never-replaces-incumbent",
    "normalization_scales": {
        "ignition_delay_sigma_log": 0.1,
        "qdot_w_m3": 100000000.0,
        "temperature_k": 20.0,
        "C3H8_mass_fraction": 0.01,
        "O2_mass_fraction": 0.01,
        "CO_mass_fraction": 0.01,
        "CO2_mass_fraction": 0.01,
        "H2O_mass_fraction": 0.01,
    },
    "component_weights": {
        "ignition_delay": 0.2,
        "qdot": 0.2,
        "temperature": 0.15,
        "C3H8": 0.1,
        "O2": 0.1,
        "CO": 0.1,
        "CO2": 0.075,
        "H2O": 0.075,
    },
}
EXPECTED_INTERPRETATION = [
    "optimizer output is an internal model-structure diagnostic only",
    "protocol training and holdout labels are not experimental validation partitions",
    "a low objective cannot cure missing radical species or grant mechanism rights",
    "holdout data cannot be read by the objective during parameter selection",
    "no parameter may enter runtime defaults or a production manifest",
]

# These globals are populated only after their source bytes pass the contract hash gate.
ObjectiveEvaluation: Any = None
ParameterSpec: Any = None
SearchConfig: Any = None
SyntheticExecutionError: Any = None
SPECIES: list[str] = []
parse_runner_output: Any = None
radical_inverse: Any = None
run_pattern_search: Any = None
validate_reference: Any = None


def strict_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ScreenError("duplicate-json-key:" + key)
        result[key] = value
    return result


def reject_constant(token: str) -> None:
    raise ScreenError("nonfinite-json-token:" + token)


def load_json(path: Path, label: str) -> tuple[dict[str, Any], bytes]:
    try:
        raw = path.read_bytes()
        document = json.loads(
            raw.decode("ascii"),
            object_pairs_hook=strict_pairs,
            parse_constant=reject_constant,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ScreenError(label + "-read-failed:" + str(error)) from error
    if not isinstance(document, dict):
        raise ScreenError(label + "-root-not-object")
    return document, raw


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def file_digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def verify_runtime() -> None:
    actual = {
        "platform_system": platform.system(),
        "platform_machine": platform.machine(),
        "python_version": platform.python_version(),
        "runner_timeout_seconds": 30,
    }
    if actual != EXPECTED_EXECUTION:
        raise ScreenError("execution-runtime-mismatch")


def verify_repository_artifacts(repository_root: Path) -> None:
    for entry in EXPECTED_ARTIFACTS["repository_artifacts"]:
        path = repository_root / entry["path"]
        if not path.is_file() or file_digest(path) != entry["sha256"]:
            raise ScreenError("repository-artifact-mismatch:" + entry["path"])


def load_bound_helpers() -> None:
    global ObjectiveEvaluation
    global ParameterSpec
    global SearchConfig
    global SyntheticExecutionError
    global SPECIES
    global parse_runner_output
    global radical_inverse
    global run_pattern_search
    global validate_reference

    protocol = importlib.import_module("run_r2_0d_rfam_protocol")
    optimizer = importlib.import_module("run_r2_joint_fit_synthetic")
    candidate = importlib.import_module("screen_r2_detailed_candidate")
    ObjectiveEvaluation = optimizer.ObjectiveEvaluation
    ParameterSpec = optimizer.ParameterSpec
    SearchConfig = optimizer.SearchConfig
    SyntheticExecutionError = optimizer.SyntheticExecutionError
    SPECIES = list(candidate.SPECIES)
    parse_runner_output = protocol.parse_runner_output
    radical_inverse = optimizer.radical_inverse
    run_pattern_search = optimizer.run_pattern_search
    validate_reference = candidate.validate_reference


def validate_partition_reference(
    document: dict[str, Any],
    expected_partition: str,
    expected_case_ids: Sequence[str],
) -> tuple[list[float], list[dict[str, Any]]]:
    expected_keys = {
        "cantera_version", "cases", "generator_sha256", "major_species",
        "mechanism_sha256", "partition_artifact", "reactor",
        "runtime_identity", "sample_times_s", "schema",
        "source_reference_sha256", "status", "units",
    }
    if set(document) != expected_keys:
        raise ScreenError("partition-reference-key-schema-mismatch:" + expected_partition)
    if (
        document["partition_artifact"] != expected_partition
        or document["source_reference_sha256"]
        != EXPECTED_ARTIFACTS["source_reference_sha256"]
    ):
        raise ScreenError("partition-reference-origin-mismatch:" + expected_partition)
    sample_times, cases = validate_reference(document)
    actual_case_ids = [case["id"] for case in cases]
    if actual_case_ids != list(expected_case_ids) or any(
        case.get("partition") != expected_partition for case in cases
    ):
        raise ScreenError("partition-reference-cases-mismatch:" + expected_partition)
    return sample_times, cases


def strip_fit_parameter_header(
    stdout: str,
    parameters: Sequence[float],
    case_id: str,
) -> str:
    lines = stdout.splitlines(keepends=True)
    expected = ["fit_parameter_schema=rfam.r2-joint-fit-parameters/v1\n"]
    expected.extend(
        "fit_" + name + "=" + format(float(value), ".17g") + "\n"
        for name, value in zip(EXPECTED_PARAMETER_NAMES, parameters)
    )
    if len(lines) < len(expected):
        raise ScreenError("runner-fit-header-truncated:" + case_id)
    for index, (actual, wanted) in enumerate(zip(lines, expected)):
        if actual != wanted:
            raise ScreenError(
                "runner-fit-header-mismatch:" + case_id + ":" + str(index)
            )
    return "".join(lines[len(expected):])


def require_sha(value: str, label: str) -> None:
    if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
        raise ScreenError("invalid-sha256:" + label)


def exact_keys(value: dict[str, Any], expected: Sequence[str], label: str) -> None:
    if set(value) != set(expected):
        raise ScreenError("key-schema-mismatch:" + label)


def object_value(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ScreenError("not-object:" + label)
    return value


def list_value(value: Any, label: str) -> list[Any]:
    if not isinstance(value, list):
        raise ScreenError("not-list:" + label)
    return value


def finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ScreenError("not-number:" + label)
    result = float(value)
    if not math.isfinite(result):
        raise ScreenError("nonfinite:" + label)
    return result


def validate_contract(document: dict[str, Any]) -> tuple[
    tuple[ParameterSpec, ...], SearchConfig, list[tuple[str, tuple[float, ...]]]
]:
    exact_keys(
        document,
        [
            "schema", "status", "execution", "artifacts", "case_partition",
            "parameters", "objective", "optimizer", "scope",
            "required_interpretation",
        ],
        "root",
    )
    if document["schema"] != "rfam.r2-ucsd-internal-fit-screen-contract/v1":
        raise ScreenError("contract-schema-mismatch")
    if document["status"] != "blocked-internal-route-screen-only":
        raise ScreenError("contract-status-mismatch")
    if document["execution"] != EXPECTED_EXECUTION:
        raise ScreenError("contract-execution-mismatch")
    if document["artifacts"] != EXPECTED_ARTIFACTS:
        raise ScreenError("contract-artifacts-mismatch")

    partition = object_value(document["case_partition"], "partition")
    exact_keys(
        partition,
        [
            "training_case_ids", "holdout_case_ids",
            "holdout_calls_during_selection", "holdout_calls_after_selection",
        ],
        "partition",
    )
    training = list_value(partition["training_case_ids"], "training-cases")
    holdout = list_value(partition["holdout_case_ids"], "holdout-cases")
    if (
        any(not isinstance(item, str) or not item for item in training + holdout)
        or len(set(training + holdout)) != len(training + holdout)
        or training
        != ["train_phi0p8_diluted", "train_phi1p0", "train_phi1p2"]
        or holdout
        != ["holdout_phi0p6_diluted", "holdout_phi1p1", "holdout_phi1p4"]
        or partition["holdout_calls_during_selection"] != 0
        or partition["holdout_calls_after_selection"] != 1
    ):
        raise ScreenError("partition-invalid")

    parameters = list_value(document["parameters"], "parameters")
    if len(parameters) != len(EXPECTED_PARAMETER_NAMES):
        raise ScreenError("parameter-count-mismatch")
    specs = []
    for index, (entry_raw, expected_parameter) in enumerate(
        zip(parameters, EXPECTED_PARAMETERS)
    ):
        expected_name, expected_transform, expected_lower, expected_upper = (
            expected_parameter
        )
        entry = object_value(entry_raw, "parameter")
        exact_keys(entry, ["name", "transform", "lower", "upper"], "parameter")
        if (
            entry["name"] != expected_name
            or entry["transform"] != expected_transform
            or entry["lower"] != expected_lower
            or entry["upper"] != expected_upper
        ):
            raise ScreenError("parameter-identity-mismatch:" + str(index))
        lower = finite(entry["lower"], expected_name + ".lower")
        upper = finite(entry["upper"], expected_name + ".upper")
        if not upper > lower or (entry["transform"] == "log10" and lower <= 0.0):
            raise ScreenError("parameter-bounds-invalid:" + expected_name)
        transformed_lower = math.log10(lower) if entry["transform"] == "log10" else lower
        transformed_upper = math.log10(upper) if entry["transform"] == "log10" else upper
        specs.append(ParameterSpec(
            expected_name, entry["transform"], lower, upper,
            transformed_lower, transformed_upper,
        ))

    objective = object_value(document["objective"], "objective")
    exact_keys(
        objective,
        [
            "observable_order", "ignition_delay_definition",
            "ignition_delay_residual", "history_residual", "case_aggregation",
            "reduction", "nonignition_policy", "normalization_scales",
            "component_weights",
        ],
        "objective",
    )
    if objective != EXPECTED_OBJECTIVE_POLICY:
        raise ScreenError("objective-policy-mismatch")
    scales = object_value(objective["normalization_scales"], "scales")
    weights = object_value(objective["component_weights"], "weights")
    expected_scales = {
        "ignition_delay_sigma_log", "qdot_w_m3", "temperature_k",
        "C3H8_mass_fraction", "O2_mass_fraction", "CO_mass_fraction",
        "CO2_mass_fraction", "H2O_mass_fraction",
    }
    expected_weights = {"ignition_delay", "qdot", "temperature", "C3H8", "O2", "CO", "CO2", "H2O"}
    if set(scales) != expected_scales or set(weights) != expected_weights:
        raise ScreenError("objective-scale-or-weight-schema-mismatch")
    if any(finite(value, "scale") <= 0.0 for value in scales.values()):
        raise ScreenError("objective-scale-invalid")
    if any(finite(value, "weight") <= 0.0 for value in weights.values()):
        raise ScreenError("objective-weight-invalid")
    if math.fsum(float(value) for value in weights.values()) != 1.0:
        raise ScreenError("objective-weights-do-not-sum-to-one")

    optimizer = object_value(document["optimizer"], "optimizer")
    exact_keys(
        optimizer,
        [
            "algorithm", "selection_split", "starts", "initial_normalized_step",
            "step_contraction", "minimum_normalized_step",
            "maximum_objective_evaluations_per_start", "parallel_case_processes",
            "global_selection",
        ],
        "optimizer",
    )
    expected_optimizer = {
        "algorithm": "bounded-coordinate-pattern-search-v1",
        "selection_split": "training-only",
        "starts": {
            "midpoint_first": True,
            "halton_indices": [1, 2, 3, 4],
            "halton_bases": [2, 3, 5, 7, 11, 13],
        },
        "initial_normalized_step": 0.25,
        "step_contraction": 0.5,
        "minimum_normalized_step": 0.015625,
        "maximum_objective_evaluations_per_start": 800,
        "parallel_case_processes": 3,
        "global_selection": (
            "lowest-training-objective-then-lexicographically-smallest-transformed-vector"
        ),
    }
    if optimizer != expected_optimizer:
        raise ScreenError("optimizer-policy-mismatch")
    starts = object_value(optimizer["starts"], "starts")
    exact_keys(starts, ["midpoint_first", "halton_indices", "halton_bases"], "starts")
    indices = list_value(starts["halton_indices"], "halton-indices")
    bases = list_value(starts["halton_bases"], "halton-bases")
    if starts["midpoint_first"] is not True or bases != [2, 3, 5, 7, 11, 13]:
        raise ScreenError("start-sequence-mismatch")
    if any(isinstance(index, bool) or not isinstance(index, int) or index <= 0 for index in indices):
        raise ScreenError("halton-index-invalid")
    start_points = [("midpoint", tuple(0.5 for _ in specs))]
    start_points.extend(
        (
            "halton-" + str(index),
            tuple(radical_inverse(index, int(base)) for base in bases),
        )
        for index in indices
    )
    config = SearchConfig(
        finite(optimizer["initial_normalized_step"], "initial-step"),
        finite(optimizer["step_contraction"], "contraction"),
        finite(optimizer["minimum_normalized_step"], "minimum-step"),
        int(optimizer["maximum_objective_evaluations_per_start"]),
        len(specs),
    )
    if (
        isinstance(optimizer["maximum_objective_evaluations_per_start"], bool)
        or not isinstance(optimizer["maximum_objective_evaluations_per_start"], int)
        or config.maximum_evaluations <= 0
        or not (0.0 < config.minimum_step <= config.initial_step <= 1.0)
        or not (0.0 < config.contraction < 1.0)
        or optimizer["parallel_case_processes"] != 3
    ):
        raise ScreenError("optimizer-numeric-policy-invalid")

    scope = object_value(document["scope"], "scope")
    exact_keys(
        scope,
        [
            "mechanism_role", "experimental_training_selected",
            "experimental_holdout_selected", "mechanism_rights_clear_for_fit",
            "parameter_vector_publishable", "physical_validation_approved",
            "production_ready", "hem_chemistry_enabled",
            "fluent_execution_approved",
        ],
        "scope",
    )
    for key, value in scope.items():
        if key == "mechanism_role":
            if value != "ucsd-converted-candidate-internal-route-screen-not-selected-for-fit":
                raise ScreenError("mechanism-role-mismatch")
        elif value is not False:
            raise ScreenError("scope-authorization-not-false:" + key)
    if document["required_interpretation"] != EXPECTED_INTERPRETATION:
        raise ScreenError("required-interpretation-mismatch")
    return tuple(specs), config, start_points


class R2PhysicalObjective:
    def __init__(
        self,
        runner: Path,
        specs: Sequence[ParameterSpec],
        training_cases: Sequence[dict[str, Any]],
        sample_times: Sequence[float],
        contract: dict[str, Any],
    ) -> None:
        self.runner = runner
        self.specs = tuple(specs)
        self.case_by_id = {case["id"]: case for case in training_cases}
        self.sample_times = list(sample_times)
        partition = contract["case_partition"]
        self.case_ids = {
            "training": tuple(partition["training_case_ids"]),
        }
        self.expected_holdout_case_ids = tuple(partition["holdout_case_ids"])
        self.runner_timeout_seconds = int(
            contract["execution"]["runner_timeout_seconds"]
        )
        objective = contract["objective"]
        self.scales = objective["normalization_scales"]
        self.weights = objective["component_weights"]
        self.selection_active = True
        self.objective_calls = Counter()
        self.case_process_attempts = Counter()
        self.case_process_successes = Counter()
        self.holdout_access_attempts_during_selection = 0
        self.invalid_reasons: Counter[str] = Counter()

    def install_holdout(
        self,
        holdout_cases: Sequence[dict[str, Any]],
        sample_times: Sequence[float],
    ) -> None:
        if self.selection_active:
            raise ScreenError("holdout-installed-during-selection")
        if "holdout" in self.case_ids or list(sample_times) != self.sample_times:
            raise ScreenError("holdout-installation-state-invalid")
        actual_case_ids = tuple(case["id"] for case in holdout_cases)
        if actual_case_ids != self.expected_holdout_case_ids:
            raise ScreenError("holdout-installation-case-mismatch")
        self.case_by_id.update((case["id"], case) for case in holdout_cases)
        self.case_ids["holdout"] = actual_case_ids

    def parameters(self, point: Sequence[float]) -> tuple[float, ...]:
        return tuple(spec.physical_from_u(float(value)) for spec, value in zip(self.specs, point))

    def _run_case(
        self,
        parameters: Sequence[float],
        case_id: str,
    ) -> tuple[dict[str, str], list[dict[str, Any]]]:
        case = self.case_by_id[case_id]
        initial = case["initial_major_species_mass_fractions"]
        command = [
            str(self.runner),
            *[format(float(value), ".17g") for value in parameters],
            case_id,
            format(float(case["temperature_k"]), ".17g"),
            format(float(case["pressure_pa"]), ".17g"),
            str(len(self.sample_times)),
            "0",
            *[format(float(initial[name]), ".17g") for name in SPECIES],
            *[format(float(value), ".17g") for value in self.sample_times],
        ]
        result = subprocess.run(
            command,
            check=False,
            text=True,
            capture_output=True,
            timeout=self.runner_timeout_seconds,
        )
        if result.returncode != 0:
            raise ScreenError("runner-rejected:" + case_id + ":rc=" + str(result.returncode))
        trajectory_stdout = strip_fit_parameter_header(
            result.stdout, parameters, case_id
        )
        headers, samples = parse_runner_output(
            trajectory_stdout, result.stderr, case_id, len(self.sample_times)
        )
        return headers, samples

    def evaluate(self, point: Sequence[float], split: str) -> ObjectiveEvaluation:
        if split == "holdout" and self.selection_active:
            self.holdout_access_attempts_during_selection += 1
            raise ScreenError("holdout-accessed-during-selection")
        if split not in self.case_ids:
            raise ScreenError("unknown-objective-split:" + split)
        self.objective_calls[split] += 1
        case_ids = self.case_ids[split]
        parameters = self.parameters(point)
        try:
            self.case_process_attempts[split] += len(case_ids)
            with ThreadPoolExecutor(max_workers=3) as executor:
                futures = [
                    executor.submit(self._run_case, parameters, case_id)
                    for case_id in case_ids
                ]
                predictions = [future.result() for future in futures]
            self.case_process_successes[split] += len(case_ids)
            residuals: list[float] = []
            case_objectives: dict[str, float] = {}
            case_count = len(case_ids)
            sample_count = len(self.sample_times)
            species_index = {name: index for index, name in enumerate(SPECIES)}
            for case_id, (headers, samples) in zip(case_ids, predictions):
                reference = self.case_by_id[case_id]
                if headers.get("ignition_detected") != "1":
                    raise ScreenError("ignition-not-detected:" + case_id)
                predicted_idt = float(headers["ignition_delay_s"])
                reference_idt = float(reference["ignition_delay_s"])
                if predicted_idt <= 0.0 or reference_idt <= 0.0:
                    raise ScreenError("ignition-delay-invalid:" + case_id)
                case_rows = [
                    math.sqrt(float(self.weights["ignition_delay"]) / case_count)
                    * math.log(predicted_idt / reference_idt)
                    / float(self.scales["ignition_delay_sigma_log"])
                ]
                reference_samples = reference["samples"]
                observable_values = {
                    "qdot": (
                        [float(sample["qdot_w_m3"]) for sample in samples],
                        [float(sample["heat_release_rate_w_m3"]) for sample in reference_samples],
                        float(self.scales["qdot_w_m3"]),
                    ),
                    "temperature": (
                        [float(sample["temperature_k"]) for sample in samples],
                        [float(sample["temperature_k"]) for sample in reference_samples],
                        float(self.scales["temperature_k"]),
                    ),
                }
                for name in ("C3H8", "O2", "CO", "CO2", "H2O"):
                    observable_values[name] = (
                        [
                            float(sample["mass_fractions"][species_index[name]])
                            for sample in samples
                        ],
                        [
                            float(sample["major_species_mass_fractions"][name])
                            for sample in reference_samples
                        ],
                        float(self.scales[name + "_mass_fraction"]),
                    )
                for observable in ("qdot", "temperature", "C3H8", "O2", "CO", "CO2", "H2O"):
                    predicted, target, scale = observable_values[observable]
                    factor = math.sqrt(
                        float(self.weights[observable]) / (case_count * sample_count)
                    )
                    case_rows.extend(
                        factor * (left - right) / scale
                        for left, right in zip(predicted, target)
                    )
                if any(not math.isfinite(value) for value in case_rows):
                    raise ScreenError("nonfinite-objective-row:" + case_id)
                residuals.extend(case_rows)
                case_objectives[case_id] = math.fsum(value * value for value in case_rows)
            objective = math.fsum(value * value for value in residuals)
            if not math.isfinite(objective):
                raise ScreenError("nonfinite-objective")
            return ObjectiveEvaluation(
                True, objective, tuple(residuals), case_objectives
            )
        except (
            ScreenError,
            ValueError,
            KeyError,
            IndexError,
            subprocess.TimeoutExpired,
        ) as error:
            reason = str(error)
            self.invalid_reasons[reason] += 1
            return ObjectiveEvaluation.invalid(reason)


def rounded(value: Any) -> Any:
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ScreenError("nonfinite-result")
        return float(format(value, ".17g"))
    if isinstance(value, dict):
        return {key: rounded(item) for key, item in value.items()}
    if isinstance(value, list):
        return [rounded(item) for item in value]
    return value


def canonical_bytes(document: dict[str, Any]) -> bytes:
    return (
        json.dumps(
            rounded(document), sort_keys=True, separators=(",", ":"),
            ensure_ascii=True, allow_nan=False,
        )
        + "\n"
    ).encode("ascii")


def write_atomic(path: Path, raw: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".tmp.", dir=str(path.parent)
    )
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise


def execute(
    contract_path: Path,
    training_reference_path: Path,
    holdout_reference_path: Path,
    runner_path: Path,
) -> dict[str, Any]:
    verify_runtime()
    repository_root = Path(__file__).resolve().parent.parent
    canonical_contract_path = (
        repository_root
        / "tests"
        / "reference"
        / "r2_ucsd_internal_fit_screen_contract_v1.json"
    ).resolve()
    if contract_path != canonical_contract_path:
        raise ScreenError("contract-path-mismatch")
    contract, contract_raw = load_json(contract_path, "contract")
    if digest(contract_raw) != EXPECTED_CONTRACT_SHA256:
        raise ScreenError("contract-sha256-mismatch")
    if contract.get("execution") != EXPECTED_EXECUTION:
        raise ScreenError("contract-execution-mismatch")
    if contract.get("artifacts") != EXPECTED_ARTIFACTS:
        raise ScreenError("contract-artifacts-mismatch")
    verify_repository_artifacts(repository_root)
    load_bound_helpers()

    training_reference, training_raw = load_json(
        training_reference_path, "training-reference"
    )
    if digest(training_raw) != EXPECTED_ARTIFACTS["training_reference"]["sha256"]:
        raise ScreenError("training-reference-sha256-mismatch")
    if (
        not runner_path.is_file()
        or file_digest(runner_path)
        != EXPECTED_ARTIFACTS["parameterized_runner"]["sha256"]
    ):
        raise ScreenError("runner-sha256-mismatch")

    specs, config, start_points = validate_contract(contract)
    sample_times, training_cases = validate_partition_reference(
        training_reference,
        "protocol-training",
        contract["case_partition"]["training_case_ids"],
    )

    objective = R2PhysicalObjective(
        runner_path, specs, training_cases, sample_times, contract
    )
    default_u = tuple(
        spec.u_from_physical(value) for spec, value in zip(specs, DEFAULT_PARAMETERS)
    )
    default_training = objective.evaluate(default_u, "training")
    if not default_training.valid:
        raise ScreenError("default-training-evaluation-invalid:" + default_training.reason)

    search_results = []
    for start_index, (kind, point) in enumerate(start_points):
        search_results.append(run_pattern_search(
            start_index,
            kind,
            point,
            specs,
            config,
            lambda candidate: objective.evaluate(candidate, "training"),
        ))
    valid_results = [result for result in search_results if result.valid]
    if not valid_results:
        raise ScreenError("all-search-starts-invalid")
    selected = min(
        valid_results,
        key=lambda result: (
            result.evaluation.objective,
            tuple(
                spec.transformed_from_u(value)
                for spec, value in zip(specs, result.selected_u)
            ),
        ),
    )
    if objective.holdout_access_attempts_during_selection != 0:
        raise ScreenError("holdout-access-accounting-mismatch")

    # The holdout artifact is intentionally unopened until selection is final.
    objective.selection_active = False
    holdout_reference, holdout_raw = load_json(
        holdout_reference_path, "holdout-reference"
    )
    if digest(holdout_raw) != EXPECTED_ARTIFACTS["holdout_reference"]["sha256"]:
        raise ScreenError("holdout-reference-sha256-mismatch")
    holdout_times, holdout_cases = validate_partition_reference(
        holdout_reference,
        "protocol-holdout",
        contract["case_partition"]["holdout_case_ids"],
    )
    metadata_keys = [
        "cantera_version", "generator_sha256", "major_species",
        "mechanism_sha256", "reactor", "runtime_identity", "sample_times_s",
        "schema", "source_reference_sha256", "status", "units",
    ]
    if any(
        training_reference[key] != holdout_reference[key]
        for key in metadata_keys
    ):
        raise ScreenError("partition-reference-metadata-mismatch")
    objective.install_holdout(holdout_cases, holdout_times)
    selected_holdout = objective.evaluate(selected.selected_u, "holdout")
    if not selected_holdout.valid:
        raise ScreenError("selected-holdout-invalid:" + selected_holdout.reason)
    if objective.objective_calls["holdout"] != 1:
        raise ScreenError("holdout-call-count-mismatch")
    expected_training_calls = 1 + sum(
        result.evaluation_count for result in search_results
    )
    if (
        objective.objective_calls["training"] != expected_training_calls
        or objective.case_process_attempts["training"]
        != 3 * expected_training_calls
        or objective.case_process_attempts["holdout"] != 3
    ):
        raise ScreenError("runner-accounting-mismatch")

    return {
        "schema": "rfam.r2-ucsd-internal-fit-screen-result/v1",
        "status": "blocked-internal-route-screen-only",
        "source": {
            "contract_path": (
                "tests/reference/r2_ucsd_internal_fit_screen_contract_v1.json"
            ),
            "contract_sha256": EXPECTED_CONTRACT_SHA256,
            "training_reference": EXPECTED_ARTIFACTS["training_reference"],
            "holdout_reference": EXPECTED_ARTIFACTS["holdout_reference"],
            "parameterized_runner": EXPECTED_ARTIFACTS["parameterized_runner"],
            "repository_artifacts": EXPECTED_ARTIFACTS["repository_artifacts"],
            "execution": EXPECTED_EXECUTION,
        },
        "baseline": {
            "parameter_vector": DEFAULT_PARAMETERS,
            "training_objective": default_training.objective,
            "training_case_objectives": default_training.case_objectives,
        },
        "selection": {
            "start_count": len(search_results),
            "valid_start_count": len(valid_results),
            "selected_start_index": selected.start_index,
            "selected_start_kind": selected.start_kind,
            "selected_normalized_coordinates": list(selected.selected_u),
            "selected_parameter_vector": list(objective.parameters(selected.selected_u)),
            "training_objective": selected.evaluation.objective,
            "training_case_objectives": selected.evaluation.case_objectives,
            "holdout_objective": selected_holdout.objective,
            "holdout_case_objectives": selected_holdout.case_objectives,
        },
        "starts": [
            {
                "index": result.start_index,
                "kind": result.start_kind,
                "valid": result.valid,
                "objective": result.evaluation.objective if result.valid else None,
                "selected_normalized_coordinates": list(result.selected_u),
                "selected_parameter_vector": list(
                    objective.parameters(result.selected_u)
                ),
                "evaluation_count": result.evaluation_count,
                "cache_hits": result.cache_hits,
                "accepted_coordinates": result.accepted_coordinates,
                "completed_sweeps": result.completed_sweeps,
                "contractions": result.contractions,
                "termination": result.termination,
                "trace_sha256": result.trace_sha256,
            }
            for result in search_results
        ],
        "runner_accounting": {
            "training_objective_calls": objective.objective_calls["training"],
            "holdout_access_attempts_during_selection": (
                objective.holdout_access_attempts_during_selection
            ),
            "holdout_objective_calls_after_selection": objective.objective_calls["holdout"],
            "training_case_process_attempts": objective.case_process_attempts["training"],
            "training_case_process_successes": objective.case_process_successes["training"],
            "holdout_case_process_attempts": objective.case_process_attempts["holdout"],
            "holdout_case_process_successes": objective.case_process_successes["holdout"],
            "invalid_candidate_reasons": dict(sorted(objective.invalid_reasons.items())),
        },
        "claims": {
            "parameterized_runner_default_byte_equivalence_is_separate_required_evidence": True,
            "optimizer_executed_against_internal_detailed_candidate": True,
            "experimental_fit_executed": False,
            "physical_validation_claimed": False,
            "parameter_vector_publishable": False,
            "production_ready": False,
        },
        "limits": list(contract["required_interpretation"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--training-reference", type=Path, required=True)
    parser.add_argument("--holdout-reference", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    contract_path = arguments.contract.resolve()
    training_reference_path = arguments.training_reference.resolve()
    holdout_reference_path = arguments.holdout_reference.resolve()
    runner_path = arguments.runner.resolve()
    output_path = arguments.output.resolve()
    try:
        if output_path in {
            contract_path,
            training_reference_path,
            holdout_reference_path,
            runner_path,
        }:
            raise ScreenError("output-aliases-input")
        if output_path.exists():
            if not output_path.is_file():
                raise ScreenError("output-path-not-file")
            output_path.unlink()
        result = execute(
            contract_path,
            training_reference_path,
            holdout_reference_path,
            runner_path,
        )
        raw = canonical_bytes(result)
        write_atomic(output_path, raw)
        print(
            "R2_UCSD_INTERNAL_FIT_SCREEN_BLOCKED training_objective="
            + format(result["selection"]["training_objective"], ".12g")
            + " holdout_objective="
            + format(result["selection"]["holdout_objective"], ".12g")
            + " starts="
            + str(result["selection"]["start_count"])
            + " physical_validation=0 production_ready=0"
        )
        return 0
    except Exception as error:
        try:
            if output_path.is_file():
                output_path.unlink()
        except OSError:
            pass
        print("R2_UCSD_INTERNAL_FIT_SCREEN_FAIL " + str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
