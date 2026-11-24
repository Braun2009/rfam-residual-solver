#!/usr/bin/env python3
"""Run the blocked C3Mech 44I0 versus RFAM two-step fit screen.

The detailed mechanism is a candidate reference, not experimental truth.  The
training composite is verified before selection.  The holdout manifest and all
holdout artifacts remain unopened until the selected training point is frozen.
"""

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
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple


BLOCKED_MARKER = "R2_C3MECH_FIT_SCREEN_BLOCKED"
FAIL_MARKER = "R2_C3MECH_FIT_SCREEN_FAIL"

CONTRACT_AUTHORITY = {
    "path": "realFluidAllMach/tests/reference/r2_c3mech_fit_screen_contract_v1.json",
    "sha256": "aa4565afa3312915cf7e7b1532dd13d06ea4ba385a328125894eea7cc36fdbb4",
    "byte_count": 5762,
}
EXPECTED_INPUT_COMPOSITES = {
    "training_and_shared": {
        "manifest_path": "realFluidAllMach/tests/reference/r2_c3mech_fit_training_inputs_v1_composite.json",
        "manifest_sha256": "b2379dc8617afb070ba329e20475223ded25eec08f4486acae3266d86274ac1e",
        "external_composite_sha256": "375cc1071029ad2ea8792964cdf1e63de3dc89de4897797a3387b1fdc506a971",
        "artifact_count": 17,
        "verified_before_selection": True,
    },
    "holdout": {
        "manifest_path": "realFluidAllMach/tests/reference/r2_c3mech_fit_holdout_inputs_v1_composite.json",
        "manifest_sha256": "aef881cb31309a8b5af40dc2507e4100a085eec0246c0a44d9d4fdab8a12a6a6",
        "external_composite_sha256": "6e526134f88f57692e19017f9e73a37aaf48f4bf3fa1db493f4a11edb209ba20",
        "artifact_count": 4,
        "verified_only_after_training_selection_frozen": True,
    },
}
BOUND_PYTHON_AUTHORITIES = {
    "artifact_composite": {
        "path": "realFluidAllMach/tools/artifact_composite.py",
        "sha256": "15dff1addf9adc0adfcb167467663e68a11871fa944601e1d25df361ee14922d",
        "size": 21759,
    },
    "run_r2_0d_rfam_protocol": {
        "path": "realFluidAllMach/tools/run_r2_0d_rfam_protocol.py",
        "sha256": "5a4642c37ce4b07ca97db82a6af218d39e039b8a73e417d20efb6094d43f58a8",
        "size": 20748,
    },
    "run_r2_joint_fit_synthetic": {
        "path": "realFluidAllMach/tools/run_r2_joint_fit_synthetic.py",
        "sha256": "0669e6d7c2d81fc42a71dfd96688fdb13d1c821dbe953c8a6f1f9bbefd40f99d",
        "size": 71874,
    },
}

SPECIES = ("C3H8", "O2", "CO", "CO2", "H2O", "N2")
HISTORY_OBSERVABLES = ("qdot", "temperature", "C3H8", "O2", "CO", "CO2", "H2O")
PARAMETER_NAMES = (
    "R1_A_SI",
    "R1_Ea_J_per_mol",
    "R1_order_C3H8",
    "R1_order_O2",
    "R2_A_SI",
    "R2_Ea_J_per_mol",
)
DEFAULT_PARAMETERS = (1.0e10, 125520.0, 0.1, 1.65, 2.0e12, 167360.0)

EXPECTED_EXECUTION = {
    "platform": "WSL-Linux-x86_64",
    "python_major_minor": "3.8",
    "runner_timeout_seconds": 30,
    "parallel_case_processes": 3,
}
EXPECTED_MECHANISM = {
    "model_id": "C3Mech-V4.0.1-44I0",
    "role": "rights-cleared-detailed-reference-candidate-not-selected-for-production",
    "license": "CC-BY-4.0",
    "training_reference_path": "tmp/a408-r2-c3mech-44i0-training-v2.json",
    "holdout_reference_path": "tmp/a408-r2-c3mech-44i0-holdout-v2.json",
    "thermo_min_k": 400.0,
    "thermo_max_k": 2000.0,
    "eligible_history_policy": (
        "score-only-samples-with-status-eligible-within-mechanism-thermo-intersection"
    ),
    "outside_endpoint_and_later_samples_scored": False,
    "runner_horizon_policy": (
        "append-reference-thermo-envelope-first-crossing-time-as-unscored-terminal-"
        "only-when-reference-ignition-delay-is-active"
    ),
    "runner_terminal_sentinel_history_scored": False,
}
EXPECTED_RUNNER = {
    "path": "tmp/rfam_r2_0d_fit_runner",
    "parameter_schema": "rfam.r2-joint-fit-parameters/v1",
    "species_order": list(SPECIES),
    "energy_basis": "sensible-enthalpy-plus-explicit-formation-heat",
    "phase": "ideal-gas-constant-pressure-rfam-runner",
}
EXPECTED_PARTITIONS = {
    "training_case_ids": ["train_phi0p8_diluted", "train_phi1p0", "train_phi1p2"],
    "holdout_case_ids": ["holdout_phi0p6_diluted", "holdout_phi1p1", "holdout_phi1p4"],
    "holdout_opened_during_selection": False,
    "holdout_objective_calls_after_selection": 1,
}
EXPECTED_PARAMETERS = (
    ("R1_A_SI", "log10", 1.0e5, 1.0e20),
    ("R1_Ea_J_per_mol", "linear", 0.0, 3.0e5),
    ("R1_order_C3H8", "linear", 0.0, 2.0),
    ("R1_order_O2", "linear", 0.0, 2.5),
    ("R2_A_SI", "log10", 1.0e5, 1.0e20),
    ("R2_Ea_J_per_mol", "linear", 0.0, 3.0e5),
)
EXPECTED_OBJECTIVE = {
    "observable_order": ["ignition_delay", *HISTORY_OBSERVABLES],
    "component_weights": {
        "ignition_delay": 0.30,
        "qdot": 0.15,
        "temperature": 0.20,
        "C3H8": 0.10,
        "O2": 0.10,
        "CO": 0.05,
        "CO2": 0.05,
        "H2O": 0.05,
    },
    "ignition_delay_residual": "log-prediction-over-reference-divided-by-0.25",
    "history_component": (
        "mean-over-eligible-samples-of-square-of-(prediction-minus-reference)-"
        "divided-by-scale"
    ),
    "qdot_scale": "maximum-absolute-eligible-reference-or-1e8-W-m3-floor",
    "temperature_scale_k": 50.0,
    "species_mass_fraction_scale": 0.02,
    "missing_reference_component_policy": (
        "omit-and-renormalize-active-component-weights-per-case"
    ),
    "case_objective": (
        "math-fsum-of-active-component-weight-times-component-divided-by-"
        "math-fsum-of-active-component-weights"
    ),
    "case_aggregation": "arithmetic-mean-in-frozen-case-order",
    "reduction": "python-math-fsum-binary64-in-case-observable-sample-order",
    "invalid_candidate_policy": "invalid-point-never-replaces-incumbent",
}
EXPECTED_OPTIMIZER = {
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
    "global_selection": (
        "lowest-training-objective-then-lexicographically-smallest-transformed-vector"
    ),
}
EXPECTED_APPROVALS = {
    "diagnostic_optimizer_execution_allowed": True,
    "mechanism_selected_for_runtime": False,
    "experimental_training_data_selected": False,
    "experimental_holdout_data_selected": False,
    "parameter_vector_publishable": False,
    "physical_validation_approved": False,
    "hem_chemistry_enabled": False,
    "fluent_execution_approved": False,
    "reacting_hot_start_approved": False,
    "production_ready": False,
}
EXPECTED_STOP_CONDITIONS = [
    "training or holdout input composite or manifest raw hash mismatch",
    "duplicate key unknown field nonfinite number or schema drift",
    (
        "holdout composite verified or any holdout artifact opened before the "
        "training selection is frozen"
    ),
    "any reference outside-envelope endpoint or later sample is scored",
    "any nonfinite out-of-bound or invalid runner result replaces an incumbent",
    (
        "fit result represented as experimental validation parameter publication "
        "HEM Fluent or production authority"
    ),
]


class ScreenError(RuntimeError):
    """A fail-closed identity, schema, runner, or objective violation."""


# Loaded only after the training/shared composite passes.
ObjectiveEvaluation: Any = None
ParameterSpec: Any = None
SearchConfig: Any = None
parse_runner_output: Any = None
radical_inverse: Any = None
run_pattern_search: Any = None


def strict_pairs(pairs: Sequence[Tuple[str, Any]]) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ScreenError("duplicate-json-key:" + key)
        result[key] = value
    return result


def reject_constant(value: str) -> None:
    raise ScreenError("nonfinite-json-token:" + value)


def load_json(path: Path, label: str) -> Tuple[Dict[str, Any], bytes]:
    try:
        raw = path.read_bytes()
        document = json.loads(
            raw.decode("ascii"),
            object_pairs_hook=strict_pairs,
            parse_constant=reject_constant,
        )
    except ScreenError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError, TypeError, ValueError) as error:
        raise ScreenError(label + "-read-failed:" + str(error)) from error
    if not isinstance(document, dict):
        raise ScreenError(label + "-root-not-object")
    return document, raw


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def exact_keys(value: Dict[str, Any], expected: Sequence[str], label: str) -> None:
    if not isinstance(value, dict) or set(value) != set(expected):
        raise ScreenError("key-schema-mismatch:" + label)


def exact(actual: Any, expected: Any, label: str) -> None:
    if type(actual) is not type(expected):
        raise ScreenError("contract-mismatch:" + label)
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise ScreenError("contract-mismatch:" + label)
        for key in expected:
            exact(actual[key], expected[key], label + "." + key)
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ScreenError("contract-mismatch:" + label)
        for index, (left, right) in enumerate(zip(actual, expected)):
            exact(left, right, label + "[" + str(index) + "]")
    elif actual != expected:
        raise ScreenError("contract-mismatch:" + label)


def finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ScreenError("not-number:" + label)
    result = float(value)
    if not math.isfinite(result):
        raise ScreenError("nonfinite:" + label)
    return result


def resolve_under(root: Path, relative: str, label: str) -> Path:
    if (
        not isinstance(relative, str)
        or not relative
        or "\\" in relative
        or relative.startswith("/")
        or any(part in ("", ".", "..") for part in relative.split("/"))
    ):
        raise ScreenError("invalid-relative-path:" + label)
    root = root.resolve()
    candidate = (root / Path(*relative.split("/"))).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as error:
        raise ScreenError("path-escapes-root:" + label) from error
    return candidate


def validate_contract(
    document: Dict[str, Any],
    build_search: bool = True,
) -> Tuple[Tuple[Any, ...], Any, List[Tuple[str, Tuple[float, ...]]]]:
    exact_keys(
        document,
        [
            "schema", "contract_version", "status", "scope", "input_composites",
            "execution", "mechanism_reference", "runner", "partitions", "parameters",
            "objective", "optimizer", "approvals", "stop_conditions",
        ],
        "contract",
    )
    exact(document["schema"], "rfam.r2-c3mech-fit-screen-contract/v1", "schema")
    exact(document["contract_version"], 1, "contract_version")
    exact(
        document["status"],
        "blocked-detailed-mechanism-candidate-fit-screen-only",
        "status",
    )
    exact(
        document["scope"],
        "c3mech-44i0-valid-prefix-training-fit-and-post-selection-holdout-no-experiment-no-publication",
        "scope",
    )
    exact(document["input_composites"], EXPECTED_INPUT_COMPOSITES, "input_composites")

    exact(document["execution"], EXPECTED_EXECUTION, "execution")
    exact(document["mechanism_reference"], EXPECTED_MECHANISM, "mechanism_reference")
    exact(document["runner"], EXPECTED_RUNNER, "runner")
    exact(document["partitions"], EXPECTED_PARTITIONS, "partitions")
    exact(document["objective"], EXPECTED_OBJECTIVE, "objective")
    exact(document["optimizer"], EXPECTED_OPTIMIZER, "optimizer")
    exact(document["approvals"], EXPECTED_APPROVALS, "approvals")
    exact(document["stop_conditions"], EXPECTED_STOP_CONDITIONS, "stop_conditions")

    parameters = document["parameters"]
    if not isinstance(parameters, list) or len(parameters) != len(EXPECTED_PARAMETERS):
        raise ScreenError("parameter-count-mismatch")
    if build_search and (ParameterSpec is None or SearchConfig is None or radical_inverse is None):
        raise ScreenError("optimizer-helpers-not-loaded")
    specs = []
    for index, (entry, expected_parameter) in enumerate(zip(parameters, EXPECTED_PARAMETERS)):
        exact_keys(entry, ["name", "transform", "lower", "upper"], "parameter")
        name, transform, lower, upper = expected_parameter
        exact(entry, {"name": name, "transform": transform, "lower": lower, "upper": upper}, "parameter[" + str(index) + "]")
        transformed_lower = math.log10(lower) if transform == "log10" else lower
        transformed_upper = math.log10(upper) if transform == "log10" else upper
        if build_search:
            specs.append(ParameterSpec(name, transform, lower, upper, transformed_lower, transformed_upper))

    if not build_search:
        return tuple(), None, []

    optimizer = document["optimizer"]
    starts = [("midpoint", tuple(0.5 for _ in specs))]
    for index in optimizer["starts"]["halton_indices"]:
        starts.append((
            "halton-" + str(index),
            tuple(radical_inverse(index, base) for base in optimizer["starts"]["halton_bases"]),
        ))
    config = SearchConfig(
        optimizer["initial_normalized_step"],
        optimizer["step_contraction"],
        optimizer["minimum_normalized_step"],
        optimizer["maximum_objective_evaluations_per_start"],
        len(specs),
    )
    return tuple(specs), config, starts


def verify_runtime() -> None:
    actual = {
        "platform": "WSL-" + platform.system() + "-" + platform.machine(),
        "python_major_minor": str(sys.version_info.major) + "." + str(sys.version_info.minor),
        "runner_timeout_seconds": 30,
        "parallel_case_processes": 3,
    }
    if actual != EXPECTED_EXECUTION:
        raise ScreenError("execution-runtime-mismatch:" + json.dumps(actual, sort_keys=True))


def verify_composite(workspace_root: Path, authority: Dict[str, Any], label: str) -> Path:
    manifest_path = resolve_under(workspace_root, authority["manifest_path"], label + ".path")
    try:
        raw = manifest_path.read_bytes()
    except OSError as error:
        raise ScreenError("composite-manifest-open-failed:" + label) from error
    if sha256(raw) != authority["manifest_sha256"]:
        raise ScreenError("composite-manifest-sha256-mismatch:" + label)
    try:
        composite = importlib.import_module("artifact_composite")
        expected_module = resolve_under(
            workspace_root,
            BOUND_PYTHON_AUTHORITIES["artifact_composite"]["path"],
            "artifact-composite-module",
        )
        if Path(composite.__file__).resolve() != expected_module:
            raise ScreenError("artifact-composite-module-path-mismatch")
        actual = composite.verify_manifest(
            workspace_root, manifest_path, authority["external_composite_sha256"]
        )
    except Exception as error:
        raise ScreenError("composite-verification-failed:" + label + ":" + str(error)) from error
    if actual != authority["external_composite_sha256"]:
        raise ScreenError("composite-result-mismatch:" + label)
    document, _ = load_json(manifest_path, label + "-manifest")
    if document.get("artifact_count") != authority["artifact_count"]:
        raise ScreenError("composite-artifact-count-mismatch:" + label)
    return manifest_path


def preflight_bound_module(
    workspace_root: Path,
    manifest_path: Path,
    module_name: str,
) -> Path:
    authority = BOUND_PYTHON_AUTHORITIES[module_name]
    document, _ = load_json(manifest_path, "training-manifest-preflight")
    artifacts = document.get("artifacts")
    matches = [
        entry
        for entry in artifacts
        if isinstance(entry, dict) and entry.get("path") == authority["path"]
    ] if isinstance(artifacts, list) else []
    if len(matches) != 1 or matches[0] != {
        "id": "artifact:" + authority["path"],
        "path": authority["path"],
        "sha256": authority["sha256"],
        "size": authority["size"],
    }:
        raise ScreenError("bound-module-manifest-mismatch:" + module_name)
    path = resolve_under(workspace_root, authority["path"], module_name + ".path")
    try:
        raw = path.read_bytes()
    except OSError as error:
        raise ScreenError("bound-module-open-failed:" + module_name) from error
    if len(raw) != authority["size"] or sha256(raw) != authority["sha256"]:
        raise ScreenError("bound-module-content-mismatch:" + module_name)
    return path


def manifest_artifact(workspace_root: Path, manifest_path: Path, relative: str, label: str) -> Path:
    document, _ = load_json(manifest_path, label + "-manifest")
    artifacts = document.get("artifacts")
    matches = [entry for entry in artifacts if isinstance(entry, dict) and entry.get("path") == relative] if isinstance(artifacts, list) else []
    if len(matches) != 1:
        raise ScreenError("manifest-artifact-identity-mismatch:" + label)
    path = resolve_under(workspace_root, relative, label + ".path")
    try:
        raw = path.read_bytes()
    except OSError as error:
        raise ScreenError("manifest-artifact-open-failed:" + label) from error
    entry = matches[0]
    if len(raw) != entry.get("size") or sha256(raw) != entry.get("sha256"):
        raise ScreenError("manifest-artifact-content-mismatch:" + label)
    return path


def load_bound_helpers(workspace_root: Path, training_manifest: Path) -> None:
    global ObjectiveEvaluation, ParameterSpec, SearchConfig
    global parse_runner_output, radical_inverse, run_pattern_search
    if ObjectiveEvaluation is not None:
        return
    expected_protocol = preflight_bound_module(
        workspace_root, training_manifest, "run_r2_0d_rfam_protocol"
    )
    expected_optimizer = preflight_bound_module(
        workspace_root, training_manifest, "run_r2_joint_fit_synthetic"
    )
    protocol = importlib.import_module("run_r2_0d_rfam_protocol")
    optimizer = importlib.import_module("run_r2_joint_fit_synthetic")
    if Path(protocol.__file__).resolve() != expected_protocol:
        raise ScreenError("protocol-helper-module-path-mismatch")
    if Path(optimizer.__file__).resolve() != expected_optimizer:
        raise ScreenError("optimizer-helper-module-path-mismatch")
    ObjectiveEvaluation = optimizer.ObjectiveEvaluation
    ParameterSpec = optimizer.ParameterSpec
    SearchConfig = optimizer.SearchConfig
    parse_runner_output = protocol.parse_runner_output
    radical_inverse = optimizer.radical_inverse
    run_pattern_search = optimizer.run_pattern_search


@dataclass(frozen=True)
class PreparedCase:
    case_id: str
    partition: str
    temperature_k: float
    pressure_pa: float
    initial_mass_fractions: Tuple[float, ...]
    eligible_samples: Tuple[Dict[str, Any], ...]
    runner_times_s: Tuple[float, ...]
    reference_ignition_delay_s: Optional[float]


def _validate_eligible_sample(sample: Any, expected_time: float, label: str) -> Dict[str, Any]:
    if not isinstance(sample, dict):
        raise ScreenError("reference-sample-not-object:" + label)
    exact_keys(
        sample,
        [
            "density_kg_m3", "heat_release_rate_w_m3", "major_species_mass_fractions",
            "pressure_pa", "status", "temperature_k", "time_s",
        ],
        label,
    )
    if sample["status"] != "eligible-within-mechanism-thermo-intersection":
        raise ScreenError("reference-sample-status-mismatch:" + label)
    if finite(sample["time_s"], label + ".time") != expected_time:
        raise ScreenError("reference-sample-time-mismatch:" + label)
    temperature = finite(sample["temperature_k"], label + ".temperature")
    if not 400.0 <= temperature <= 2000.0:
        raise ScreenError("reference-sample-temperature-out-of-envelope:" + label)
    if finite(sample["pressure_pa"], label + ".pressure") <= 0.0 or finite(sample["density_kg_m3"], label + ".density") <= 0.0:
        raise ScreenError("reference-sample-state-invalid:" + label)
    finite(sample["heat_release_rate_w_m3"], label + ".qdot")
    composition = sample["major_species_mass_fractions"]
    exact_keys(composition, list(SPECIES), label + ".composition")
    values = [finite(composition[name], label + ".Y." + name) for name in SPECIES]
    # C3Mech publishes this six-species projection, not a closed composition;
    # intermediates and radicals account for the missing mass after ignition.
    projected_sum = math.fsum(values)
    if any(value < 0.0 for value in values) or projected_sum <= 0.0 or projected_sum > 1.0 + 2.0e-8:
        raise ScreenError("reference-sample-composition-invalid:" + label)
    return sample


def validate_reference(
    document: Dict[str, Any],
    expected_partition: str,
    expected_case_ids: Sequence[str],
) -> Tuple[List[PreparedCase], Tuple[float, ...]]:
    exact_keys(
        document,
        [
            "approvals", "cases", "major_species", "model", "partition",
            "prefix_diagnostic", "reactor", "runtime_identity", "sample_times_s",
            "schema", "source", "status",
        ],
        "reference",
    )
    if (
        document["schema"] != "rfam.r2-c3mech-cantera-reference/v2"
        or document["partition"] != expected_partition
        or document["major_species"] != list(SPECIES)
        or document["status"] != "OUT_OF_MECHANISM_THERMO_ENVELOPE"
        or document["model"].get("model_id") != "44I0"
    ):
        raise ScreenError("reference-origin-mismatch:" + expected_partition)
    exact_keys(
        document["model"],
        [
            "author_header_reaction_count", "cantera_loaded_reaction_count",
            "count_discrepancy_resolved_as_equivalence", "counting_scope", "model_id",
            "role", "species_count", "thermo_max_k", "thermo_min_k",
            "yaml_reaction_entry_count",
        ],
        "reference.model",
    )
    if (
        document["model"]["thermo_min_k"] != 400.0
        or document["model"]["thermo_max_k"] != 2000.0
        or document["model"]["role"] != "primary-lt-ht-reference-candidate-not-selected"
    ):
        raise ScreenError("reference-model-policy-mismatch:" + expected_partition)
    exact_keys(
        document["source"],
        [
            "case_manifest_path", "case_manifest_sha256", "generator_sha256",
            "mechanism_path", "mechanism_sha256", "mechanism_size",
        ],
        "reference.source",
    )
    exact_keys(
        document["reactor"],
        [
            "atol", "energy", "ignition_temperature_advance_limit_k", "model",
            "outside_envelope_policy", "phase", "rtol", "temperature_advance_limit_k",
        ],
        "reference.reactor",
    )
    if document["reactor"] != {
        "atol": 1.0e-18,
        "energy": "on",
        "ignition_temperature_advance_limit_k": 1.0,
        "model": "IdealGasConstPressureReactor",
        "outside_envelope_policy": (
            "stop-at-first-detected-exit; publish-no-state; leave-later-samples-unevaluated"
        ),
        "phase": "ideal-gas",
        "rtol": 1.0e-10,
        "temperature_advance_limit_k": 10.0,
    }:
        raise ScreenError("reference-reactor-policy-mismatch:" + expected_partition)
    exact_keys(
        document["runtime_identity"],
        [
            "cantera_extension", "cantera_extension_sha256", "cantera_extension_size",
            "cantera_module", "cantera_module_sha256", "cantera_module_size",
            "cantera_version", "cantera_wheel", "cantera_wheel_sha256",
            "cantera_wheel_size", "execution_authority", "python_executable",
            "python_sha256", "python_size", "wheel_to_installed_tree_equivalence_proven",
        ],
        "reference.runtime_identity",
    )
    if (
        document["runtime_identity"]["cantera_version"] != "3.2.0"
        or document["runtime_identity"]["execution_authority"]
        != "bound-installed-module-and-extension"
        or document["runtime_identity"]["wheel_to_installed_tree_equivalence_proven"] is not False
    ):
        raise ScreenError("reference-runtime-policy-mismatch:" + expected_partition)
    approvals = document["approvals"]
    exact_keys(
        approvals,
        [
            "fit_executed", "holdout_used_for_selection", "mechanism_selected_for_fit",
            "physical_validation_approved", "production_ready",
        ],
        "reference.approvals",
    )
    if any(value is not False for value in approvals.values()):
        raise ScreenError("reference-approval-enabled:" + expected_partition)
    if document["prefix_diagnostic"] != {
        "eligible_samples_only": True,
        "later_samples_evaluated": False,
        "outside_endpoint_state_published": False,
    }:
        raise ScreenError("reference-prefix-policy-mismatch:" + expected_partition)
    times = tuple(finite(value, "reference.sample-time") for value in document["sample_times_s"])
    if not times or times[0] != 0.0 or any(right <= left for left, right in zip(times, times[1:])):
        raise ScreenError("reference-sample-times-invalid:" + expected_partition)
    cases = document["cases"]
    if not isinstance(cases, list) or [case.get("id") for case in cases if isinstance(case, dict)] != list(expected_case_ids):
        raise ScreenError("reference-case-order-mismatch:" + expected_partition)

    prepared = []
    for case in cases:
        case_id = case["id"]
        exact_keys(
            case,
            [
                "equivalence_ratio", "fuel", "full_protocol_status", "id",
                "ignition_delay", "ignition_delay_definition",
                "initial_major_species_mass_fractions", "oxidizer", "partition",
                "pressure_pa", "samples", "temperature_k", "thermo_envelope_exit",
            ],
            "case." + case_id,
        )
        if case["partition"] != expected_partition or case["full_protocol_status"] != "OUT_OF_MECHANISM_THERMO_ENVELOPE":
            raise ScreenError("reference-case-partition-mismatch:" + case_id)
        initial = case["initial_major_species_mass_fractions"]
        exact_keys(initial, list(SPECIES), case_id + ".initial")
        initial_values = tuple(finite(initial[name], case_id + ".initial." + name) for name in SPECIES)
        if any(value < 0.0 for value in initial_values) or abs(math.fsum(initial_values) - 1.0) > 2.0e-9:
            raise ScreenError("reference-initial-composition-invalid:" + case_id)
        raw_samples = case["samples"]
        if not isinstance(raw_samples, list) or len(raw_samples) != len(times):
            raise ScreenError("reference-sample-count-mismatch:" + case_id)
        eligible: List[Dict[str, Any]] = []
        noneligible_seen = False
        outside_seen = False
        for index, (sample, expected_time) in enumerate(zip(raw_samples, times)):
            if not isinstance(sample, dict) or sample.get("time_s") != expected_time:
                raise ScreenError("reference-sample-time-mismatch:" + case_id + ":" + str(index))
            status = sample.get("status")
            if status == "eligible-within-mechanism-thermo-intersection":
                if noneligible_seen:
                    raise ScreenError("reference-eligible-prefix-broken:" + case_id)
                eligible.append(_validate_eligible_sample(sample, expected_time, case_id + ".sample[" + str(index) + "]"))
            elif status == "outside-thermo-envelope-unscored":
                if outside_seen or not eligible:
                    raise ScreenError("reference-outside-endpoint-invalid:" + case_id)
                exact_keys(
                    sample,
                    ["exit_detection", "status", "time_s"],
                    case_id + ".outside-endpoint",
                )
                if sample["exit_detection"] != case["thermo_envelope_exit"]:
                    raise ScreenError("reference-exit-detection-mismatch:" + case_id)
                noneligible_seen = True
                outside_seen = True
            elif status == "not-evaluated-after-thermo-envelope-exit":
                if not outside_seen:
                    raise ScreenError("reference-later-sample-before-exit:" + case_id)
                exact_keys(
                    sample,
                    ["status", "time_s"],
                    case_id + ".later-sample",
                )
                noneligible_seen = True
            else:
                raise ScreenError("reference-sample-status-invalid:" + case_id)
        if not eligible or not outside_seen:
            raise ScreenError("reference-prefix-or-exit-missing:" + case_id)

        envelope_exit = case["thermo_envelope_exit"]
        exact_keys(
            envelope_exit,
            [
                "first_crossing_estimate", "requested_sample_time_s", "status",
                "thermo_max_k", "thermo_min_k",
            ],
            case_id + ".thermo-envelope-exit",
        )
        crossing = envelope_exit["first_crossing_estimate"]
        exact_keys(
            crossing,
            ["method", "pressure_pa", "publishable_chemistry_state", "temperature_k", "time_s"],
            case_id + ".first-crossing",
        )
        if (
            envelope_exit["status"] != "thermo-envelope-exit-detected-no-endpoint-state-published"
            or envelope_exit["thermo_min_k"] != 400.0
            or envelope_exit["thermo_max_k"] != 2000.0
            or crossing["method"] != "10-K-advance-limited-linear-temperature-interpolation"
            or crossing["publishable_chemistry_state"] is not False
            or crossing["temperature_k"] != 2000.0
            or crossing["pressure_pa"] != case["pressure_pa"]
        ):
            raise ScreenError("reference-thermo-exit-policy-mismatch:" + case_id)

        ignition = case["ignition_delay"]
        if not isinstance(ignition, dict) or "status" not in ignition:
            raise ScreenError("reference-ignition-invalid:" + case_id)
        reference_idt: Optional[float]
        runner_times = [float(sample["time_s"]) for sample in eligible]
        if ignition["status"] == "eligible-within-mechanism-thermo-intersection":
            exact_keys(
                ignition,
                ["ignition_delay_s", "method", "status", "threshold_temperature_k"],
                case_id + ".ignition",
            )
            if ignition["method"] != "1-K-advance-limited-linear-threshold-interpolation":
                raise ScreenError("reference-ignition-method-mismatch:" + case_id)
            reference_idt = finite(ignition.get("ignition_delay_s"), case_id + ".ignition-delay")
            if reference_idt <= 0.0:
                raise ScreenError("reference-ignition-delay-invalid:" + case_id)
            try:
                sentinel = finite(
                    case["thermo_envelope_exit"]["first_crossing_estimate"]["time_s"],
                    case_id + ".terminal-sentinel",
                )
            except (KeyError, TypeError) as error:
                raise ScreenError("reference-terminal-sentinel-missing:" + case_id) from error
            if (
                sentinel <= runner_times[-1]
                or sentinel > times[len(eligible)]
                or reference_idt >= sentinel
            ):
                raise ScreenError("reference-terminal-sentinel-invalid:" + case_id)
            runner_times.append(sentinel)
        elif ignition["status"] == "threshold-at-or-outside-thermo-upper-bound-unscored":
            exact_keys(
                ignition,
                ["status", "thermo_max_k", "threshold_temperature_k"],
                case_id + ".ignition",
            )
            if "ignition_delay_s" in ignition:
                raise ScreenError("reference-unscored-ignition-has-value:" + case_id)
            reference_idt = None
        else:
            raise ScreenError("reference-ignition-status-invalid:" + case_id)

        prepared.append(PreparedCase(
            case_id,
            expected_partition,
            finite(case["temperature_k"], case_id + ".temperature"),
            finite(case["pressure_pa"], case_id + ".pressure"),
            initial_values,
            tuple(eligible),
            tuple(runner_times),
            reference_idt,
        ))
    return prepared, times


def strip_fit_parameter_header(stdout: str, parameters: Sequence[float], case_id: str) -> str:
    lines = stdout.splitlines(keepends=True)
    expected = ["fit_parameter_schema=rfam.r2-joint-fit-parameters/v1\n"]
    expected.extend(
        "fit_" + name + "=" + format(float(value), ".17g") + "\n"
        for name, value in zip(PARAMETER_NAMES, parameters)
    )
    if len(lines) < len(expected) or lines[: len(expected)] != expected:
        raise ScreenError("runner-fit-header-mismatch:" + case_id)
    return "".join(lines[len(expected) :])


def score_case(
    case: PreparedCase,
    headers: Dict[str, str],
    samples: Sequence[Dict[str, Any]],
    weights: Dict[str, float],
) -> Tuple[float, Tuple[float, ...], Dict[str, float]]:
    eligible_count = len(case.eligible_samples)
    if len(samples) != len(case.runner_times_s) or eligible_count > len(samples):
        raise ScreenError("prediction-sample-count-mismatch:" + case.case_id)
    components: Dict[str, float] = {}
    if case.reference_ignition_delay_s is not None:
        if headers.get("ignition_detected") != "1":
            raise ScreenError("prediction-ignition-not-detected:" + case.case_id)
        predicted_idt = float(headers["ignition_delay_s"])
        if not math.isfinite(predicted_idt) or predicted_idt <= 0.0:
            raise ScreenError("prediction-ignition-delay-invalid:" + case.case_id)
        residual = math.log(predicted_idt / case.reference_ignition_delay_s) / 0.25
        components["ignition_delay"] = residual * residual

    scored_predictions = samples[:eligible_count]
    references = case.eligible_samples
    qdot_scale = max(
        max(abs(float(sample["heat_release_rate_w_m3"])) for sample in references),
        1.0e8,
    )
    species_index = {name: index for index, name in enumerate(SPECIES)}
    for observable in HISTORY_OBSERVABLES:
        if observable == "qdot":
            predicted = [float(sample["qdot_w_m3"]) for sample in scored_predictions]
            target = [float(sample["heat_release_rate_w_m3"]) for sample in references]
            scale = qdot_scale
        elif observable == "temperature":
            predicted = [float(sample["temperature_k"]) for sample in scored_predictions]
            target = [float(sample["temperature_k"]) for sample in references]
            scale = 50.0
        else:
            predicted = [float(sample["mass_fractions"][species_index[observable]]) for sample in scored_predictions]
            target = [float(sample["major_species_mass_fractions"][observable]) for sample in references]
            scale = 0.02
        squared = []
        for left, right in zip(predicted, target):
            residual = (left - right) / scale
            if not math.isfinite(residual):
                raise ScreenError("prediction-component-nonfinite:" + case.case_id + ":" + observable)
            squared.append(residual * residual)
        components[observable] = math.fsum(squared) / float(eligible_count)

    active_weight = math.fsum(weights[name] for name in components)
    if not math.isfinite(active_weight) or active_weight <= 0.0:
        raise ScreenError("active-component-weight-invalid:" + case.case_id)
    weighted_rows = tuple(
        math.sqrt(weights[name] / active_weight * components[name])
        for name in EXPECTED_OBJECTIVE["observable_order"]
        if name in components
    )
    objective = math.fsum(value * value for value in weighted_rows)
    return objective, weighted_rows, components


class R2C3MechObjective:
    def __init__(
        self,
        runner: Path,
        specs: Sequence[Any],
        training_cases: Sequence[PreparedCase],
        contract: Dict[str, Any],
    ) -> None:
        self.runner = runner
        self.specs = tuple(specs)
        self.case_by_id = {case.case_id: case for case in training_cases}
        self.case_ids = {"training": tuple(contract["partitions"]["training_case_ids"])}
        self.expected_holdout_ids = tuple(contract["partitions"]["holdout_case_ids"])
        self.timeout = int(contract["execution"]["runner_timeout_seconds"])
        self.workers = int(contract["execution"]["parallel_case_processes"])
        self.weights = dict(contract["objective"]["component_weights"])
        self.selection_active = True
        self.objective_calls: Counter[str] = Counter()
        self.case_process_attempts: Counter[str] = Counter()
        self.case_process_successes: Counter[str] = Counter()
        self.invalid_reasons: Counter[str] = Counter()
        self.holdout_access_attempts_during_selection = 0
        self.case_component_cache: Dict[
            Tuple[str, Tuple[str, ...]], Dict[str, Dict[str, float]]
        ] = {}

    def parameters(self, point: Sequence[float]) -> Tuple[float, ...]:
        if len(point) != len(self.specs):
            raise ScreenError("parameter-point-size-mismatch")
        return tuple(spec.physical_from_u(float(value)) for spec, value in zip(self.specs, point))

    def freeze_selection(self) -> None:
        if not self.selection_active:
            raise ScreenError("selection-already-frozen")
        self.selection_active = False

    def case_components(
        self,
        split: str,
        point: Sequence[float],
    ) -> Dict[str, Dict[str, float]]:
        key = (split, tuple(float(value).hex() for value in point))
        if key not in self.case_component_cache:
            raise ScreenError("case-component-cache-miss:" + split)
        return {
            case_id: dict(components)
            for case_id, components in self.case_component_cache[key].items()
        }

    def install_holdout(self, cases: Sequence[PreparedCase]) -> None:
        if self.selection_active:
            raise ScreenError("holdout-installed-during-selection")
        if "holdout" in self.case_ids:
            raise ScreenError("holdout-installed-more-than-once")
        ids = tuple(case.case_id for case in cases)
        if ids != self.expected_holdout_ids:
            raise ScreenError("holdout-case-order-mismatch")
        self.case_by_id.update((case.case_id, case) for case in cases)
        self.case_ids["holdout"] = ids

    def _run_case(self, parameters: Sequence[float], case: PreparedCase) -> Tuple[Dict[str, str], List[Dict[str, Any]]]:
        command = [
            str(self.runner),
            *[format(float(value), ".17g") for value in parameters],
            case.case_id,
            format(case.temperature_k, ".17g"),
            format(case.pressure_pa, ".17g"),
            str(len(case.runner_times_s)),
            "0",
            *[format(value, ".17g") for value in case.initial_mass_fractions],
            *[format(value, ".17g") for value in case.runner_times_s],
        ]
        result = subprocess.run(
            command,
            check=False,
            text=True,
            capture_output=True,
            timeout=self.timeout,
        )
        if result.returncode != 0:
            raise ScreenError("runner-rejected:" + case.case_id + ":rc=" + str(result.returncode))
        trajectory = strip_fit_parameter_header(result.stdout, parameters, case.case_id)
        headers, samples = parse_runner_output(
            trajectory, result.stderr, case.case_id, len(case.runner_times_s)
        )
        for index, (sample, expected_time) in enumerate(zip(samples, case.runner_times_s)):
            if abs(float(sample["time_s"]) - expected_time) > 2.0e-18:
                raise ScreenError("runner-time-mismatch:" + case.case_id + ":" + str(index))
        return headers, samples

    def evaluate(self, point: Sequence[float], split: str) -> Any:
        if split == "holdout" and self.selection_active:
            self.holdout_access_attempts_during_selection += 1
            raise ScreenError("holdout-accessed-during-selection")
        if split not in self.case_ids:
            raise ScreenError("unknown-objective-split:" + split)
        self.objective_calls[split] += 1
        ids = self.case_ids[split]
        parameters = self.parameters(point)
        try:
            self.case_process_attempts[split] += len(ids)
            with ThreadPoolExecutor(max_workers=self.workers) as executor:
                futures = [
                    executor.submit(self._run_case, parameters, self.case_by_id[case_id])
                    for case_id in ids
                ]
                predictions = [future.result() for future in futures]
            self.case_process_successes[split] += len(ids)
            residuals: List[float] = []
            case_objectives: Dict[str, float] = {}
            case_components: Dict[str, Dict[str, float]] = {}
            case_count = len(ids)
            for case_id, (headers, samples) in zip(ids, predictions):
                case_objective, case_rows, components = score_case(
                    self.case_by_id[case_id], headers, samples, self.weights
                )
                case_objectives[case_id] = case_objective
                case_components[case_id] = components
                residuals.extend(value / math.sqrt(float(case_count)) for value in case_rows)
            objective = math.fsum(value * value for value in residuals)
            if not math.isfinite(objective):
                raise ScreenError("objective-nonfinite")
            key = (split, tuple(float(value).hex() for value in point))
            self.case_component_cache[key] = case_components
            return ObjectiveEvaluation(True, objective, tuple(residuals), case_objectives)
        except (
            ScreenError,
            KeyError,
            IndexError,
            TypeError,
            ValueError,
            OverflowError,
            subprocess.TimeoutExpired,
        ) as error:
            reason = str(error)
            self.invalid_reasons[reason] += 1
            return ObjectiveEvaluation.invalid(reason)


def _transformed_key(point: Sequence[float], specs: Sequence[Any]) -> Tuple[float, ...]:
    return tuple(spec.transformed_from_u(value) for spec, value in zip(specs, point))


def rounded(value: Any) -> Any:
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ScreenError("result-nonfinite")
        return float(format(value, ".15g"))
    if isinstance(value, dict):
        return {key: rounded(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [rounded(item) for item in value]
    return value


def canonical_bytes(document: Dict[str, Any]) -> bytes:
    return (json.dumps(rounded(document), indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def admit_output_path(path: Path) -> Path:
    output = Path(os.path.abspath(os.fspath(path)))
    if not output.parent.is_dir():
        raise ScreenError("output-parent-not-directory")
    if os.path.lexists(str(output)):
        raise ScreenError("output-path-already-exists")
    return output


def write_atomic(path: Path, raw: bytes) -> None:
    if os.path.lexists(str(path)):
        raise ScreenError("output-path-already-exists")
    descriptor, temporary_name = tempfile.mkstemp(prefix=path.name + ".tmp.", dir=str(path.parent))
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise
    try:
        os.unlink(temporary_name)
    except OSError:
        pass


def execute(contract_path: Path) -> Dict[str, Any]:
    verify_runtime()
    repository_root = Path(__file__).resolve().parent.parent
    workspace_root = repository_root.parent
    canonical_contract = repository_root / "tests/reference/r2_c3mech_fit_screen_contract_v1.json"
    if contract_path.resolve() != canonical_contract.resolve():
        raise ScreenError("contract-path-mismatch")
    contract, contract_raw = load_json(contract_path, "contract")
    if (
        len(contract_raw) != CONTRACT_AUTHORITY["byte_count"]
        or sha256(contract_raw) != CONTRACT_AUTHORITY["sha256"]
    ):
        raise ScreenError("contract-content-mismatch")
    # Reject schema and authority redirection before opening any path named by
    # the contract or loading the bound optimizer/protocol helpers.
    validate_contract(contract, build_search=False)

    tools_text = str(repository_root / "tools")
    if tools_text not in sys.path:
        sys.path.insert(0, tools_text)
    training_manifest_path = resolve_under(
        workspace_root,
        contract["input_composites"]["training_and_shared"]["manifest_path"],
        "training-manifest-preflight",
    )
    preflight_bound_module(
        workspace_root, training_manifest_path, "artifact_composite"
    )
    training_manifest = verify_composite(
        workspace_root, contract["input_composites"]["training_and_shared"], "training"
    )
    load_bound_helpers(workspace_root, training_manifest)
    specs, config, starts = validate_contract(contract)

    training_path = manifest_artifact(
        workspace_root,
        training_manifest,
        contract["mechanism_reference"]["training_reference_path"],
        "training-reference",
    )
    runner_path = manifest_artifact(
        workspace_root, training_manifest, contract["runner"]["path"], "runner"
    )
    training_document, _ = load_json(training_path, "training-reference")
    training_cases, _ = validate_reference(
        training_document,
        "protocol-training",
        contract["partitions"]["training_case_ids"],
    )
    objective = R2C3MechObjective(runner_path, specs, training_cases, contract)

    default_u = tuple(
        spec.u_from_physical(value) for spec, value in zip(specs, DEFAULT_PARAMETERS)
    )
    baseline = objective.evaluate(default_u, "training")
    if not baseline.valid:
        raise ScreenError("baseline-training-invalid:" + baseline.reason)

    results = []
    for start_index, (kind, start) in enumerate(starts):
        results.append(run_pattern_search(
            start_index,
            kind,
            start,
            specs,
            config,
            lambda candidate: objective.evaluate(candidate, "training"),
        ))
    valid_results = [result for result in results if result.valid]
    if not valid_results:
        raise ScreenError("all-search-starts-invalid")
    selected = min(
        valid_results,
        key=lambda result: (result.evaluation.objective, _transformed_key(result.selected_u, specs)),
    )
    if objective.holdout_access_attempts_during_selection != 0:
        raise ScreenError("holdout-access-accounting-mismatch")

    objective.freeze_selection()
    holdout_manifest = verify_composite(
        workspace_root, contract["input_composites"]["holdout"], "holdout"
    )
    holdout_path = manifest_artifact(
        workspace_root,
        holdout_manifest,
        contract["mechanism_reference"]["holdout_reference_path"],
        "holdout-reference",
    )
    holdout_document, _ = load_json(holdout_path, "holdout-reference")
    holdout_cases, _ = validate_reference(
        holdout_document,
        "protocol-holdout",
        contract["partitions"]["holdout_case_ids"],
    )
    for key in ("major_species", "model", "reactor", "runtime_identity", "schema"):
        if training_document[key] != holdout_document[key]:
            raise ScreenError("training-holdout-metadata-mismatch:" + key)
    for key in ("generator_sha256", "mechanism_path", "mechanism_sha256", "mechanism_size"):
        if training_document["source"].get(key) != holdout_document["source"].get(key):
            raise ScreenError("training-holdout-source-mismatch:" + key)
    objective.install_holdout(holdout_cases)
    holdout = objective.evaluate(selected.selected_u, "holdout")
    if not holdout.valid:
        raise ScreenError("selected-holdout-invalid:" + holdout.reason)
    if objective.objective_calls["holdout"] != 1:
        raise ScreenError("holdout-objective-call-count-mismatch")
    verify_composite(
        workspace_root, contract["input_composites"]["training_and_shared"], "training-final"
    )
    verify_composite(
        workspace_root, contract["input_composites"]["holdout"], "holdout-final"
    )

    expected_training_calls = 1 + sum(result.evaluation_count for result in results)
    if (
        objective.objective_calls["training"] != expected_training_calls
        or objective.case_process_attempts["training"] != 3 * expected_training_calls
        or objective.case_process_attempts["holdout"] != 3
    ):
        raise ScreenError("runner-accounting-mismatch")

    return {
        "schema": "rfam.r2-c3mech-fit-screen-result/v1",
        "status": "blocked-detailed-mechanism-candidate-fit-screen-only",
        "source": {
            "contract_path": "realFluidAllMach/tests/reference/r2_c3mech_fit_screen_contract_v1.json",
            "contract_sha256": sha256(contract_raw),
            "training_composite": contract["input_composites"]["training_and_shared"],
            "holdout_composite": contract["input_composites"]["holdout"],
        },
        "baseline": {
            "parameter_vector": DEFAULT_PARAMETERS,
            "training_objective": baseline.objective,
            "training_case_objectives": baseline.case_objectives,
        },
        "selection": {
            "split": "training-only",
            "start_count": len(results),
            "valid_start_count": len(valid_results),
            "selected_start_index": selected.start_index,
            "selected_start_kind": selected.start_kind,
            "selected_normalized_coordinates": selected.selected_u,
            "selected_parameter_vector": objective.parameters(selected.selected_u),
            "training_objective": selected.evaluation.objective,
            "training_case_objectives": selected.evaluation.case_objectives,
            "training_case_components": objective.case_components(
                "training", selected.selected_u
            ),
            "holdout_objective": holdout.objective,
            "holdout_case_objectives": holdout.case_objectives,
            "holdout_case_components": objective.case_components(
                "holdout", selected.selected_u
            ),
        },
        "starts": [
            {
                "index": result.start_index,
                "kind": result.start_kind,
                "valid": result.valid,
                "objective": result.evaluation.objective if result.valid else None,
                "selected_normalized_coordinates": result.selected_u,
                "selected_parameter_vector": objective.parameters(result.selected_u),
                "evaluation_count": result.evaluation_count,
                "cache_hits": result.cache_hits,
                "accepted_coordinates": result.accepted_coordinates,
                "termination": result.termination,
                "trace_sha256": result.trace_sha256,
            }
            for result in results
        ],
        "runner_accounting": {
            "training_objective_calls": objective.objective_calls["training"],
            "holdout_access_attempts_during_selection": objective.holdout_access_attempts_during_selection,
            "holdout_objective_calls_after_selection": objective.objective_calls["holdout"],
            "training_case_process_attempts": objective.case_process_attempts["training"],
            "training_case_process_successes": objective.case_process_successes["training"],
            "holdout_case_process_attempts": objective.case_process_attempts["holdout"],
            "holdout_case_process_successes": objective.case_process_successes["holdout"],
            "invalid_candidate_reasons": dict(sorted(objective.invalid_reasons.items())),
        },
        "case_protocol": {
            case.case_id: {
                "eligible_history_samples": len(case.eligible_samples),
                "runner_samples": len(case.runner_times_s),
                "terminal_sentinel_appended": case.reference_ignition_delay_s is not None,
                "reference_ignition_delay_active": case.reference_ignition_delay_s is not None,
            }
            for case in training_cases + holdout_cases
        },
        "approvals": EXPECTED_APPROVALS,
        "claims": {
            "detailed_mechanism_candidate_screen_executed": True,
            "experimental_fit_executed": False,
            "parameter_vector_publishable": False,
            "physical_validation_approved": False,
            "hem_chemistry_enabled": False,
            "fluent_execution_approved": False,
            "production_ready": False,
        },
        "stop_conditions": EXPECTED_STOP_CONDITIONS,
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args(argv)
    try:
        output = admit_output_path(arguments.output)
        result = execute(arguments.contract.resolve())
        write_atomic(output, canonical_bytes(result))
        print(
            BLOCKED_MARKER
            + " training_objective="
            + format(result["selection"]["training_objective"], ".12g")
            + " holdout_objective="
            + format(result["selection"]["holdout_objective"], ".12g")
            + " starts="
            + str(result["selection"]["start_count"])
            + " publishable=0 production_ready=0"
        )
        return 0
    except Exception as error:
        reason = str(error).replace("\r", "\\r").replace("\n", "\\n")
        print(FAIL_MARKER + " " + reason, file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
