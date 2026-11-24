#!/usr/bin/env python3
"""Compare bounded C3Mech partitions with the fixed RFAM two-step runner."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

from run_r2_0d_rfam_protocol import ProtocolError, parse_runner_output


class ScreenError(RuntimeError):
    pass


SPECIES = ["C3H8", "O2", "CO", "CO2", "H2O", "N2"]
PLAN_SHA256 = "ccef1959d057a5618114511cda9ef3762979be7e897dc9d092af6c1c90ccb95b"
PROTOCOL_MODULE_SHA256 = "5a4642c37ce4b07ca97db82a6af218d39e039b8a73e417d20efb6094d43f58a8"
REFERENCE_BINDINGS = {
    ("44I0", "protocol-training"): {
        "path": "F:\\tmp\\a408-r2-c3mech-44i0-training-v2.json",
        "sha256": "a1298082963b4adafc5d02de38450e663f2483704847bc827c1f4585939ee853",
        "size": 16668,
    },
    ("3R00", "protocol-training"): {
        "path": "F:\\tmp\\a408-r2-c3mech-3r00-training-v2.json",
        "sha256": "42e7aeedb90b9ceadc782fa5680eb633314da2c227cc424dce395b4e19238a52",
        "size": 16651,
    },
    ("44I0", "protocol-holdout"): {
        "path": "F:\\tmp\\a408-r2-c3mech-44i0-holdout-v2.json",
        "sha256": "4961e2e76f0a7fe4b54ba1ea67d9dace7896416ed2b8633b56474a7046f342ac",
        "size": 16323,
    },
    ("3R00", "protocol-holdout"): {
        "path": "F:\\tmp\\a408-r2-c3mech-3r00-holdout-v2.json",
        "sha256": "bfda66faacc8335eff0c76e313ade9994d19b37ae800bc4214548c7da232f98d",
        "size": 16300,
    },
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ScreenError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        result = json.loads(
            path.read_text(encoding="utf-8"), object_pairs_hook=strict_object
        )
    except ScreenError:
        raise
    except Exception as exc:
        raise ScreenError(f"invalid-json:{label}:{type(exc).__name__}") from exc
    if not isinstance(result, dict):
        raise ScreenError(f"not-an-object:{label}")
    return result


def finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ScreenError(f"not-a-number:{label}")
    result = float(value)
    if not math.isfinite(result):
        raise ScreenError(f"nonfinite:{label}")
    return result


def stable(value: float, label: str) -> float:
    if not math.isfinite(value):
        raise ScreenError(f"nonfinite:{label}")
    return float(format(float(value), ".12g"))


def exact_path(path: Path, expected: str, label: str) -> Path:
    resolved = path.resolve()
    if str(resolved).casefold() != str(Path(expected).resolve()).casefold():
        raise ScreenError(f"path-mismatch:{label}")
    return resolved


def validate_bound_file(path: Path, binding: dict[str, Any], label: str) -> Path:
    resolved = exact_path(path, binding["path"], label)
    if not resolved.is_file():
        raise ScreenError(f"missing:{label}")
    if resolved.stat().st_size != binding["size"]:
        raise ScreenError(f"size-mismatch:{label}")
    if sha256(resolved) != binding["sha256"]:
        raise ScreenError(f"sha256-mismatch:{label}")
    return resolved


def model_from_plan(plan: dict[str, Any], model_id: str) -> dict[str, Any]:
    matches = [item for item in plan["mechanisms"] if item.get("model_id") == model_id]
    if len(matches) != 1:
        raise ScreenError(f"plan-model-missing:{model_id}")
    return matches[0]


def validate_plan(path: Path) -> dict[str, Any]:
    if sha256(path) != PLAN_SHA256:
        raise ScreenError("plan-sha256-mismatch")
    plan = load_json(path, "plan")
    expected_root = {
        "schema", "contract_version", "status", "scope",
        "repository_provenance", "runtime_identity", "protocol", "role_lock",
        "mechanisms", "reaction_count_semantics", "thermo_envelope_policy",
        "rfam_runner", "approvals",
    }
    if set(plan) != expected_root:
        raise ScreenError("plan-root-schema-mismatch")
    if plan.get("schema") != "rfam.r2-c3mech-candidate-screen-plan/v2":
        raise ScreenError("plan-schema-mismatch")
    if plan.get("contract_version") != 2:
        raise ScreenError("plan-version-mismatch")
    if plan.get("status") != "blocked-offline-prefix-diagnostic-only":
        raise ScreenError("plan-status-mismatch")
    role_lock = plan.get("role_lock")
    expected_role_lock = {
        "locked": True,
        "basis": "pre-trajectory rights and intended-scope screen only",
        "trajectory_based_selection_executed": False,
        "fit_executed": False,
        "training_results_may_change_roles": False,
        "holdout_results_may_change_roles": False,
        "holdout_access_stage": "after immutable role lock and training comparison completion",
        "primary_model_id": "44I0",
        "comparator_model_id": "3R00",
        "selected_model_id": None,
    }
    if role_lock != expected_role_lock:
        raise ScreenError("role-lock-mismatch")
    approvals = plan.get("approvals")
    if not isinstance(approvals, dict) or not approvals or any(
        type(value) is not bool or value for value in approvals.values()
    ):
        raise ScreenError("plan-approval-enabled")
    if plan.get("reaction_count_semantics", {}).get("header_and_yaml_counts_interchangeable") is not False:
        raise ScreenError("reaction-count-authority-open")
    runtime = plan.get("runtime_identity")
    if not isinstance(runtime, dict):
        raise ScreenError("runtime-identity-schema-mismatch")
    if runtime.get("execution_authority") != "bound-installed-module-and-extension":
        raise ScreenError("runtime-execution-authority-mismatch")
    if runtime.get("wheel_to_installed_tree_equivalence_proven") is not False:
        raise ScreenError("wheel-installed-tree-equivalence-enabled")
    for prefix in ("python", "cantera_wheel", "cantera_module", "cantera_extension"):
        if not isinstance(runtime.get(f"{prefix}_sha256"), str):
            raise ScreenError(f"runtime-sha256-missing:{prefix}")
        if not isinstance(runtime.get(f"{prefix}_size"), int) or runtime[f"{prefix}_size"] <= 0:
            raise ScreenError(f"runtime-size-invalid:{prefix}")
    policy = plan.get("thermo_envelope_policy", {})
    if policy.get("full_fixed_protocol_status") != "OUT_OF_MECHANISM_THERMO_ENVELOPE":
        raise ScreenError("full-protocol-status-open")
    if policy.get("prefix_diagnostic_status") != "blocked-valid-prefix-only-no-extrapolation":
        raise ScreenError("prefix-diagnostic-status-open")
    if policy.get("first_crossing_publishable_chemistry_state") is not False:
        raise ScreenError("crossing-state-publication-enabled")
    if policy.get("outside_samples_scored") is not False:
        raise ScreenError("outside-envelope-scoring-enabled")
    if policy.get("extrapolation_approved") is not False:
        raise ScreenError("extrapolation-enabled")
    return plan


def validate_reference(
    path: Path,
    model_id: str,
    partition: str,
    plan: dict[str, Any],
) -> dict[str, Any]:
    binding = REFERENCE_BINDINGS[(model_id, partition)]
    resolved = validate_bound_file(path, binding, f"reference:{model_id}:{partition}")
    reference = load_json(resolved, f"reference:{model_id}:{partition}")
    allowed_root = {
        "schema", "status", "model", "partition", "runtime_identity", "source",
        "reactor", "sample_times_s", "major_species", "prefix_diagnostic",
        "cases", "approvals",
    }
    if set(reference) != allowed_root:
        raise ScreenError(f"reference-root-schema-mismatch:{model_id}:{partition}")
    if reference["schema"] != "rfam.r2-c3mech-cantera-reference/v2":
        raise ScreenError("reference-schema-mismatch")
    if reference["status"] != "OUT_OF_MECHANISM_THERMO_ENVELOPE":
        raise ScreenError("reference-status-mismatch")
    if reference["partition"] != partition or reference["major_species"] != SPECIES:
        raise ScreenError("reference-partition-or-species-mismatch")
    if reference["runtime_identity"] != plan["runtime_identity"]:
        raise ScreenError("reference-runtime-identity-mismatch")
    if reference.get("prefix_diagnostic") != {
        "eligible_samples_only": True,
        "outside_endpoint_state_published": False,
        "later_samples_evaluated": False,
    }:
        raise ScreenError("reference-prefix-policy-mismatch")
    expected_model = model_from_plan(plan, model_id)
    model = reference["model"]
    for key in (
        "model_id", "role", "species_count", "author_header_reaction_count",
        "yaml_reaction_entry_count", "cantera_loaded_reaction_count",
        "thermo_min_k", "thermo_max_k",
    ):
        if model.get(key) != expected_model[key]:
            raise ScreenError(f"reference-model-mismatch:{model_id}:{key}")
    if model.get("count_discrepancy_resolved_as_equivalence") is not False:
        raise ScreenError("reference-reaction-count-equivalence-enabled")
    if model["yaml_reaction_entry_count"] != model["cantera_loaded_reaction_count"]:
        raise ScreenError("reference-execution-reaction-count-mismatch")
    source = reference["source"]
    if source.get("mechanism_sha256") != expected_model["sha256"]:
        raise ScreenError("reference-mechanism-sha256-mismatch")
    if source.get("generator_sha256") != plan["protocol"]["reference_generator_sha256"]:
        raise ScreenError("reference-generator-sha256-mismatch")
    manifest_key = "training_manifest_sha256" if partition == "protocol-training" else "holdout_manifest_sha256"
    if source.get("case_manifest_sha256") != plan["protocol"][manifest_key]:
        raise ScreenError("reference-case-manifest-sha256-mismatch")
    expected_cases = plan["protocol"][
        "training_case_ids" if partition == "protocol-training" else "holdout_case_ids"
    ]
    if reference["sample_times_s"] != plan["protocol"]["sample_times_s"]:
        raise ScreenError("reference-sample-times-mismatch")
    if [case.get("id") for case in reference["cases"]] != expected_cases:
        raise ScreenError("reference-case-order-mismatch")
    approvals = reference["approvals"]
    if not isinstance(approvals, dict) or any(
        type(value) is not bool or value for value in approvals.values()
    ):
        raise ScreenError("reference-approval-enabled")
    thermo_min = float(expected_model["thermo_min_k"])
    thermo_max = float(expected_model["thermo_max_k"])
    for case in reference["cases"]:
        validate_reference_case(case, partition, reference["sample_times_s"], thermo_min, thermo_max)
    reference["_bound_path"] = str(resolved)
    reference["_bound_sha256"] = binding["sha256"]
    return reference


def validate_reference_case(
    case: dict[str, Any],
    partition: str,
    sample_times: list[float],
    thermo_min: float,
    thermo_max: float,
) -> None:
    expected_case_keys = {
        "id", "temperature_k", "pressure_pa", "equivalence_ratio", "fuel",
        "oxidizer", "partition", "full_protocol_status",
        "initial_major_species_mass_fractions", "ignition_delay_definition",
        "ignition_delay", "thermo_envelope_exit", "samples",
    }
    if set(case) != expected_case_keys:
        raise ScreenError("case-root-schema-mismatch")
    if case.get("partition") != partition:
        raise ScreenError("case-partition-mismatch")
    if case.get("full_protocol_status") != "OUT_OF_MECHANISM_THERMO_ENVELOPE":
        raise ScreenError("case-full-protocol-status-mismatch")
    initial = case.get("initial_major_species_mass_fractions")
    if not isinstance(initial, dict) or set(initial) != set(SPECIES):
        raise ScreenError("case-initial-composition-schema-mismatch")
    initial_values = [finite(initial[name], f"{case['id']}.initial.{name}") for name in SPECIES]
    if any(value < 0.0 for value in initial_values) or abs(sum(initial_values) - 1.0) > 2.0e-9:
        raise ScreenError("case-initial-composition-invalid")
    samples = case.get("samples")
    if not isinstance(samples, list) or len(samples) != len(sample_times):
        raise ScreenError("case-sample-count-mismatch")
    statuses = [sample.get("status") if isinstance(sample, dict) else None for sample in samples]
    outside = [index for index, status in enumerate(statuses) if status == "outside-thermo-envelope-unscored"]
    if len(outside) != 1:
        raise ScreenError("case-envelope-exit-count-mismatch")
    exit_index = outside[0]
    if any(status != "eligible-within-mechanism-thermo-intersection" for status in statuses[:exit_index]):
        raise ScreenError("case-pre-exit-status-mismatch")
    if any(status != "not-evaluated-after-thermo-envelope-exit" for status in statuses[exit_index + 1:]):
        raise ScreenError("case-post-exit-status-mismatch")
    for index, (sample, expected_time) in enumerate(zip(samples, sample_times)):
        if finite(sample.get("time_s"), "sample-time") != float(expected_time):
            raise ScreenError("case-sample-time-mismatch")
        if index < exit_index:
            expected_keys = {
                "time_s", "status", "temperature_k", "pressure_pa", "density_kg_m3",
                "heat_release_rate_w_m3", "major_species_mass_fractions",
            }
            if set(sample) != expected_keys:
                raise ScreenError("eligible-sample-schema-mismatch")
            temperature = finite(sample["temperature_k"], "sample-temperature")
            if temperature < thermo_min or temperature > thermo_max:
                raise ScreenError("eligible-sample-outside-thermo-envelope")
            finite(sample["pressure_pa"], "sample-pressure")
            finite(sample["density_kg_m3"], "sample-density")
            finite(sample["heat_release_rate_w_m3"], "sample-qdot")
            composition = sample["major_species_mass_fractions"]
            if not isinstance(composition, dict) or set(composition) != set(SPECIES):
                raise ScreenError("eligible-sample-composition-schema-mismatch")
            for name in SPECIES:
                value = finite(composition[name], f"sample-Y.{name}")
                if value < 0.0 or value > 1.0:
                    raise ScreenError("eligible-sample-composition-out-of-range")
        elif index == exit_index:
            if set(sample) != {"time_s", "status", "exit_detection"}:
                raise ScreenError("outside-sample-schema-mismatch")
            detection = sample["exit_detection"]
            if detection != case.get("thermo_envelope_exit"):
                raise ScreenError("case-envelope-exit-rebind-mismatch")
            if detection.get("status") != "thermo-envelope-exit-detected-no-endpoint-state-published":
                raise ScreenError("case-envelope-exit-status-mismatch")
            if set(detection) != {
                "requested_sample_time_s", "thermo_min_k", "thermo_max_k",
                "first_crossing_estimate", "status",
            }:
                raise ScreenError("case-envelope-exit-schema-mismatch")
            crossing = detection.get("first_crossing_estimate")
            if not isinstance(crossing, dict) or set(crossing) != {
                "time_s", "temperature_k", "pressure_pa", "method",
                "publishable_chemistry_state",
            }:
                raise ScreenError("case-first-crossing-schema-mismatch")
            if crossing.get("method") != "10-K-advance-limited-linear-temperature-interpolation":
                raise ScreenError("case-first-crossing-method-mismatch")
            if crossing.get("publishable_chemistry_state") is not False:
                raise ScreenError("case-first-crossing-state-published")
            crossing_time = finite(crossing.get("time_s"), "crossing-time")
            if crossing_time <= 0.0 or crossing_time > float(expected_time):
                raise ScreenError("case-first-crossing-time-invalid")
            if finite(crossing.get("temperature_k"), "crossing-temperature") != thermo_max:
                raise ScreenError("case-first-crossing-temperature-mismatch")
            crossing_pressure = finite(crossing.get("pressure_pa"), "crossing-pressure")
            initial_pressure = finite(case.get("pressure_pa"), "case-pressure")
            if abs(crossing_pressure - initial_pressure) > max(1.0e-8 * initial_pressure, 1.0e-6):
                raise ScreenError("case-first-crossing-pressure-mismatch")
            if finite(detection.get("thermo_max_k"), "exit-thermo-max") != thermo_max:
                raise ScreenError("case-envelope-exit-bound-mismatch")
            if finite(detection.get("thermo_min_k"), "exit-thermo-min") != thermo_min:
                raise ScreenError("case-envelope-min-bound-mismatch")
            if finite(detection.get("requested_sample_time_s"), "exit-request-time") != float(expected_time):
                raise ScreenError("case-envelope-request-time-mismatch")
        elif set(sample) != {"time_s", "status"}:
            raise ScreenError("unevaluated-sample-schema-mismatch")
    ignition = case.get("ignition_delay")
    if not isinstance(ignition, dict):
        raise ScreenError("ignition-schema-mismatch")
    if ignition.get("status") == "eligible-within-mechanism-thermo-intersection":
        if set(ignition) != {
            "status", "ignition_delay_s", "threshold_temperature_k", "method"
        }:
            raise ScreenError("eligible-ignition-schema-mismatch")
        if finite(ignition["ignition_delay_s"], "ignition-delay") <= 0.0:
            raise ScreenError("ignition-delay-invalid")
        threshold = finite(ignition["threshold_temperature_k"], "ignition-threshold")
        if threshold < thermo_min or threshold >= thermo_max:
            raise ScreenError("eligible-ignition-threshold-outside-envelope")
        if ignition["method"] != "1-K-advance-limited-linear-threshold-interpolation":
            raise ScreenError("eligible-ignition-method-mismatch")
    elif ignition.get("status") == "threshold-at-or-outside-thermo-upper-bound-unscored":
        if "ignition_delay_s" in ignition:
            raise ScreenError("outside-ignition-delay-published")
        if set(ignition) != {"status", "threshold_temperature_k", "thermo_max_k"}:
            raise ScreenError("outside-ignition-schema-mismatch")
        if finite(ignition.get("threshold_temperature_k"), "ignition-threshold") < thermo_max:
            raise ScreenError("outside-ignition-threshold-inside")
        if finite(ignition.get("thermo_max_k"), "ignition-thermo-max") != thermo_max:
            raise ScreenError("outside-ignition-bound-mismatch")
    else:
        raise ScreenError("ignition-status-mismatch")


def to_wsl(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive
    if len(drive) != 2 or drive[1] != ":":
        raise ScreenError("unsupported-runner-path")
    relative = resolved.relative_to(Path(drive + "\\"))
    return f"/mnt/{drive[0].lower()}/{relative.as_posix()}"


def run_rfam_case(runner: Path, case: dict[str, Any], sample_times: list[float]) -> tuple[dict[str, str], list[dict[str, Any]]]:
    initial = case["initial_major_species_mass_fractions"]
    command = [
        "wsl", "-e", to_wsl(runner), case["id"],
        format(float(case["temperature_k"]), ".17g"),
        format(float(case["pressure_pa"]), ".17g"),
        str(len(sample_times)), "0",
        *[format(float(initial[name]), ".17g") for name in SPECIES],
        *[format(float(value), ".17g") for value in sample_times],
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise ScreenError(f"rfam-runner-rejected:{case['id']}:rc={result.returncode}")
    try:
        headers, samples = parse_runner_output(
            result.stdout, result.stderr, case["id"], len(sample_times)
        )
    except ProtocolError as exc:
        raise ScreenError(f"rfam-runner-output-invalid:{case['id']}:{exc}") from exc
    for index, (sample, expected_time) in enumerate(zip(samples, sample_times)):
        if sample["time_s"] != float(expected_time):
            raise ScreenError(f"rfam-sample-time-mismatch:{case['id']}:{index}")
        if not all(math.isfinite(value) for value in sample["mass_fractions"]):
            raise ScreenError(f"rfam-sample-composition-nonfinite:{case['id']}:{index}")
    return headers, samples


def compare_case(
    case: dict[str, Any],
    headers: dict[str, str],
    rfam_samples: list[dict[str, Any]],
    rfam_thermo_min: float,
    rfam_thermo_max: float,
) -> dict[str, Any]:
    detailed_samples = case["samples"]
    detailed_eligible_indices = [
        index for index, sample in enumerate(detailed_samples)
        if sample["status"] == "eligible-within-mechanism-thermo-intersection"
    ]
    eligible_indices = [
        index for index in detailed_eligible_indices
        if rfam_thermo_min <= float(rfam_samples[index]["temperature_k"]) <= rfam_thermo_max
    ]
    if not eligible_indices:
        raise ScreenError(f"no-paired-thermo-envelope-samples:{case['id']}")
    paired_noninitial_indices = [index for index in eligible_indices if index != 0]
    temperature_errors: list[float] = []
    species_l1_errors: list[float] = []
    for index in paired_noninitial_indices:
        detailed = detailed_samples[index]
        rfam = rfam_samples[index]
        temperature_errors.append(abs(rfam["temperature_k"] - float(detailed["temperature_k"])))
        species_l1_errors.append(sum(
            abs(rfam["mass_fractions"][species_index] - float(detailed["major_species_mass_fractions"][name]))
            for species_index, name in enumerate(SPECIES)
        ))
    detailed_qdot0 = float(detailed_samples[0]["heat_release_rate_w_m3"])
    rfam_qdot0 = float(rfam_samples[0]["qdot_w_m3"])
    ignition = case["ignition_delay"]
    ignition_result: dict[str, Any]
    if ignition["status"] == "eligible-within-mechanism-thermo-intersection":
        if headers.get("ignition_detected") != "1":
            raise ScreenError(f"rfam-ignition-not-detected:{case['id']}")
        detailed_idt = float(ignition["ignition_delay_s"])
        rfam_idt = float(headers["ignition_delay_s"])
        ignition_result = {
            "status": "eligible-within-detailed-mechanism-thermo-intersection",
            "detailed_ignition_delay_s": stable(detailed_idt, "detailed-idt"),
            "rfam_ignition_delay_s": stable(rfam_idt, "rfam-idt"),
            "rfam_to_detailed_ratio": stable(rfam_idt / detailed_idt, "idt-ratio"),
        }
    else:
        ignition_result = {
            "status": "unscored-detailed-threshold-at-thermo-upper-bound",
            "detailed_ignition_delay_s": None,
            "rfam_ignition_delay_s": stable(float(headers["ignition_delay_s"]), "rfam-idt"),
            "rfam_to_detailed_ratio": None,
        }
    last = eligible_indices[-1]
    last_noninitial = (
        paired_noninitial_indices[-1] if paired_noninitial_indices else None
    )
    return {
        "case_id": case["id"],
        "partition": case["partition"],
        "status": "blocked-valid-prefix-only-no-extrapolation",
        "full_fixed_protocol_status": case["full_protocol_status"],
        "eligible_sample_count": len(eligible_indices),
        "paired_noninitial_sample_count": len(paired_noninitial_indices),
        "detailed_mechanism_eligible_sample_count": len(detailed_eligible_indices),
        "unscored_paired_rfam_outside_count": len(detailed_eligible_indices) - len(eligible_indices),
        "outside_sample_count": 1,
        "not_evaluated_sample_count": len(detailed_samples) - len(detailed_eligible_indices) - 1,
        "last_eligible_sample_time_s": stable(float(detailed_samples[last]["time_s"]), "last-time"),
        "ignition": ignition_result,
        "maximum_temperature_absolute_error_k": (
            stable(max(temperature_errors), "temperature-error")
            if temperature_errors else None
        ),
        "maximum_major_species_l1": (
            stable(max(species_l1_errors), "species-l1")
            if species_l1_errors else None
        ),
        "last_eligible_temperature_delta_k": (
            stable(
                float(rfam_samples[last_noninitial]["temperature_k"])
                - float(detailed_samples[last_noninitial]["temperature_k"]),
                "last-temperature-delta",
            )
            if last_noninitial is not None else None
        ),
        "detailed_initial_qdot_w_m3": stable(detailed_qdot0, "detailed-qdot"),
        "rfam_initial_qdot_w_m3": stable(rfam_qdot0, "rfam-qdot"),
        "initial_qdot_absolute_ratio": stable(
            abs(rfam_qdot0) / max(abs(detailed_qdot0), 1.0), "qdot-ratio"
        ),
        "initial_qdot_sign_matches": math.copysign(1.0, rfam_qdot0) == math.copysign(1.0, detailed_qdot0),
    }


def compare_reference(
    reference: dict[str, Any],
    runner: Path,
    rfam_thermo_min: float,
    rfam_thermo_max: float,
) -> list[dict[str, Any]]:
    return [
        compare_case(
            case,
            *run_rfam_case(runner, case, reference["sample_times_s"]),
            rfam_thermo_min,
            rfam_thermo_max,
        )
        for case in reference["cases"]
    ]


def aggregate(comparisons: list[dict[str, Any]]) -> dict[str, Any]:
    ratios = [
        item["ignition"]["rfam_to_detailed_ratio"]
        for item in comparisons
        if item["ignition"]["rfam_to_detailed_ratio"] is not None
    ]
    temperature_errors = [
        item["maximum_temperature_absolute_error_k"]
        for item in comparisons
        if item["maximum_temperature_absolute_error_k"] is not None
    ]
    species_errors = [
        item["maximum_major_species_l1"]
        for item in comparisons
        if item["maximum_major_species_l1"] is not None
    ]
    noninitial_count = sum(
        item["paired_noninitial_sample_count"] for item in comparisons
    )
    return {
        "case_count": len(comparisons),
        "eligible_sample_count": sum(item["eligible_sample_count"] for item in comparisons),
        "noninitial_eligible_sample_count": noninitial_count,
        "trajectory_error_metrics_available": noninitial_count > 0,
        "detailed_mechanism_eligible_sample_count": sum(
            item["detailed_mechanism_eligible_sample_count"] for item in comparisons
        ),
        "unscored_paired_rfam_outside_count": sum(
            item["unscored_paired_rfam_outside_count"] for item in comparisons
        ),
        "outside_sample_count": sum(item["outside_sample_count"] for item in comparisons),
        "not_evaluated_sample_count": sum(item["not_evaluated_sample_count"] for item in comparisons),
        "scored_ignition_delay_count": len(ratios),
        "minimum_rfam_to_detailed_ignition_delay_ratio": stable(min(ratios), "min-idt-ratio") if ratios else None,
        "maximum_rfam_to_detailed_ignition_delay_ratio": stable(max(ratios), "max-idt-ratio") if ratios else None,
        "maximum_temperature_absolute_error_k": (
            stable(max(temperature_errors), "max-temperature-error")
            if temperature_errors else None
        ),
        "maximum_major_species_l1": (
            stable(max(species_errors), "max-species-l1")
            if species_errors else None
        ),
        "minimum_initial_qdot_absolute_ratio": stable(
            min(item["initial_qdot_absolute_ratio"] for item in comparisons), "min-qdot-ratio"
        ),
        "maximum_initial_qdot_absolute_ratio": stable(
            max(item["initial_qdot_absolute_ratio"] for item in comparisons), "max-qdot-ratio"
        ),
        "all_initial_qdot_signs_match": all(
            item["initial_qdot_sign_matches"] for item in comparisons
        ),
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    workspace = args.workspace_root.resolve()
    plan_path = args.plan.resolve()
    plan = validate_plan(plan_path)
    access_trace = ["plan-and-immutable-role-lock-validated"]
    runner = workspace / plan["rfam_runner"]["executable_path"]
    runner_source = workspace / plan["rfam_runner"]["source_path"]
    protocol_module = Path(__file__).resolve().with_name("run_r2_0d_rfam_protocol.py")
    if sha256(runner) != plan["rfam_runner"]["executable_sha256"] or runner.stat().st_size != plan["rfam_runner"]["executable_size"]:
        raise ScreenError("rfam-runner-identity-mismatch")
    if sha256(runner_source) != plan["rfam_runner"]["source_sha256"]:
        raise ScreenError("rfam-runner-source-identity-mismatch")
    if sha256(protocol_module) != PROTOCOL_MODULE_SHA256:
        raise ScreenError("protocol-module-identity-mismatch")

    training_references: dict[str, dict[str, Any]] = {}
    for model_id, path in (("44I0", args.training_44i0), ("3R00", args.training_3r00)):
        training_references[model_id] = validate_reference(
            path, model_id, "protocol-training", plan
        )
        access_trace.append(f"training-reference-{model_id}-validated")
    training_comparisons = {
        model_id: compare_reference(
            reference,
            runner,
            float(plan["rfam_runner"]["thermo_comparison_min_k"]),
            float(plan["rfam_runner"]["thermo_comparison_max_k"]),
        )
        for model_id, reference in training_references.items()
    }
    access_trace.append("all-training-comparisons-completed-without-role-change")
    access_trace.append("holdout-access-enabled-after-role-lock-and-training-completion")

    holdout_references: dict[str, dict[str, Any]] = {}
    for model_id, path in (("44I0", args.holdout_44i0), ("3R00", args.holdout_3r00)):
        holdout_references[model_id] = validate_reference(
            path, model_id, "protocol-holdout", plan
        )
        access_trace.append(f"holdout-reference-{model_id}-validated")
    holdout_comparisons = {
        model_id: compare_reference(
            reference,
            runner,
            float(plan["rfam_runner"]["thermo_comparison_min_k"]),
            float(plan["rfam_runner"]["thermo_comparison_max_k"]),
        )
        for model_id, reference in holdout_references.items()
    }
    access_trace.append("all-holdout-comparisons-completed-without-role-change")

    full_protocol_cases = []
    for partition, references in (
        ("protocol-training", training_references),
        ("protocol-holdout", holdout_references),
    ):
        for model_id in ("44I0", "3R00"):
            for case in references[model_id]["cases"]:
                full_protocol_cases.append({
                    "model_id": model_id,
                    "partition": partition,
                    "case_id": case["id"],
                    "status": case["full_protocol_status"],
                    "first_crossing_estimate": case["thermo_envelope_exit"][
                        "first_crossing_estimate"
                    ],
                })
    if any(
        item["status"] != "OUT_OF_MECHANISM_THERMO_ENVELOPE"
        for item in full_protocol_cases
    ):
        raise ScreenError("full-protocol-case-status-open")

    paired_noninitial_count = sum(
        item["paired_noninitial_sample_count"]
        for model in (training_comparisons, holdout_comparisons)
        for comparisons in model.values()
        for item in comparisons
    )
    if paired_noninitial_count != 0:
        raise ScreenError("unexpected-paired-noninitial-trajectory-samples")

    reference_sources = []
    for model_id in ("44I0", "3R00"):
        for partition, references in (
            ("protocol-training", training_references),
            ("protocol-holdout", holdout_references),
        ):
            reference = references[model_id]
            reference_sources.append({
                "model_id": model_id,
                "partition": partition,
                "path": reference.pop("_bound_path"),
                "sha256": reference.pop("_bound_sha256"),
            })
    return {
        "schema": "rfam.r2-c3mech-candidate-screen/v2",
        "status": "blocked-offline-valid-prefix-diagnostic-only",
        "source": {
            "plan_path": str(plan_path),
            "plan_sha256": sha256(plan_path),
            "screen_tool_path": str(Path(__file__).resolve()),
            "screen_tool_sha256": sha256(Path(__file__).resolve()),
            "protocol_module_path": str(protocol_module),
            "protocol_module_sha256": sha256(protocol_module),
            "rfam_runner_path": str(runner.resolve()),
            "rfam_runner_sha256": sha256(runner),
            "rfam_runner_source_path": str(runner_source.resolve()),
            "rfam_runner_source_sha256": sha256(runner_source),
            "references": reference_sources,
        },
        "role_lock": plan["role_lock"],
        "access_trace": access_trace,
        "reaction_count_semantics": plan["reaction_count_semantics"],
        "thermo_envelope_policy": plan["thermo_envelope_policy"],
        "full_fixed_protocol": {
            "status": "OUT_OF_MECHANISM_THERMO_ENVELOPE",
            "case_count": len(full_protocol_cases),
            "all_cases_outside_mechanism_thermo_envelope": True,
            "cases": full_protocol_cases,
        },
        "prefix_diagnostic": {
            "status": "blocked-valid-prefix-only-no-extrapolation",
            "results": {
                "44I0": {
                    "role": model_from_plan(plan, "44I0")["role"],
                    "training": training_comparisons["44I0"],
                    "training_aggregate": aggregate(training_comparisons["44I0"]),
                    "holdout": holdout_comparisons["44I0"],
                    "holdout_aggregate": aggregate(holdout_comparisons["44I0"]),
                },
                "3R00": {
                    "role": model_from_plan(plan, "3R00")["role"],
                    "training": training_comparisons["3R00"],
                    "training_aggregate": aggregate(training_comparisons["3R00"]),
                    "holdout": holdout_comparisons["3R00"],
                    "holdout_aggregate": aggregate(holdout_comparisons["3R00"]),
                },
            },
        },
        "cross_model_selection": {
            "ranking_computed": False,
            "winner_selected": False,
            "training_or_holdout_used_to_change_roles": False,
            "reason": "roles were frozen from the pre-trajectory rights/scope screen; this output is diagnostic only",
        },
        "evidence_boundary": {
            "paired_noninitial_trajectory_sample_count": paired_noninitial_count,
            "trajectory_shape_comparison_available": False,
            "trajectory_error_metrics": None,
            "ignition_delay_comparison_is_diagnostic_only": True,
            "initial_qdot_comparison_is_diagnostic_only": True,
            "reason": "RFAM leaves its accepted 300-3400 K thermo envelope before the first nonzero fixed sample time",
        },
        "approvals": plan["approvals"],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--training-44i0", type=Path, required=True)
    parser.add_argument("--training-3r00", type=Path, required=True)
    parser.add_argument("--holdout-44i0", type=Path, required=True)
    parser.add_argument("--holdout-3r00", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    output = args.output.resolve()
    temporary = output.with_name(output.name + ".tmp")
    try:
        document = run(args)
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary.write_text(
            json.dumps(document, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, output)
    except Exception as exc:
        try:
            temporary.unlink(missing_ok=True)
        except OSError:
            pass
        print(
            "R2_C3MECH_CANDIDATE_SCREEN_FAIL reason=" + str(exc).replace("\n", " | "),
            file=sys.stderr,
        )
        return 2
    print("R2_C3MECH_CANDIDATE_SCREEN_BLOCKED models=2 cases=12")
    print("fit_executed=0 physical_validation_approved=0 production_ready=0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
