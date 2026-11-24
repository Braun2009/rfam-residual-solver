#!/usr/bin/env python3
"""Fail-closed provenance gate for the offline C3Mech candidate screen."""

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
from typing import Any


class GateError(RuntimeError):
    pass


CONTRACT_SHA256 = "4ffbaa735277cab1c1d2257887692fba75db3e5703ea762a0a54718099633f29"
CONTRACT_SIZE = 10167
ARTIFACT_IDS = [
    "rights_report", "c3mech_license", "cantera_wheel", "mechanism_44i0",
    "mechanism_3r00",
    "reference_generator", "training_cases", "holdout_cases", "screen_plan",
    "reference_44i0_training", "reference_44i0_holdout",
    "reference_3r00_training", "reference_3r00_holdout", "screen_tool",
    "rfam_protocol_module", "rfam_runner_source", "rfam_runner_executable",
    "screen_output",
]
REFERENCE_JOBS = {
    "44I0-training": (
        "mechanism_44i0", "44I0", "training_cases", "protocol-training",
        "reference_44i0_training",
    ),
    "44I0-holdout": (
        "mechanism_44i0", "44I0", "holdout_cases", "protocol-holdout",
        "reference_44i0_holdout",
    ),
    "3R00-training": (
        "mechanism_3r00", "3R00", "training_cases", "protocol-training",
        "reference_3r00_training",
    ),
    "3R00-holdout": (
        "mechanism_3r00", "3R00", "holdout_cases", "protocol-holdout",
        "reference_3r00_holdout",
    ),
}
EXPECTED_ACCESS_TRACE = [
    "plan-and-immutable-role-lock-validated",
    "training-reference-44I0-validated",
    "training-reference-3R00-validated",
    "all-training-comparisons-completed-without-role-change",
    "holdout-access-enabled-after-role-lock-and-training-completion",
    "holdout-reference-44I0-validated",
    "holdout-reference-3R00-validated",
    "all-holdout-comparisons-completed-without-role-change",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise GateError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        result = json.loads(
            path.read_text(encoding="utf-8"), object_pairs_hook=strict_object
        )
    except GateError:
        raise
    except Exception as exc:
        raise GateError(f"invalid-json:{label}:{type(exc).__name__}") from exc
    if not isinstance(result, dict):
        raise GateError(f"not-an-object:{label}")
    return result


def resolve_artifact(root: Path, relative: str) -> Path:
    path_fragment = Path(relative)
    if path_fragment.is_absolute() or ".." in path_fragment.parts:
        raise GateError(f"invalid-artifact-path:{relative}")
    path = (root / path_fragment).resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise GateError(f"artifact-path-escapes-root:{relative}") from exc
    return path


def verify_contract(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.stat().st_size != CONTRACT_SIZE:
        raise GateError("contract-size-mismatch")
    if sha256(path) != CONTRACT_SHA256:
        raise GateError("contract-sha256-mismatch")
    contract = load_json(path, "contract")
    expected_root_keys = {
        "schema", "contract_version", "status", "scope", "artifacts",
        "runtime_identity", "partition_access", "mechanism_roles",
        "reaction_count_semantics", "evidence_boundary", "approvals",
        "limitations", "supersession",
    }
    if set(contract) != expected_root_keys:
        raise GateError("contract-root-schema-mismatch")
    if contract["schema"] != "rfam.r2-c3mech-candidate-screen-contract/v2":
        raise GateError("contract-schema-mismatch")
    if type(contract["contract_version"]) is not int or contract["contract_version"] != 2:
        raise GateError("contract-version-mismatch")
    if contract["status"] != "blocked-offline-valid-prefix-diagnostic-only":
        raise GateError("contract-status-mismatch")
    artifacts = contract["artifacts"]
    if not isinstance(artifacts, list) or [item.get("id") for item in artifacts] != ARTIFACT_IDS:
        raise GateError("contract-artifact-order-mismatch")
    approvals = contract["approvals"]
    if not isinstance(approvals, dict) or not approvals or any(
        type(value) is not bool or value for value in approvals.values()
    ):
        raise GateError("contract-approval-enabled")
    if contract["mechanism_roles"].get("selected_model_id") is not None:
        raise GateError("contract-mechanism-selected")
    if contract["partition_access"].get("selection_executed") is not False:
        raise GateError("contract-selection-enabled")
    if contract["partition_access"].get("fit_executed") is not False:
        raise GateError("contract-fit-enabled")
    if contract["reaction_count_semantics"].get("header_and_execution_counts_interchangeable") is not False:
        raise GateError("contract-reaction-count-equivalence-enabled")
    if contract["evidence_boundary"].get("outside_envelope_samples_scored") is not False:
        raise GateError("contract-outside-envelope-scoring-enabled")
    boundary = contract["evidence_boundary"]
    if boundary.get("full_fixed_protocol_status") != "OUT_OF_MECHANISM_THERMO_ENVELOPE":
        raise GateError("contract-full-protocol-status-open")
    if boundary.get("paired_noninitial_trajectory_sample_count") != 0:
        raise GateError("contract-noninitial-trajectory-claim")
    if boundary.get("trajectory_error_metrics") is not None:
        raise GateError("contract-trajectory-metrics-published")
    if boundary.get("trajectory_shape_comparison_available") is not False:
        raise GateError("contract-trajectory-shape-enabled")
    if boundary.get("first_crossing_publishable_chemistry_state") is not False:
        raise GateError("contract-crossing-state-publication-enabled")
    supersession = contract.get("supersession")
    if not isinstance(supersession, dict) or supersession.get("disposition") != (
        "v1-superseded-before-acceptance-do-not-review-or-concatenate"
    ):
        raise GateError("contract-supersession-mismatch")
    return contract


def verify_artifacts(root: Path, contract: dict[str, Any]) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for artifact in contract["artifacts"]:
        if set(artifact) != {"id", "path", "sha256", "size", "kind"}:
            raise GateError(f"artifact-schema-mismatch:{artifact.get('id')}")
        artifact_id = artifact["id"]
        path = resolve_artifact(root, artifact["path"])
        if not path.is_file():
            raise GateError(f"artifact-missing:{artifact_id}")
        if type(artifact["size"]) is not int or path.stat().st_size != artifact["size"]:
            raise GateError(f"artifact-size-mismatch:{artifact_id}")
        if sha256(path) != artifact["sha256"]:
            raise GateError(f"artifact-sha256-mismatch:{artifact_id}")
        paths[artifact_id] = path
    return paths


def verify_runtime(contract: dict[str, Any], paths: dict[str, Path]) -> None:
    runtime = contract["runtime_identity"]
    expected_keys = {
        "python_executable", "python_sha256", "python_size", "cantera_version",
        "cantera_wheel", "cantera_wheel_sha256", "cantera_wheel_size",
        "cantera_module", "cantera_module_sha256", "cantera_module_size",
        "cantera_extension", "cantera_extension_sha256", "cantera_extension_size",
        "execution_authority", "wheel_to_installed_tree_equivalence_proven",
    }
    if set(runtime) != expected_keys or runtime["cantera_version"] != "3.2.0":
        raise GateError("runtime-schema-or-version-mismatch")
    if runtime["execution_authority"] != "bound-installed-module-and-extension":
        raise GateError("runtime-execution-authority-mismatch")
    if runtime["wheel_to_installed_tree_equivalence_proven"] is not False:
        raise GateError("runtime-wheel-installed-tree-equivalence-enabled")
    checks = (
        ("python", Path(runtime["python_executable"]), runtime["python_sha256"], runtime["python_size"]),
        ("cantera-wheel", Path(runtime["cantera_wheel"]), runtime["cantera_wheel_sha256"], runtime["cantera_wheel_size"]),
        ("cantera-module", Path(runtime["cantera_module"]), runtime["cantera_module_sha256"], runtime["cantera_module_size"]),
        ("cantera-extension", Path(runtime["cantera_extension"]), runtime["cantera_extension_sha256"], runtime["cantera_extension_size"]),
    )
    for label, path, digest, size in checks:
        if not path.is_file() or path.stat().st_size != size or sha256(path) != digest:
            raise GateError(f"runtime-identity-mismatch:{label}")
    if Path(sys.executable).resolve() != Path(runtime["python_executable"]).resolve():
        raise GateError("validator-python-executable-mismatch")
    wheel_artifact = artifact_by_id(contract, "cantera_wheel")
    if (
        wheel_artifact["sha256"] != runtime["cantera_wheel_sha256"]
        or wheel_artifact["size"] != runtime["cantera_wheel_size"]
        or Path(runtime["cantera_wheel"]).resolve()
        != paths["cantera_wheel"].resolve()
    ):
        raise GateError("runtime-wheel-artifact-binding-mismatch")


def artifact_by_id(contract: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [item for item in contract["artifacts"] if item["id"] == artifact_id]
    if len(matches) != 1:
        raise GateError(f"artifact-id-not-unique:{artifact_id}")
    return matches[0]


def verify_plan(contract: dict[str, Any], paths: dict[str, Path]) -> None:
    plan = load_json(paths["screen_plan"], "screen-plan")
    if plan.get("schema") != "rfam.r2-c3mech-candidate-screen-plan/v2":
        raise GateError("plan-schema-mismatch")
    if plan.get("contract_version") != 2:
        raise GateError("plan-version-mismatch")
    if plan.get("status") != "blocked-offline-prefix-diagnostic-only":
        raise GateError("plan-status-mismatch")
    runtime = plan.get("runtime_identity")
    if not isinstance(runtime, dict):
        raise GateError("plan-runtime-schema-mismatch")
    if runtime.get("wheel_to_installed_tree_equivalence_proven") is not False:
        raise GateError("plan-wheel-installed-tree-equivalence-enabled")
    if runtime != contract["runtime_identity"]:
        raise GateError("plan-runtime-mismatch")
    protocol = plan.get("protocol", {})
    bindings = {
        "reference_generator_sha256": "reference_generator",
        "training_manifest_sha256": "training_cases",
        "holdout_manifest_sha256": "holdout_cases",
    }
    for field, artifact_id in bindings.items():
        if protocol.get(field) != artifact_by_id(contract, artifact_id)["sha256"]:
            raise GateError(f"plan-protocol-binding-mismatch:{field}")
    if protocol.get("reference_generator_size") != artifact_by_id(
        contract, "reference_generator"
    )["size"]:
        raise GateError("plan-generator-size-mismatch")
    role_lock = plan.get("role_lock", {})
    roles = contract["mechanism_roles"]
    if (
        role_lock.get("locked") is not True
        or role_lock.get("primary_model_id") != roles["primary_model_id"]
        or role_lock.get("comparator_model_id") != roles["comparator_model_id"]
        or role_lock.get("selected_model_id") is not None
        or role_lock.get("training_results_may_change_roles") is not False
        or role_lock.get("holdout_results_may_change_roles") is not False
    ):
        raise GateError("plan-role-lock-mismatch")
    mechanisms = plan.get("mechanisms")
    if not isinstance(mechanisms, list) or [item.get("model_id") for item in mechanisms] != ["44I0", "3R00"]:
        raise GateError("plan-mechanism-order-mismatch")
    for item, artifact_id in zip(mechanisms, ("mechanism_44i0", "mechanism_3r00")):
        artifact = artifact_by_id(contract, artifact_id)
        if item.get("sha256") != artifact["sha256"] or item.get("size") != artifact["size"]:
            raise GateError(f"plan-mechanism-binding-mismatch:{item.get('model_id')}")
        if item.get("yaml_reaction_entry_count") != item.get("cantera_loaded_reaction_count"):
            raise GateError(f"plan-execution-count-mismatch:{item.get('model_id')}")
    policy = plan.get("thermo_envelope_policy", {})
    if policy.get("full_fixed_protocol_status") != "OUT_OF_MECHANISM_THERMO_ENVELOPE":
        raise GateError("plan-full-protocol-status-open")
    if policy.get("prefix_diagnostic_status") != "blocked-valid-prefix-only-no-extrapolation":
        raise GateError("plan-prefix-diagnostic-status-open")
    if policy.get("first_crossing_publishable_chemistry_state") is not False:
        raise GateError("plan-crossing-state-publication-enabled")
    if policy.get("outside_samples_scored") is not False or policy.get("extrapolation_approved") is not False:
        raise GateError("plan-envelope-policy-open")


def verify_references(contract: dict[str, Any], paths: dict[str, Path]) -> None:
    runtime = dict(contract["runtime_identity"])
    expected = {
        "reference_44i0_training": ("44I0", "protocol-training", "mechanism_44i0", "training_cases"),
        "reference_44i0_holdout": ("44I0", "protocol-holdout", "mechanism_44i0", "holdout_cases"),
        "reference_3r00_training": ("3R00", "protocol-training", "mechanism_3r00", "training_cases"),
        "reference_3r00_holdout": ("3R00", "protocol-holdout", "mechanism_3r00", "holdout_cases"),
    }
    for artifact_id, (model_id, partition, mechanism_id, cases_id) in expected.items():
        reference = load_json(paths[artifact_id], artifact_id)
        if set(reference) != {
            "schema", "status", "model", "partition", "runtime_identity",
            "source", "reactor", "sample_times_s", "major_species",
            "prefix_diagnostic", "cases", "approvals",
        }:
            raise GateError(f"reference-root-schema:{artifact_id}")
        if reference.get("schema") != "rfam.r2-c3mech-cantera-reference/v2":
            raise GateError(f"reference-schema-mismatch:{artifact_id}")
        if reference.get("status") != "OUT_OF_MECHANISM_THERMO_ENVELOPE":
            raise GateError(f"reference-status-mismatch:{artifact_id}")
        if reference.get("partition") != partition or reference.get("model", {}).get("model_id") != model_id:
            raise GateError(f"reference-model-or-partition-mismatch:{artifact_id}")
        if reference.get("runtime_identity") != runtime:
            raise GateError(f"reference-runtime-mismatch:{artifact_id}")
        if reference.get("prefix_diagnostic") != {
            "eligible_samples_only": True,
            "later_samples_evaluated": False,
            "outside_endpoint_state_published": False,
        }:
            raise GateError(f"reference-prefix-policy:{artifact_id}")
        source = reference.get("source", {})
        if source.get("generator_sha256") != artifact_by_id(contract, "reference_generator")["sha256"]:
            raise GateError(f"reference-generator-mismatch:{artifact_id}")
        if source.get("mechanism_sha256") != artifact_by_id(contract, mechanism_id)["sha256"]:
            raise GateError(f"reference-mechanism-mismatch:{artifact_id}")
        if source.get("case_manifest_sha256") != artifact_by_id(contract, cases_id)["sha256"]:
            raise GateError(f"reference-cases-mismatch:{artifact_id}")
        approvals = reference.get("approvals")
        if not isinstance(approvals, dict) or any(
            type(value) is not bool or value for value in approvals.values()
        ):
            raise GateError(f"reference-approval-enabled:{artifact_id}")
        cases = reference.get("cases")
        if not isinstance(cases, list) or len(cases) != 3:
            raise GateError(f"reference-case-count:{artifact_id}")
        for case in cases:
            if case.get("full_protocol_status") != "OUT_OF_MECHANISM_THERMO_ENVELOPE":
                raise GateError(f"reference-full-protocol-status:{artifact_id}:{case.get('id')}")
            statuses = [sample.get("status") for sample in case.get("samples", [])]
            if statuses.count("outside-thermo-envelope-unscored") != 1:
                raise GateError(f"reference-envelope-exit-count:{artifact_id}:{case.get('id')}")
            exit_index = statuses.index("outside-thermo-envelope-unscored")
            if any(
                status != "not-evaluated-after-thermo-envelope-exit"
                for status in statuses[exit_index + 1:]
            ):
                raise GateError(f"reference-post-exit-publication:{artifact_id}:{case.get('id')}")
            for sample in case["samples"][exit_index + 1:]:
                if set(sample) != {"status", "time_s"}:
                    raise GateError(f"reference-post-exit-payload:{artifact_id}:{case.get('id')}")
            exit_sample = case["samples"][exit_index]
            if set(exit_sample) != {"status", "time_s", "exit_detection"}:
                raise GateError(f"reference-exit-sample-schema:{artifact_id}:{case.get('id')}")
            detection = exit_sample.get("exit_detection")
            if detection != case.get("thermo_envelope_exit"):
                raise GateError(f"reference-exit-rebind:{artifact_id}:{case.get('id')}")
            crossing = detection.get("first_crossing_estimate") if isinstance(detection, dict) else None
            if not isinstance(crossing, dict) or set(crossing) != {
                "method", "pressure_pa", "publishable_chemistry_state",
                "temperature_k", "time_s",
            }:
                raise GateError(f"reference-crossing-schema:{artifact_id}:{case.get('id')}")
            if crossing.get("method") != "10-K-advance-limited-linear-temperature-interpolation":
                raise GateError(f"reference-crossing-method:{artifact_id}:{case.get('id')}")
            if crossing.get("publishable_chemistry_state") is not False:
                raise GateError(f"reference-crossing-state-published:{artifact_id}:{case.get('id')}")
            for field in ("pressure_pa", "temperature_k", "time_s"):
                value = crossing.get(field)
                if (
                    isinstance(value, bool)
                    or not isinstance(value, (int, float))
                    or not math.isfinite(float(value))
                ):
                    raise GateError(f"reference-crossing-numeric:{artifact_id}:{case.get('id')}:{field}")
            ignition = case.get("ignition_delay", {})
            if ignition.get("status") == "eligible-within-mechanism-thermo-intersection":
                if ignition.get("method") != "1-K-advance-limited-linear-threshold-interpolation":
                    raise GateError(f"reference-ignition-method:{artifact_id}:{case.get('id')}")
            elif ignition.get("status") != "threshold-at-or-outside-thermo-upper-bound-unscored":
                raise GateError(f"reference-ignition-status:{artifact_id}:{case.get('id')}")


def verify_screen_output(contract: dict[str, Any], paths: dict[str, Path]) -> None:
    output = load_json(paths["screen_output"], "screen-output")
    expected_root = {
        "schema", "status", "source", "role_lock", "access_trace",
        "reaction_count_semantics", "thermo_envelope_policy",
        "full_fixed_protocol", "prefix_diagnostic", "cross_model_selection",
        "evidence_boundary", "approvals",
    }
    if set(output) != expected_root:
        raise GateError("screen-output-root-schema-mismatch")
    if output.get("schema") != "rfam.r2-c3mech-candidate-screen/v2":
        raise GateError("screen-output-schema-mismatch")
    if output.get("status") != "blocked-offline-valid-prefix-diagnostic-only":
        raise GateError("screen-output-status-mismatch")
    if output.get("access_trace") != EXPECTED_ACCESS_TRACE:
        raise GateError("screen-output-access-order-mismatch")
    plan = load_json(paths["screen_plan"], "screen-plan-for-output")
    for field in ("role_lock", "reaction_count_semantics", "thermo_envelope_policy"):
        if output.get(field) != plan.get(field):
            raise GateError(f"screen-output-plan-rebind:{field}")
    if output.get("cross_model_selection") != {
        "ranking_computed": False,
        "winner_selected": False,
        "training_or_holdout_used_to_change_roles": False,
        "reason": "roles were frozen from the pre-trajectory rights/scope screen; this output is diagnostic only",
    }:
        raise GateError("screen-output-selection-mismatch")
    boundary = output.get("evidence_boundary", {})
    if boundary.get("paired_noninitial_trajectory_sample_count") != 0:
        raise GateError("screen-output-noninitial-trajectory-claim")
    if boundary.get("trajectory_shape_comparison_available") is not False:
        raise GateError("screen-output-trajectory-shape-enabled")
    if boundary.get("trajectory_error_metrics") is not None:
        raise GateError("screen-output-trajectory-metrics-published")
    approvals = output.get("approvals")
    if approvals != contract["approvals"]:
        raise GateError("screen-output-approval-mismatch")
    source = output.get("source", {})
    bindings = {
        "plan_sha256": "screen_plan",
        "screen_tool_sha256": "screen_tool",
        "protocol_module_sha256": "rfam_protocol_module",
        "rfam_runner_sha256": "rfam_runner_executable",
        "rfam_runner_source_sha256": "rfam_runner_source",
    }
    for field, artifact_id in bindings.items():
        if source.get(field) != artifact_by_id(contract, artifact_id)["sha256"]:
            raise GateError(f"screen-output-source-mismatch:{field}")
    reference_pairs = [
        (item.get("model_id"), item.get("partition"), item.get("sha256"))
        for item in source.get("references", [])
    ]
    expected_pairs = [
        ("44I0", "protocol-training", artifact_by_id(contract, "reference_44i0_training")["sha256"]),
        ("44I0", "protocol-holdout", artifact_by_id(contract, "reference_44i0_holdout")["sha256"]),
        ("3R00", "protocol-training", artifact_by_id(contract, "reference_3r00_training")["sha256"]),
        ("3R00", "protocol-holdout", artifact_by_id(contract, "reference_3r00_holdout")["sha256"]),
    ]
    if reference_pairs != expected_pairs:
        raise GateError("screen-output-reference-binding-mismatch")

    expected_full_cases = []
    for partition, reference_ids in (
        ("protocol-training", ("reference_44i0_training", "reference_3r00_training")),
        ("protocol-holdout", ("reference_44i0_holdout", "reference_3r00_holdout")),
    ):
        for model_id, reference_id in zip(("44I0", "3R00"), reference_ids):
            reference = load_json(paths[reference_id], reference_id)
            for case in reference["cases"]:
                expected_full_cases.append({
                    "model_id": model_id,
                    "partition": partition,
                    "case_id": case["id"],
                    "status": case["full_protocol_status"],
                    "first_crossing_estimate": case["thermo_envelope_exit"][
                        "first_crossing_estimate"
                    ],
                })
    full_protocol = output.get("full_fixed_protocol")
    if full_protocol != {
        "status": "OUT_OF_MECHANISM_THERMO_ENVELOPE",
        "case_count": 12,
        "all_cases_outside_mechanism_thermo_envelope": True,
        "cases": expected_full_cases,
    }:
        raise GateError("screen-output-full-protocol-evidence-mismatch")
    for item in full_protocol["cases"]:
        if item["first_crossing_estimate"].get("publishable_chemistry_state") is not False:
            raise GateError("screen-output-crossing-state-published")

    prefix = output.get("prefix_diagnostic")
    if not isinstance(prefix, dict) or set(prefix) != {"status", "results"}:
        raise GateError("screen-output-prefix-schema-mismatch")
    if prefix.get("status") != "blocked-valid-prefix-only-no-extrapolation":
        raise GateError("screen-output-prefix-status-mismatch")
    results = prefix.get("results")
    if not isinstance(results, dict) or set(results) != {"44I0", "3R00"}:
        raise GateError("screen-output-results-schema-mismatch")
    for model_id, model in results.items():
        for partition in ("training", "holdout"):
            cases = model.get(partition)
            aggregate = model.get(partition + "_aggregate")
            if not isinstance(cases, list) or len(cases) != 3 or not isinstance(aggregate, dict):
                raise GateError(f"screen-output-partition-shape:{model_id}:{partition}")
            if aggregate.get("noninitial_eligible_sample_count") != 0:
                raise GateError(f"screen-output-noninitial-scored:{model_id}:{partition}")
            if aggregate.get("trajectory_error_metrics_available") is not False:
                raise GateError(f"screen-output-trajectory-metrics-available:{model_id}:{partition}")
            if aggregate.get("maximum_temperature_absolute_error_k") is not None:
                raise GateError(f"screen-output-temperature-metric-published:{model_id}:{partition}")
            if aggregate.get("maximum_major_species_l1") is not None:
                raise GateError(f"screen-output-species-metric-published:{model_id}:{partition}")
            if aggregate.get("all_initial_qdot_signs_match") is not False:
                raise GateError(f"screen-output-qdot-sign-claim:{model_id}:{partition}")
            if any(case.get("eligible_sample_count") != 1 for case in cases):
                raise GateError(f"screen-output-paired-sample-count:{model_id}:{partition}")
            for case in cases:
                if case.get("paired_noninitial_sample_count") != 0:
                    raise GateError(f"screen-output-case-noninitial:{model_id}:{partition}")
                for field in (
                    "maximum_temperature_absolute_error_k",
                    "maximum_major_species_l1",
                    "last_eligible_temperature_delta_k",
                ):
                    if case.get(field) is not None:
                        raise GateError(
                            f"screen-output-case-trajectory-metric:{model_id}:{partition}:{field}"
                        )


def run_checked(command: list[str], env: dict[str, str] | None, label: str) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, capture_output=True, text=True, check=False, env=env)
    if result.returncode != 0:
        detail = (result.stdout + " | " + result.stderr).strip().replace("\n", " | ")
        raise GateError(f"{label}-failed:rc={result.returncode}:{detail}")
    return result


def verify_screen_regeneration(root: Path, paths: dict[str, Path]) -> None:
    with tempfile.TemporaryDirectory(prefix="rfam-c3mech-screen-", dir=root / "tmp") as directory:
        output = Path(directory) / "screen.json"
        command = [
            sys.executable, str(paths["screen_tool"]),
            "--workspace-root", str(root),
            "--plan", str(paths["screen_plan"]),
            "--training-44i0", str(paths["reference_44i0_training"]),
            "--training-3r00", str(paths["reference_3r00_training"]),
            "--holdout-44i0", str(paths["reference_44i0_holdout"]),
            "--holdout-3r00", str(paths["reference_3r00_holdout"]),
            "--output", str(output),
        ]
        run_checked(command, None, "screen-regeneration")
        if output.read_bytes() != paths["screen_output"].read_bytes():
            raise GateError("screen-regeneration-byte-mismatch")


def verify_reference_regeneration(
    root: Path,
    contract: dict[str, Any],
    paths: dict[str, Path],
    job_name: str,
) -> None:
    mechanism_id, model_id, cases_id, partition, reference_id = REFERENCE_JOBS[job_name]
    runtime = contract["runtime_identity"]
    environment = dict(os.environ)
    package_root = str(Path(runtime["cantera_module"]).resolve().parents[1])
    environment["PYTHONPATH"] = package_root + os.pathsep + environment.get("PYTHONPATH", "")
    with tempfile.TemporaryDirectory(prefix="rfam-c3mech-reference-", dir=root / "tmp") as directory:
        output = Path(directory) / "reference.json"
        command = [
            runtime["python_executable"], str(paths["reference_generator"]),
            "--mechanism", str(paths[mechanism_id]),
            "--model-id", model_id,
            "--case-manifest", str(paths[cases_id]),
            "--partition", partition,
            "--output", str(output),
        ]
        run_checked(command, environment, f"reference-regeneration:{job_name}")
        if output.read_bytes() != paths[reference_id].read_bytes():
            raise GateError(f"reference-regeneration-byte-mismatch:{job_name}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument(
        "--contract", type=Path,
        default=Path("realFluidAllMach/tests/reference/r2_c3mech_candidate_screen_contract_v2.json"),
    )
    parser.add_argument("--verify-screen-regeneration", action="store_true")
    parser.add_argument(
        "--verify-reference-regeneration", action="append", choices=sorted(REFERENCE_JOBS), default=[]
    )
    args = parser.parse_args(argv)
    try:
        root = args.workspace_root.resolve()
        contract_path = args.contract if args.contract.is_absolute() else root / args.contract
        contract = verify_contract(contract_path.resolve())
        paths = verify_artifacts(root, contract)
        verify_runtime(contract, paths)
        verify_plan(contract, paths)
        verify_references(contract, paths)
        verify_screen_output(contract, paths)
        checks = len(paths) + 7
        if args.verify_screen_regeneration:
            verify_screen_regeneration(root, paths)
            checks += 1
        for job_name in args.verify_reference_regeneration:
            verify_reference_regeneration(root, contract, paths, job_name)
            checks += 1
    except Exception as exc:
        print(
            "R2_C3MECH_CANDIDATE_SCREEN_GATE_FAIL reason=" + str(exc).replace("\n", " | "),
            file=sys.stderr,
        )
        return 2
    print(f"R2_C3MECH_CANDIDATE_SCREEN_GATE_BLOCKED checks={checks}")
    print("fit_executed=0 physical_validation_approved=0 production_ready=0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
