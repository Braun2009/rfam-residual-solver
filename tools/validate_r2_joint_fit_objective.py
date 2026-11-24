#!/usr/bin/env python3
"""Validate and evaluate the frozen synthetic joint R1/R2 objective contract.

The successful state is still blocked.  This tool evaluates only a static,
repository-authored fixture and contains no mechanism runner or optimizer.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any


SCHEMA = "rfam.r2-joint-fit-objective-contract/v1"
FIXTURE_SCHEMA = "rfam.r2-joint-fit-synthetic-fixture/v1"
BLOCKED_MARKER = "R2_JOINT_FIT_OBJECTIVE_SYNTHETIC_ONLY_BLOCKED"
FAIL_MARKER = "R2_JOINT_FIT_OBJECTIVE_FAIL"

EXPECTED_CONTRACT_SHA256 = (
    "4e15ed2c7cca709b5bd54f04faa6cd2db50347e31dacc9676213b3fe375679e7"
)
EXPECTED_ARTIFACT = {
    "id": "synthetic_objective_fixture",
    "path": "tests/reference/r2_joint_fit_synthetic_fixture_v1.json",
    "sha256": "6cef95b7fd3d56eee4a009214a5205d539c69681712e3088be6b7418a1b2199c",
    "byte_count": 6652,
    "kind": "repository-authored-synthetic-nonphysical-objective-wiring-fixture",
}
EXPECTED_APPROVALS = {
    "detailed_mechanism_selected_for_fit": False,
    "experimental_training_data_selected": False,
    "experimental_holdout_data_selected": False,
    "physical_right_to_fit_resolved": False,
    "physical_right_to_redistribute_resolved": False,
    "optimizer_execution_approved": False,
    "parameter_publication_approved": False,
    "physical_validation_approved": False,
    "hem_chemistry_enabled": False,
    "fluent_execution_approved": False,
    "reacting_hot_start_approved": False,
    "production_ready": False,
}
EXPECTED_STATE_SPECIES = ["C3H8", "O2", "CO", "CO2", "H2O", "N2"]
EXPECTED_OBSERVABLE_SPECIES = ["C3H8", "O2", "CO", "CO2", "H2O"]
EXPECTED_PARAMETERS = [
    "R1_A_SI",
    "R1_Ea_J_per_mol",
    "R1_order_C3H8",
    "R1_order_O2",
    "R2_A_SI",
    "R2_Ea_J_per_mol",
]
EXPECTED_OBSERVABLES = [
    "ignition_delay", "qdot", "temperature", "C3H8", "O2", "CO", "CO2", "H2O"
]
EXPECTED_CASE_ORDER = [
    "SYNTHETIC_TRAIN_001",
    "SYNTHETIC_TRAIN_002",
    "SYNTHETIC_HOLDOUT_001",
    "SYNTHETIC_HOLDOUT_002",
]


class ObjectiveContractError(ValueError):
    """A fail-closed contract violation."""


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ObjectiveContractError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def _reject_nonfinite_json(token: str) -> None:
    raise ObjectiveContractError(f"nonfinite-json-number:{token}")


def _load_json(path: Path, label: str) -> tuple[dict[str, Any], bytes]:
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise ObjectiveContractError(f"artifact-open-failed:{label}") from exc
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ObjectiveContractError(f"json-not-utf8:{label}") from exc
    try:
        value = json.loads(
            text,
            object_pairs_hook=_reject_duplicate_pairs,
            parse_constant=_reject_nonfinite_json,
        )
    except ObjectiveContractError:
        raise
    except json.JSONDecodeError as exc:
        raise ObjectiveContractError(f"json-invalid:{label}:{exc.msg}") from exc
    if not isinstance(value, dict):
        raise ObjectiveContractError(f"not-object:{label}")
    return value, raw


def _sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _object(value: Any, label: str, keys: set[str]) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ObjectiveContractError(f"not-object:{label}")
    missing = keys - set(value)
    unknown = set(value) - keys
    if missing:
        raise ObjectiveContractError(f"missing-keys:{label}:" + ",".join(sorted(missing)))
    if unknown:
        raise ObjectiveContractError(f"unknown-keys:{label}:" + ",".join(sorted(unknown)))
    return value


def _array(value: Any, label: str) -> list[Any]:
    if not isinstance(value, list):
        raise ObjectiveContractError(f"not-array:{label}")
    return value


def _finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ObjectiveContractError(f"invalid-number:{label}")
    result = float(value)
    if not math.isfinite(result):
        raise ObjectiveContractError(f"nonfinite-number:{label}")
    return result


def _exact(actual: Any, expected: Any, label: str) -> None:
    if type(actual) is not type(expected) or actual != expected:
        raise ObjectiveContractError(f"contract-mismatch:{label}")


def _resolve_under(root: Path, relative: Any, label: str) -> Path:
    if (
        not isinstance(relative, str)
        or not relative
        or "\\" in relative
        or relative.startswith("/")
        or re.match(r"^[A-Za-z]:", relative)
    ):
        raise ObjectiveContractError(f"invalid-relative-path:{label}")
    root = root.resolve()
    path = (root / relative).resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise ObjectiveContractError(f"path-escapes-root:{label}") from exc
    return path


def _verify_contract(document: dict[str, Any], raw: bytes) -> None:
    _object(
        document,
        "root",
        {
            "schema", "contract_version", "status", "scope", "approvals", "artifact",
            "mechanism_adapter", "domain", "partitions", "observables_and_basis",
            "metrics_weights_tolerances", "parameter_bounds", "optimizer_protocol",
            "identifiability_uncertainty", "provenance_rights", "synthetic_expectation",
            "cross_scope", "stop_conditions",
        },
    )
    _exact(document["schema"], SCHEMA, "root.schema")
    if type(document["contract_version"]) is not int or document["contract_version"] != 1:
        raise ObjectiveContractError("invalid-int:root.contract_version")
    _exact(document["status"], "blocked-synthetic-objective-wiring-only", "root.status")
    _exact(
        document["scope"],
        "mechanism-agnostic-joint-r1-r2-objective-schema-no-fit-no-physical-validation",
        "root.scope",
    )
    approvals = _object(document["approvals"], "approvals", set(EXPECTED_APPROVALS))
    for key, expected in EXPECTED_APPROVALS.items():
        if type(approvals[key]) is not bool or approvals[key] is not expected:
            raise ObjectiveContractError(f"approval-must-remain-false:{key}")
    if _sha256(raw) != EXPECTED_CONTRACT_SHA256:
        raise ObjectiveContractError("contract-sha256-contract-mismatch")


def _verify_static_sections(document: dict[str, Any]) -> None:
    adapter = _object(
        document["mechanism_adapter"],
        "mechanism_adapter",
        {
            "reference_backend", "candidate_backend", "reactor_model", "species_order",
            "reaction_topology", "required_runner_inputs", "required_runner_outputs",
            "runner_failure_policy", "posthoc_reference_or_partition_changes_allowed",
        },
    )
    _exact(adapter["species_order"], EXPECTED_STATE_SPECIES,
           "mechanism_adapter.species_order")
    _exact(adapter["posthoc_reference_or_partition_changes_allowed"], False,
           "mechanism_adapter.posthoc_reference_or_partition_changes_allowed")
    topology = _object(
        adapter["reaction_topology"], "mechanism_adapter.reaction_topology",
        {"R1", "R2", "R2_reverse_rate", "R2_forward_orders",
         "independently_fitted_reverse_parameters"},
    )
    _exact(topology["R2_forward_orders"], {"CO": 1.0, "O2": 0.5},
           "mechanism_adapter.reaction_topology.R2_forward_orders")
    _exact(topology["independently_fitted_reverse_parameters"], False,
           "mechanism_adapter.reaction_topology.independent_reverse")

    domain = _object(
        document["domain"], "domain",
        {"phase_scope", "initial_temperature_k", "initial_pressure_pa", "equivalence_ratio",
         "diluent_mass_fraction", "initial_composition", "sample_time_s",
         "out_of_domain_policy", "rfam_hem_owner"},
    )
    _exact(domain["phase_scope"], "gas-phase-only-no-vle-no-quality-no-solid", "domain.phase_scope")
    _exact(domain["out_of_domain_policy"], "fail-closed-no-extrapolation", "domain.out_of_domain_policy")
    composition = _object(
        domain["initial_composition"], "domain.initial_composition",
        {"basis", "species_order", "sum_tolerance", "diluent_species",
         "diluent_definition", "equivalence_ratio_definition",
         "runner_state_authority", "molecular_weights_kg_per_kmol",
         "metadata_consistency_absolute_tolerance"},
    )
    _exact(composition["basis"], "mass-fraction", "domain.initial_composition.basis")
    _exact(composition["species_order"], EXPECTED_STATE_SPECIES,
           "domain.initial_composition.species_order")
    _exact(composition["diluent_species"], "N2",
           "domain.initial_composition.diluent_species")
    _exact(composition["runner_state_authority"],
           "ordered-initial-mass-fraction-vector-no-reconstruction-or-override-from-metadata",
           "domain.initial_composition.runner_state_authority")

    partitions = _object(
        document["partitions"], "partitions",
        {"training_case_ids", "holdout_case_ids", "disjoint_required", "frozen_before_selection",
         "holdout_objective_calls_during_parameter_selection", "holdout_role"},
    )
    training_ids = _array(partitions["training_case_ids"], "partitions.training_case_ids")
    holdout_ids = _array(partitions["holdout_case_ids"], "partitions.holdout_case_ids")
    _exact(training_ids, EXPECTED_CASE_ORDER[:2], "partitions.training_case_ids")
    _exact(holdout_ids, EXPECTED_CASE_ORDER[2:], "partitions.holdout_case_ids")
    if set(training_ids) & set(holdout_ids):
        raise ObjectiveContractError("training-holdout-overlap")
    _exact(partitions["disjoint_required"], True, "partitions.disjoint_required")
    _exact(partitions["frozen_before_selection"], True, "partitions.frozen_before_selection")
    _exact(partitions["holdout_objective_calls_during_parameter_selection"], 0,
           "partitions.holdout_objective_calls_during_parameter_selection")

    basis = _object(
        document["observables_and_basis"], "observables_and_basis",
        {"observable_order", "ignition_delay", "qdot", "temperature", "species", "energy_basis"},
    )
    _exact(basis["observable_order"], EXPECTED_OBSERVABLES,
           "observables_and_basis.observable_order")
    _exact(basis["species"]["names"], EXPECTED_OBSERVABLE_SPECIES,
           "observables_and_basis.species.names")
    _exact(basis["energy_basis"]["id"],
           "formation-inclusive-total-species-enthalpy-298.15K-v1",
           "observables_and_basis.energy_basis.id")
    _exact(basis["energy_basis"]["common_reference_and_prediction_basis_required"], True,
           "observables_and_basis.energy_basis.common_reference_and_prediction_basis_required")
    _exact(basis["energy_basis"]["basis_offset_fitting_allowed"], False,
           "observables_and_basis.energy_basis.basis_offset_fitting_allowed")

    parameters = _object(
        document["parameter_bounds"], "parameter_bounds",
        {"rate_definition", "concentration_unit", "A_unit_rule", "fixed_parameters",
         "ordered_fit_parameters", "bounds_status"},
    )
    _exact(parameters["concentration_unit"], "mol/m3", "parameter_bounds.concentration_unit")
    entries = _array(parameters["ordered_fit_parameters"], "parameter_bounds.ordered_fit_parameters")
    if [entry.get("name") if isinstance(entry, dict) else None for entry in entries] != EXPECTED_PARAMETERS:
        raise ObjectiveContractError("parameter-order-mismatch")
    if len(entries) != 6:
        raise ObjectiveContractError("parameter-count-mismatch")
    for index, entry in enumerate(entries):
        item = _object(entry, f"parameter_bounds.ordered_fit_parameters[{index}]",
                       {"name", "reaction", "role", "lower", "upper", "transform", "units"})
        lower = _finite(item["lower"], f"parameter[{index}].lower")
        upper = _finite(item["upper"], f"parameter[{index}].upper")
        if not lower < upper:
            raise ObjectiveContractError(f"parameter-bounds-not-ordered:{item['name']}")
        if item["transform"] == "log10" and lower <= 0.0:
            raise ObjectiveContractError(f"log-parameter-lower-not-positive:{item['name']}")
        if item["role"] == "activation-energy":
            _exact(item["units"], "J/mol", f"parameter[{index}].activation_energy_units")
    fixed = _object(
        parameters["fixed_parameters"], "parameter_bounds.fixed_parameters",
        {"R1_beta", "R2_beta", "R2_order_CO", "R2_order_O2", "R2_reverse_policy"},
    )
    _exact(fixed["R2_order_CO"], 1.0, "parameter_bounds.fixed_parameters.R2_order_CO")
    _exact(fixed["R2_order_O2"], 0.5, "parameter_bounds.fixed_parameters.R2_order_O2")

    optimizer = _object(
        document["optimizer_protocol"], "optimizer_protocol",
        {"status", "selection_split", "coordinate_space", "start_sequence", "search",
         "failed_candidate_policy", "holdout_evaluation", "optimizer_executed"},
    )
    _exact(optimizer["status"], "specified-not-executed", "optimizer_protocol.status")
    _exact(optimizer["selection_split"], "training-only", "optimizer_protocol.selection_split")
    start_sequence = _object(
        optimizer["start_sequence"], "optimizer_protocol.start_sequence",
        {"first", "remaining", "halton_indices", "halton_radical_inverse_index_origin",
         "halton_definition", "midpoint_is_separate_from_halton_sequence",
         "bases_in_parameter_order", "transformed_bound_mapping", "random_seed"},
    )
    search = _object(
        optimizer["search"], "optimizer_protocol.search",
        {"algorithm", "coordinate_order", "direction_order", "initial_normalized_step",
         "step_contraction", "minimum_normalized_step",
         "maximum_objective_evaluations_per_start", "parallel_objective_evaluations",
         "candidate_generation", "out_of_bound_candidate_policy",
         "coordinate_acceptance", "sweep_and_contraction",
         "post_acceptance_progression", "evaluation_count", "cache_scope",
         "evaluation_limit_policy", "termination", "improvement_rule",
         "invalid_start_policy", "global_start_selection", "final_tie_break"},
    )
    _exact(search["coordinate_order"], EXPECTED_PARAMETERS,
           "optimizer_protocol.search.coordinate_order")
    _exact(search["parallel_objective_evaluations"], 1,
           "optimizer_protocol.search.parallel_objective_evaluations")
    _exact(start_sequence["halton_radical_inverse_index_origin"], 1,
           "optimizer_protocol.start_sequence.halton_index_origin")
    _exact(start_sequence["halton_indices"], [1, 16],
           "optimizer_protocol.start_sequence.halton_indices")
    _exact(start_sequence["midpoint_is_separate_from_halton_sequence"], True,
           "optimizer_protocol.start_sequence.midpoint_is_separate")
    _exact(search["out_of_bound_candidate_policy"],
           "skip-without-runner-call-or-evaluation-count",
           "optimizer_protocol.search.out_of_bound_candidate_policy")
    _exact(search["cache_scope"], "per-start-cleared-before-each-start",
           "optimizer_protocol.search.cache_scope")
    _exact(optimizer["optimizer_executed"], False, "optimizer_protocol.optimizer_executed")

    ident = _object(
        document["identifiability_uncertainty"], "identifiability_uncertainty",
        {"status", "scaled_residual_order", "scaled_residual_definition", "local_sensitivity",
         "identifiability_requirements", "uncertainty_requirements"},
    )
    _exact(ident["status"], "not-evaluated-blocking", "identifiability_uncertainty.status")
    residual_definition = _object(
        ident["scaled_residual_definition"], "identifiability.scaled_residual_definition",
        {"ignition_delay_row", "history_row", "identity", "holdout_rows_excluded"},
    )
    _exact(residual_definition["holdout_rows_excluded"], True,
           "identifiability.scaled_residual_definition.holdout_rows_excluded")
    sensitivity = _object(
        ident["local_sensitivity"], "identifiability.local_sensitivity",
        {"method", "step", "one_sided_at_bounds", "boundary_rule",
         "failed_base_or_perturbation_policy"},
    )
    _exact(sensitivity["one_sided_at_bounds"], True,
           "identifiability.local_sensitivity.one_sided_at_bounds")
    requirements = _object(
        ident["identifiability_requirements"], "identifiability.requirements",
        {"svd_relative_rank_tolerance", "singular_value_order", "rank_definition",
         "required_rank", "condition_number_definition", "maximum_condition_number",
         "parameter_correlation_definition", "maximum_absolute_parameter_correlation",
         "rank_deficiency_policy"},
    )
    uncertainty = _object(
        ident["uncertainty_requirements"], "identifiability.uncertainty_requirements",
        {"confidence_level", "covariance_definition", "residual_variance_definition",
         "rank_or_variance_failure_policy", "required_outputs",
         "synthetic-fixture-uncertainty-is-physical-evidence"},
    )
    _exact(requirements["required_rank"], 6,
           "identifiability_uncertainty.identifiability_requirements.required_rank")
    _exact(uncertainty["synthetic-fixture-uncertainty-is-physical-evidence"],
           False,
           "identifiability_uncertainty.uncertainty_requirements.synthetic-evidence")

    rights = _object(
        document["provenance_rights"], "provenance_rights",
        {"input_class", "synthetic_fixture_local_contract_test_use",
         "physical_mechanism_artifact_bound", "physical_experimental_artifact_bound",
         "physical_fit_use_permission_bound", "physical_redistribution_permission_bound",
         "required_before_physical_fit",
         "synthetic_fixture_must_not_be_relabelled_as_reference_mechanism_or_experiment"},
    )
    _exact(rights["input_class"], "SYNTHETIC_NONPHYSICAL_TEST_ONLY",
           "provenance_rights.input_class")
    for key in (
        "physical_mechanism_artifact_bound", "physical_experimental_artifact_bound",
        "physical_fit_use_permission_bound", "physical_redistribution_permission_bound",
    ):
        _exact(rights[key], False, f"provenance_rights.{key}")
    _exact(rights["synthetic_fixture_must_not_be_relabelled_as_reference_mechanism_or_experiment"],
           True, "provenance_rights.synthetic_fixture_relabel_policy")


def _domain_interval(domain: dict[str, Any], key: str) -> tuple[float, float]:
    values = _array(domain[key], f"domain.{key}")
    if len(values) != 2:
        raise ObjectiveContractError(f"domain-interval-length:{key}")
    lower = _finite(values[0], f"domain.{key}[0]")
    upper = _finite(values[1], f"domain.{key}[1]")
    if lower > upper:
        raise ObjectiveContractError(f"domain-interval-not-ordered:{key}")
    return lower, upper


def _history(value: Any, label: str, count: int) -> list[float]:
    items = _array(value, label)
    if len(items) != count:
        raise ObjectiveContractError(f"history-length-mismatch:{label}")
    return [_finite(item, f"{label}[{index}]") for index, item in enumerate(items)]


def _verify_fixture_static(fixture: dict[str, Any], document: dict[str, Any]) -> None:
    _object(
        fixture, "fixture",
        {"schema", "watermark", "status", "origin", "units", "energy_basis_id",
         "species_order", "parameter_order", "synthetic_trial_parameter_vector",
         "prediction_origin", "cases"},
    )
    _exact(fixture["schema"], FIXTURE_SCHEMA, "fixture.schema")
    _exact(fixture["watermark"], "SYNTHETIC_NONPHYSICAL_TEST_ONLY", "fixture.watermark")
    _exact(fixture["status"],
           "static-objective-wiring-fixture-no-mechanism-evaluation-no-fit",
           "fixture.status")
    _exact(fixture["origin"], {
        "kind": "repository-authored-static-synthetic-fixture",
        "physical_measurements_present": False,
        "detailed_mechanism_output_present": False,
        "cantera_output_present": False,
        "rfam_model_output_present": False,
        "fit_executed": False,
    }, "fixture.origin")
    _exact(fixture["units"], {
        "time": "s", "temperature": "K", "pressure": "Pa",
        "heat_release_rate": "W/m3", "species_composition": "mass-fraction",
    }, "fixture.units")
    _exact(fixture["energy_basis_id"],
           document["observables_and_basis"]["energy_basis"]["id"],
           "fixture.energy_basis_id")
    _exact(fixture["species_order"], EXPECTED_STATE_SPECIES, "fixture.species_order")
    _exact(fixture["parameter_order"], EXPECTED_PARAMETERS, "fixture.parameter_order")
    _exact(fixture["prediction_origin"],
           "hand-authored-objective-wiring-values-not-derived-from-parameter-vector",
           "fixture.prediction_origin")

    vector = _array(fixture["synthetic_trial_parameter_vector"],
                    "fixture.synthetic_trial_parameter_vector")
    entries = document["parameter_bounds"]["ordered_fit_parameters"]
    if len(vector) != len(entries):
        raise ObjectiveContractError("synthetic-parameter-vector-length-mismatch")
    for index, (value, entry) in enumerate(zip(vector, entries)):
        number = _finite(value, f"fixture.synthetic_trial_parameter_vector[{index}]")
        if number < float(entry["lower"]) or number > float(entry["upper"]):
            raise ObjectiveContractError(f"synthetic-parameter-out-of-bounds:{entry['name']}")


def _evaluate_fixture(fixture: dict[str, Any], document: dict[str, Any]) -> tuple[dict[str, float], dict[str, float], int]:
    domain = document["domain"]
    bounds = {
        "initial_temperature_k": _domain_interval(domain, "initial_temperature_k"),
        "initial_pressure_pa": _domain_interval(domain, "initial_pressure_pa"),
        "equivalence_ratio": _domain_interval(domain, "equivalence_ratio"),
        "diluent_mass_fraction": _domain_interval(domain, "diluent_mass_fraction"),
        "sample_time_s": _domain_interval(domain, "sample_time_s"),
    }
    metric = document["metrics_weights_tolerances"]
    scales = metric["normalization_scales"]
    weights = metric["component_weights"]
    if abs(math.fsum(_finite(weights[key], f"weight.{key}") for key in EXPECTED_OBSERVABLES) - 1.0) > 1e-15:
        raise ObjectiveContractError("component-weights-do-not-sum-to-one")

    cases = _array(fixture["cases"], "fixture.cases")
    if len(cases) != len(EXPECTED_CASE_ORDER):
        raise ObjectiveContractError("fixture-case-count-mismatch")
    case_objectives: dict[str, float] = {}
    split_values: dict[str, list[float]] = {"training": [], "holdout": []}

    for case_index, case_value in enumerate(cases):
        case = _object(
            case_value, f"fixture.cases[{case_index}]",
            {"case_id", "split", "condition", "sample_times_s", "reference", "prediction"},
        )
        case_id = case["case_id"]
        _exact(case_id, EXPECTED_CASE_ORDER[case_index], f"fixture.cases[{case_index}].case_id")
        expected_split = "training" if case_index < 2 else "holdout"
        _exact(case["split"], expected_split, f"fixture.cases[{case_index}].split")

        condition = _object(
            case["condition"], f"fixture.cases[{case_index}].condition",
            {"initial_temperature_k", "initial_pressure_pa", "equivalence_ratio",
             "diluent_mass_fraction", "initial_mass_fractions"},
        )
        for key in ("initial_temperature_k", "initial_pressure_pa", "equivalence_ratio",
                    "diluent_mass_fraction"):
            number = _finite(condition[key], f"fixture.{case_id}.condition.{key}")
            lower, upper = bounds[key]
            if number < lower or number > upper:
                raise ObjectiveContractError(f"fixture-condition-out-of-domain:{case_id}:{key}")
        initial_y = _history(
            condition["initial_mass_fractions"],
            f"fixture.{case_id}.condition.initial_mass_fractions",
            len(EXPECTED_STATE_SPECIES),
        )
        if any(value < 0.0 or value > 1.0 for value in initial_y):
            raise ObjectiveContractError(f"fixture-initial-mass-fraction-out-of-range:{case_id}")
        composition_contract = domain["initial_composition"]
        sum_tolerance = _finite(composition_contract["sum_tolerance"],
                                "domain.initial_composition.sum_tolerance")
        if sum_tolerance <= 0.0 or abs(math.fsum(initial_y) - 1.0) > sum_tolerance:
            raise ObjectiveContractError(f"fixture-initial-mass-fraction-sum:{case_id}")
        metadata_tolerance = _finite(
            composition_contract["metadata_consistency_absolute_tolerance"],
            "domain.initial_composition.metadata_consistency_absolute_tolerance",
        )
        n2_index = EXPECTED_STATE_SPECIES.index("N2")
        if abs(initial_y[n2_index] - float(condition["diluent_mass_fraction"])) > metadata_tolerance:
            raise ObjectiveContractError(f"fixture-diluent-metadata-mismatch:{case_id}")
        fuel_index = EXPECTED_STATE_SPECIES.index("C3H8")
        oxygen_index = EXPECTED_STATE_SPECIES.index("O2")
        if initial_y[oxygen_index] <= 0.0:
            raise ObjectiveContractError(f"fixture-initial-oxygen-not-positive:{case_id}")
        molecular_weights = composition_contract["molecular_weights_kg_per_kmol"]
        stoichiometric_mass_ratio = (
            _finite(molecular_weights["C3H8"], "domain.initial_composition.mw.C3H8")
            / (3.5 * _finite(molecular_weights["O2"], "domain.initial_composition.mw.O2"))
        )
        computed_phi = (initial_y[fuel_index] / initial_y[oxygen_index]) / stoichiometric_mass_ratio
        if abs(computed_phi - float(condition["equivalence_ratio"])) > metadata_tolerance:
            raise ObjectiveContractError(f"fixture-equivalence-ratio-mismatch:{case_id}")

        times = _history(case["sample_times_s"], f"fixture.{case_id}.sample_times_s",
                         len(_array(case["sample_times_s"], f"fixture.{case_id}.sample_times_s")))
        if len(times) < 2 or any(right <= left for left, right in zip(times, times[1:])):
            raise ObjectiveContractError(f"fixture-time-grid-not-strictly-increasing:{case_id}")
        if times[0] < bounds["sample_time_s"][0] or times[-1] > bounds["sample_time_s"][1]:
            raise ObjectiveContractError(f"fixture-time-grid-out-of-domain:{case_id}")

        trajectories: dict[str, dict[str, Any]] = {}
        for role in ("reference", "prediction"):
            trajectory = _object(
                case[role], f"fixture.{case_id}.{role}",
                {"ignition_delay_s", "qdot_w_m3", "temperature_k", "mass_fractions"},
            )
            ignition = _finite(trajectory["ignition_delay_s"], f"fixture.{case_id}.{role}.ignition")
            if ignition <= 0.0 or ignition > bounds["sample_time_s"][1]:
                raise ObjectiveContractError(f"fixture-ignition-out-of-domain:{case_id}:{role}")
            qdot = _history(trajectory["qdot_w_m3"], f"fixture.{case_id}.{role}.qdot", len(times))
            temperature = _history(trajectory["temperature_k"], f"fixture.{case_id}.{role}.temperature", len(times))
            if temperature[0] != float(condition["initial_temperature_k"]):
                raise ObjectiveContractError(f"fixture-initial-temperature-mismatch:{case_id}:{role}")
            if any(value <= 0.0 for value in temperature):
                raise ObjectiveContractError(f"fixture-temperature-not-positive:{case_id}:{role}")
            species = _object(trajectory["mass_fractions"], f"fixture.{case_id}.{role}.mass_fractions",
                              set(EXPECTED_OBSERVABLE_SPECIES))
            histories = {name: _history(species[name], f"fixture.{case_id}.{role}.{name}", len(times))
                         for name in EXPECTED_OBSERVABLE_SPECIES}
            for name, values in histories.items():
                if any(value < 0.0 or value > 1.0 for value in values):
                    raise ObjectiveContractError(f"fixture-mass-fraction-out-of-range:{case_id}:{role}:{name}")
            for sample in range(len(times)):
                if math.fsum(histories[name][sample] for name in EXPECTED_OBSERVABLE_SPECIES) > 1.0 + 1e-12:
                    raise ObjectiveContractError(f"fixture-selected-mass-fraction-sum:{case_id}:{role}:{sample}")
            for species_index, name in enumerate(EXPECTED_OBSERVABLE_SPECIES):
                if abs(histories[name][0] - initial_y[species_index]) > metadata_tolerance:
                    raise ObjectiveContractError(
                        f"fixture-initial-trajectory-composition-mismatch:{case_id}:{role}:{name}"
                    )
            trajectories[role] = {
                "ignition_delay": ignition,
                "qdot": qdot,
                "temperature": temperature,
                **histories,
            }

        reference = trajectories["reference"]
        prediction = trajectories["prediction"]
        components: dict[str, float] = {
            "ignition_delay": (
                math.log(prediction["ignition_delay"] / reference["ignition_delay"])
                / _finite(scales["ignition_delay_sigma_log"], "scale.ignition_delay_sigma_log")
            ) ** 2
        }
        scale_keys = {
            "qdot": "qdot_w_m3",
            "temperature": "temperature_k",
            "C3H8": "C3H8_mass_fraction",
            "O2": "O2_mass_fraction",
            "CO": "CO_mass_fraction",
            "CO2": "CO2_mass_fraction",
            "H2O": "H2O_mass_fraction",
        }
        for observable, scale_key in scale_keys.items():
            scale = _finite(scales[scale_key], f"scale.{scale_key}")
            if scale <= 0.0:
                raise ObjectiveContractError(f"normalization-scale-not-positive:{scale_key}")
            components[observable] = math.fsum(
                ((predicted - expected) / scale) ** 2
                for predicted, expected in zip(prediction[observable], reference[observable])
            ) / len(times)
        objective = math.fsum(
            _finite(weights[name], f"weight.{name}") * components[name]
            for name in EXPECTED_OBSERVABLES
        )
        if not math.isfinite(objective) or objective < 0.0:
            raise ObjectiveContractError(f"invalid-case-objective:{case_id}")
        case_objectives[case_id] = objective
        split_values[expected_split].append(objective)

    split_objectives = {
        split: math.fsum(values) / len(values) for split, values in split_values.items()
    }
    return split_objectives, case_objectives, len(cases)


def validate(root: Path, contract_path: Path) -> tuple[int, dict[str, float]]:
    root = root.resolve()
    if not root.is_dir():
        raise ObjectiveContractError("root-not-directory")
    document, contract_raw = _load_json(contract_path, "contract")
    _verify_contract(document, contract_raw)
    _verify_static_sections(document)

    artifact = _object(document["artifact"], "artifact", set(EXPECTED_ARTIFACT))
    for key in ("id", "path", "sha256", "byte_count", "kind"):
        _exact(artifact[key], EXPECTED_ARTIFACT[key], f"artifact.{key}")
    fixture_path = _resolve_under(root, artifact["path"], "artifact.path")
    fixture, fixture_raw = _load_json(fixture_path, "synthetic-fixture")
    if len(fixture_raw) != EXPECTED_ARTIFACT["byte_count"]:
        raise ObjectiveContractError("artifact-byte-count-mismatch:synthetic_objective_fixture")
    if _sha256(fixture_raw) != EXPECTED_ARTIFACT["sha256"]:
        raise ObjectiveContractError("artifact-sha256-mismatch:synthetic_objective_fixture")

    _verify_fixture_static(fixture, document)
    split_objectives, case_objectives, case_count = _evaluate_fixture(fixture, document)
    expected = _object(
        document["synthetic_expectation"], "synthetic_expectation",
        {"training_cases", "holdout_cases", "training_objective", "holdout_objective",
         "case_objectives", "fit_executed", "physical_validation_claimed",
         "parameter_vector_publishable"},
    )
    tolerance = _finite(
        document["metrics_weights_tolerances"]["synthetic_regression_tolerances"]
        ["recomputed_objective_absolute"],
        "synthetic-regression-tolerance",
    )
    if expected["training_cases"] != 2 or expected["holdout_cases"] != 2:
        raise ObjectiveContractError("synthetic-expected-split-count-mismatch")
    for split in ("training", "holdout"):
        expected_value = _finite(expected[f"{split}_objective"], f"expected.{split}_objective")
        if abs(split_objectives[split] - expected_value) > tolerance:
            raise ObjectiveContractError(f"synthetic-{split}-objective-mismatch")
        maximum = _finite(
            document["metrics_weights_tolerances"]["synthetic_regression_tolerances"]
            [f"{split}_objective_maximum"],
            f"synthetic-{split}-objective-maximum",
        )
        if split_objectives[split] > maximum:
            raise ObjectiveContractError(f"synthetic-{split}-objective-over-limit")
    expected_cases = _object(expected["case_objectives"], "synthetic_expectation.case_objectives",
                             set(EXPECTED_CASE_ORDER))
    for case_id in EXPECTED_CASE_ORDER:
        if abs(case_objectives[case_id] - _finite(expected_cases[case_id], f"expected.{case_id}")) > tolerance:
            raise ObjectiveContractError(f"synthetic-case-objective-mismatch:{case_id}")
    for key in ("fit_executed", "physical_validation_claimed", "parameter_vector_publishable"):
        _exact(expected[key], False, f"synthetic_expectation.{key}")
    return 25 + case_count * 12, split_objectives


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--contract", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    contract = args.contract or root / "tests/reference/r2_joint_fit_objective_contract_v1.json"
    if not contract.is_absolute():
        contract = root / contract
    try:
        checks, objectives = validate(root, contract)
    except ObjectiveContractError as exc:
        print(f"{FAIL_MARKER} {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # Malformed inputs must never emit a traceback.
        print(f"{FAIL_MARKER} unexpected:{type(exc).__name__}", file=sys.stderr)
        return 2
    print(
        f"{BLOCKED_MARKER} synthetic_only=1 training_cases=2 holdout_cases=2 "
        f"training_objective={objectives['training']:.17g} "
        f"holdout_objective={objectives['holdout']:.17g} "
        f"fit_executed=0 physical_validation_approved=0 production_ready=0 "
        f"checks={checks}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
