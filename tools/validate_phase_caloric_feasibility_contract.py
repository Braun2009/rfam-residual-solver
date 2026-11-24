#!/usr/bin/env python3
"""Fail-closed validator for the offline phase-local caloric feasibility contract.

This tool validates a schema only. It imports no RFAM thermo implementation,
does not recover a phase, and never authorizes two-phase chemistry.
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


SCHEMA = "rfam.phase-caloric-feasibility-contract/v1"
BLOCKED_MARKER = "PHASE_CALORIC_FEASIBILITY_CONTRACT_BLOCKED"
FAIL_MARKER = "PHASE_CALORIC_FEASIBILITY_CONTRACT_FAIL"
EXPECTED_CONTRACT_SHA256 = (
    "fbef1eed655b28f5dc786856c3ec3e96b72765761422ecf0a4bc30ed9122853d"
)
EXPECTED_APPROVALS = {
    "phase_caloric_implementation_approved": False,
    "two_phase_chemistry_enabled": False,
    "phase_aware_chemistry_enabled": False,
    "equilibrium_result_admission_relaxed": False,
    "solver_consumer_admission": False,
    "checkpoint_runtime_integration": False,
    "restart_runtime_integration": False,
    "physical_validation_approved": False,
    "fluent_execution_approved": False,
    "reacting_hot_start_approved": False,
    "production_ready": False,
}
EXPECTED_RESULT_FORBIDDEN = [
    "quality", "phase_fraction", "vapor_volume_fraction",
    "hem_acoustic_closure", "sound_speed", "global_stability",
    "generic_solver_admission",
]
EXPECTED_RESULT_REQUIRED = [
    "pressure_pa", "temperature_k", "compressibility_z", "root_identity",
    "density_kg_m3", "specific_internal_energy_j_kg",
    "specific_enthalpy_j_kg", "mass_fractions", "cv_j_kg_k",
    "dedt_rho_j_kg_k", "dpdt_rho_pa_k", "dpdrho_t_pa_m3_kg", "model_id",
    "ordered_species", "energy_convention", "phase_caloric_schema_id",
    "phase_caloric_options_fingerprint", "parent_generation",
    "parent_phase_identity", "continuation_token", "residuals",
]
EXPECTED_NEGATIVE_CASES = [
    "liquid-payload-under-vapor-request",
    "vapor-supercritical-identity-mismatch",
    "seed-forced-alternate-root",
    "root-disappearance-or-spinodal-crossing",
    "near-critical-root-ambiguity",
    "composition-normalization-drift",
    "energy-reference-or-departure-caloric-drift",
    "parent-q-or-generation-mismatch",
    "token-reuse-across-rollback-restart-or-final-reflash",
    "initial-mid-rhs-or-final-failure-exact-rollback",
    "final-whole-mixture-reflash-enters-unsupported-three-phase-topology",
    "generic-solver-conversion-or-admission-attempt",
    "concurrent-or-reentrant-nondeterminism",
]
EXPECTED_ALLOWED_FILES = [
    "tests/reference/phase_caloric_feasibility_contract_v1.json",
    "tools/validate_phase_caloric_feasibility_contract.py",
    "tests/test_phase_caloric_feasibility_contract.py",
    "docs/phase_caloric_feasibility_contract.md",
]


class ContractError(ValueError):
    """A structured fail-closed contract error."""


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ContractError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def _reject_nonfinite_json(token: str) -> None:
    raise ContractError(f"nonfinite-json-number:{token}")


def _load(path: Path) -> tuple[dict[str, Any], bytes]:
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise ContractError("contract-open-failed") from exc
    try:
        value = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=_reject_duplicate_pairs,
            parse_constant=_reject_nonfinite_json,
        )
    except ContractError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ContractError("contract-json-invalid") from exc
    if not isinstance(value, dict):
        raise ContractError("contract-not-object")
    return value, raw


def _object(value: Any, label: str, keys: set[str]) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ContractError(f"not-object:{label}")
    missing = keys - set(value)
    unknown = set(value) - keys
    if missing:
        raise ContractError(f"missing-keys:{label}:" + ",".join(sorted(missing)))
    if unknown:
        raise ContractError(f"unknown-keys:{label}:" + ",".join(sorted(unknown)))
    return value


def _exact(actual: Any, expected: Any, label: str) -> None:
    if type(actual) is not type(expected) or actual != expected:
        raise ContractError(f"contract-mismatch:{label}")


def _finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"invalid-number:{label}")
    result = float(value)
    if not math.isfinite(result):
        raise ContractError(f"nonfinite-number:{label}")
    return result


def _resolve_under(root: Path, path: Any) -> Path:
    if (not isinstance(path, str) or not path or "\\" in path
            or path.startswith("/") or re.match(r"^[A-Za-z]:", path)):
        raise ContractError("invalid-contract-path")
    resolved_root = root.resolve()
    resolved = (resolved_root / path).resolve()
    try:
        resolved.relative_to(resolved_root)
    except ValueError as exc:
        raise ContractError("contract-path-escapes-root") from exc
    return resolved


def _verify(document: dict[str, Any], raw: bytes) -> int:
    _object(document, "root", {
        "schema", "contract_version", "status", "scope", "approvals",
        "phase_caloric_schema", "immutable_request_binding",
        "phase_local_validity", "chemistry_substep",
        "conservative_writeback_and_publication", "checkpoint_restart",
        "negative_cases_required", "evidence_boundary", "stop_conditions",
    })
    _exact(document["schema"], SCHEMA, "root.schema")
    _exact(document["contract_version"], 1, "root.contract_version")
    _exact(document["status"], "blocked-offline-phase-caloric-feasibility-contract",
           "root.status")
    _exact(document["scope"],
           "phase-local-vapor-or-supercritical-caloric-feasibility-no-runtime-implementation-no-two-phase-chemistry-success",
           "root.scope")
    approvals = _object(document["approvals"], "approvals", set(EXPECTED_APPROVALS))
    for key, expected in EXPECTED_APPROVALS.items():
        _exact(approvals[key], expected, f"approval-must-remain-false:{key}")
    if hashlib.sha256(raw).hexdigest() != EXPECTED_CONTRACT_SHA256:
        raise ContractError("contract-sha256-mismatch")

    schema = _object(document["phase_caloric_schema"], "phase_caloric_schema", {
        "id", "result_type", "result_is_equilibrium_state",
        "conversion_to_equilibrium_state", "forbidden_result_fields",
        "required_result_fields",
    })
    _exact(schema["id"], "rfam.phase-caloric/v1", "phase-caloric-schema-id")
    _exact(schema["result_type"], "PhaseCaloricResult", "result-type")
    _exact(schema["result_is_equilibrium_state"], False, "result-is-equilibrium-state")
    _exact(schema["conversion_to_equilibrium_state"], "forbidden", "result-conversion")
    _exact(schema["forbidden_result_fields"], EXPECTED_RESULT_FORBIDDEN,
           "result-forbidden-fields")
    _exact(schema["required_result_fields"], EXPECTED_RESULT_REQUIRED,
           "result-required-fields")

    binding = _object(document["immutable_request_binding"], "immutable_request_binding", {
        "required_parent_fields", "required_phase_identities",
        "automatic_phase_transition", "automatic_vapor_to_supercritical",
        "automatic_vapor_to_liquid", "multi_root_fallback", "first_recovery",
        "continuation",
    })
    _exact(binding["required_phase_identities"], ["VAPOR", "SUPERCRITICAL"],
           "required-phase-identities")
    for key in ("automatic_phase_transition", "automatic_vapor_to_supercritical",
                "automatic_vapor_to_liquid", "multi_root_fallback"):
        _exact(binding[key], "forbidden", f"branch-policy:{key}")
    first = _object(binding["first_recovery"], "first_recovery", {
        "must_reproduce_parent_payload", "pressure_relative_residual_max",
        "temperature_relative_residual_max", "density_relative_residual_max",
        "specific_internal_energy_relative_residual_max",
        "mass_fraction_l1_residual_max", "root_identity_must_match",
    })
    _exact(first["must_reproduce_parent_payload"], True, "first-recovery-parent-replay")
    _exact(first["root_identity_must_match"], True, "first-recovery-root-match")
    for key in ("pressure_relative_residual_max", "temperature_relative_residual_max",
                "density_relative_residual_max", "specific_internal_energy_relative_residual_max",
                "mass_fraction_l1_residual_max"):
        if _finite(first[key], f"first-recovery.{key}") <= 0.0:
            raise ContractError(f"nonpositive-residual-limit:{key}")
    continuation = _object(binding["continuation"], "continuation", {
        "seed_role", "previous_accepted_root_token_required_after_first",
        "root_continuity_required", "token_revalidation_against_parent_required",
        "parent_generation_mismatch", "parent_q_fingerprint_mismatch",
        "model_schema_options_species_or_energy_mismatch",
    })
    _exact(continuation["seed_role"], "numerical-hint-only-never-branch-authority",
           "seed-role")
    for key in ("previous_accepted_root_token_required_after_first",
                "root_continuity_required", "token_revalidation_against_parent_required"):
        _exact(continuation[key], True, f"continuation:{key}")
    for key in ("parent_generation_mismatch", "parent_q_fingerprint_mismatch",
                "model_schema_options_species_or_energy_mismatch"):
        _exact(continuation[key], "fail-closed", f"continuation:{key}")

    validity = _object(document["phase_local_validity"], "phase_local_validity", {
        "positive_cv_required", "positive_dedt_rho_required", "positive_dpdrho_t_required",
        "finite_pressure_temperature_density_energy_derivatives_required",
        "composition_sum_absolute_tolerance", "composition_nonnegative_required",
        "departure_caloric_consistency_required", "root_loss", "spinodal_crossing",
        "near_critical_ambiguity", "metastable_vapor_statement",
    })
    for key in ("positive_cv_required", "positive_dedt_rho_required",
                "positive_dpdrho_t_required",
                "finite_pressure_temperature_density_energy_derivatives_required",
                "composition_nonnegative_required", "departure_caloric_consistency_required"):
        _exact(validity[key], True, f"phase-validity:{key}")
    if _finite(validity["composition_sum_absolute_tolerance"], "composition-sum-tolerance") <= 0.0:
        raise ContractError("nonpositive-composition-sum-tolerance")
    for key in ("root_loss", "spinodal_crossing", "near_critical_ambiguity"):
        _exact(validity[key], "fail-closed", f"phase-validity:{key}")

    substep = _object(document["chemistry_substep"], "chemistry_substep", {
        "basis", "constant_pressure_chemistry", "phase_density_relative_displacement_max",
        "specific_internal_energy_relative_displacement_max",
        "mass_fraction_l1_displacement_max", "envelope_exit",
        "initial_mid_rhs_final_recovery_required", "concurrent_reentrant_determinism_required",
    })
    _exact(substep["basis"], "isochoric-fixed-phase-density-vapor-or-supercritical-chemistry",
           "chemistry-basis")
    _exact(substep["constant_pressure_chemistry"], "forbidden", "constant-pressure")
    _exact(substep["phase_density_relative_displacement_max"], 0.0, "fixed-density")
    for key in ("specific_internal_energy_relative_displacement_max",
                "mass_fraction_l1_displacement_max"):
        if _finite(substep[key], f"substep:{key}") <= 0.0:
            raise ContractError(f"nonpositive-substep-limit:{key}")
    _exact(substep["envelope_exit"], "fail-closed-with-exact-rollback", "envelope-exit")
    _exact(substep["initial_mid_rhs_final_recovery_required"], True,
           "substep-initial-mid-final")
    _exact(substep["concurrent_reentrant_determinism_required"], True,
           "substep-determinism")

    publication = _object(document["conservative_writeback_and_publication"],
                          "conservative_writeback_and_publication", {
        "vapor_volume_to_mixture_volume_rhoY_and_rhoE_required",
        "post_flash_radical_or_energy_fix", "only_publisher",
        "final_reflash_requirements", "three_phase_or_unsupported_topology",
    })
    _exact(publication["vapor_volume_to_mixture_volume_rhoY_and_rhoE_required"], True,
           "mixture-volume-writeback")
    _exact(publication["post_flash_radical_or_energy_fix"], "forbidden",
           "post-flash-fix")
    _exact(publication["only_publisher"],
           "final-whole-mixture-EquilibriumThermoModel-recoverRhoE", "only-publisher")
    _exact(publication["three_phase_or_unsupported_topology"], "fail-closed",
           "unsupported-topology")

    restart = _object(document["checkpoint_restart"], "checkpoint_restart", {
        "generic_equilibrium_state_serialization", "continuation_token_persistence",
        "restored_token_requires_parent_revalidation_before_rhs",
        "token_reuse_across_rollback_restart_or_final_reflash",
    })
    _exact(restart["generic_equilibrium_state_serialization"], "forbidden",
           "restart-generic-state")
    _exact(restart["continuation_token_persistence"], "opaque-hint-only",
           "restart-token-persistence")
    _exact(restart["restored_token_requires_parent_revalidation_before_rhs"], True,
           "restart-token-revalidation")
    _exact(restart["token_reuse_across_rollback_restart_or_final_reflash"],
           "fail-closed", "restart-token-reuse")

    _exact(document["negative_cases_required"], EXPECTED_NEGATIVE_CASES,
           "negative-case-matrix")
    boundary = _object(document["evidence_boundary"], "evidence_boundary", {
        "artifact_class", "physical_pre_validation", "runtime_implementation_present",
        "two_phase_chemistry_success_claimed", "source_files_allowed",
    })
    _exact(boundary["artifact_class"], "OFFLINE_SCHEMA_AND_TEST_CONTRACT_ONLY",
           "artifact-class")
    for key in ("physical_pre_validation", "runtime_implementation_present",
                "two_phase_chemistry_success_claimed"):
        _exact(boundary[key], False, f"evidence-boundary:{key}")
    _exact(boundary["source_files_allowed"], EXPECTED_ALLOWED_FILES,
           "allowed-source-files")
    if not isinstance(document["stop_conditions"], list) or not document["stop_conditions"]:
        raise ContractError("missing-stop-conditions")
    return 49


def validate(root: Path, contract_path: Path) -> int:
    resolved_root = root.resolve()
    if not resolved_root.is_dir():
        raise ContractError("root-not-directory")
    resolved_contract = (
        contract_path.resolve() if contract_path.is_absolute()
        else (resolved_root / contract_path).resolve()
    )
    document, raw = _load(resolved_contract)
    return _verify(document, raw)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--contract", type=Path,
                        default=Path("tests/reference/phase_caloric_feasibility_contract_v1.json"))
    args = parser.parse_args()
    root = args.root.resolve()
    contract = args.contract if args.contract.is_absolute() else root / args.contract
    try:
        checks = validate(root, contract)
    except ContractError as exc:
        print(f"{FAIL_MARKER} {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # Defensive fail-close without a traceback.
        print(f"{FAIL_MARKER} unexpected:{type(exc).__name__}", file=sys.stderr)
        return 2
    print(
        f"{BLOCKED_MARKER} offline_schema=1 runtime_implementation=0 "
        f"two_phase_chemistry_enabled=0 phase_aware_chemistry_enabled=0 "
        f"production_ready=0 checks={checks}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
