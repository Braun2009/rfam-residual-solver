#!/usr/bin/env python3
"""Fail-closed validator for the offline phase-local caloric v2 contract.

This validates a static feasibility schema only. It imports no RFAM runtime
thermo code, does not recover a phase, and cannot authorize chemistry.
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


SCHEMA = "rfam.phase-caloric-feasibility-contract/v2"
BLOCKED_MARKER = "PHASE_CALORIC_FEASIBILITY_CONTRACT_V2_BLOCKED"
FAIL_MARKER = "PHASE_CALORIC_FEASIBILITY_CONTRACT_V2_FAIL"
EXPECTED_CONTRACT_SHA256 = "7e33a853310b4af2f0c846576fc0c3f83ee68bc9f2fa50f0978034df3d5ff3f4"

APPROVALS = {
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
FORBIDDEN_RESULT_FIELDS = [
    "quality", "phase_fraction", "vapor_volume_fraction",
    "hem_acoustic_closure", "sound_speed", "global_stability",
    "generic_solver_admission",
]
SUCCESS_RESULT_FIELDS = [
    "pressure_pa", "temperature_k", "compressibility_z", "root_identity",
    "density_kg_m3", "specific_internal_energy_j_kg",
    "specific_enthalpy_j_kg", "mass_fractions", "cv_j_kg_k",
    "dedt_rho_j_kg_k", "dpdt_rho_pa_k", "dpdrho_t_pa_m3_kg",
    "model_id", "ordered_species", "energy_convention",
    "phase_caloric_schema_id", "phase_caloric_options_fingerprint",
    "parent_q_fingerprint", "parent_phase_payload_fingerprint",
    "parent_generation", "runtime_epoch", "accepted_root_token", "residuals",
]
FAILURE_FORBIDDEN_FIELDS = [
    "pressure_pa", "temperature_k", "compressibility_z", "root_identity",
    "density_kg_m3", "specific_internal_energy_j_kg",
    "specific_enthalpy_j_kg", "mass_fractions", "cv_j_kg_k",
    "dedt_rho_j_kg_k", "dpdt_rho_pa_k", "dpdrho_t_pa_m3_kg",
    "accepted_root_token",
]
FAILURE_DIAGNOSTICS = [
    "status", "failure_stage", "parent_binding_digest", "finite_residual_summary",
]
PARENT_FIELDS = [
    "committed_parent_q_fingerprint", "parent_phase_payload_fingerprint",
    "parent_generation", "runtime_epoch", "model_id", "phase_caloric_schema_id",
    "phase_caloric_options_fingerprint", "ordered_species", "energy_convention",
    "parent_phase_identity", "parent_phase_pressure_pa",
    "parent_phase_temperature_k", "parent_phase_compressibility_z",
    "parent_phase_density_kg_m3", "parent_phase_specific_internal_energy_j_kg",
    "parent_phase_specific_enthalpy_j_kg", "parent_phase_mass_fractions",
    "parent_phase_root_identity",
]
GENERATION_ADVANCE = [
    "accepted-parent-publication", "rollback-invalidation",
    "final-whole-mixture-reflash-invalidation", "restored-parent-revalidation",
]
RUNTIME_EPOCH_ADVANCE = [
    "runtime-start", "restart-restore-before-parent-revalidation",
]
TOKEN_BINDS = [
    "committed_parent_q_fingerprint", "parent_phase_payload_fingerprint",
    "parent_generation", "runtime_epoch", "token_issue_sequence", "model_id",
    "phase_caloric_schema_id", "phase_caloric_options_fingerprint",
    "ordered_species", "energy_convention", "parent_phase_identity", "root_identity",
]
REPLAY_FIELDS = [
    "pressure_pa", "temperature_k", "compressibility_z", "density_kg_m3",
    "specific_internal_energy_j_kg", "specific_enthalpy_j_kg",
    "mass_fractions", "root_identity",
]
RELATIVE_RESIDUALS = [
    "pressure_pa", "temperature_k", "compressibility_z", "density_kg_m3",
    "specific_internal_energy_j_kg", "specific_enthalpy_j_kg",
]
REQUIRED_RESIDUALS = [
    "pressure_relative", "temperature_relative", "compressibility_relative",
    "density_relative", "specific_internal_energy_relative",
    "specific_enthalpy_relative", "mass_fraction_l1", "eos_root",
    "departure_caloric", "cv_minus_dedt_rho", "root_continuity",
]
RESIDUAL_ACCEPTANCE_MAX = {
    "pressure_relative": 1e-8,
    "temperature_relative": 1e-8,
    "compressibility_relative": 1e-8,
    "density_relative": 1e-8,
    "specific_internal_energy_relative": 1e-8,
    "specific_enthalpy_relative": 1e-8,
    "mass_fraction_l1": 1e-10,
    "eos_root": 1e-8,
    "departure_caloric": 1e-8,
    "cv_minus_dedt_rho": 1e-8,
    "root_continuity": 0.0,
}
FINAL_REFLASH_REQUIREMENTS = [
    "same-model-schema-options-species-energy-identity",
    "complete-rho-e-y-input-and-conservative-ledger-closure",
    "complete-equilibrium-state", "phase-and-acoustic-closure",
    "topology-admission", "final-parent-fingerprint-generation-and-token-invalidation",
]
NEGATIVE_CASES = [
    "liquid-payload-under-vapor-request",
    "vapor-supercritical-identity-mismatch",
    "seed-forced-alternate-root",
    "root-disappearance-or-spinodal-crossing",
    "near-critical-root-ambiguity",
    "composition-normalization-drift",
    "energy-reference-or-departure-caloric-drift",
    "parent-q-phase-payload-generation-or-runtime-epoch-mismatch",
    "live-token-reuse-across-rollback-restart-or-final-reflash",
    "restart-hint-direct-use-without-parent-revalidation-and-fresh-token",
    "initial-mid-rhs-or-final-failure-diagnostics-only-exact-rollback",
    "final-whole-mixture-reflash-enters-unsupported-three-phase-topology",
    "generic-solver-conversion-or-admission-attempt",
    "concurrent-or-reentrant-nondeterminism",
]
STOP_CONDITIONS = [
    "any approval flag becomes true",
    "PhaseCaloricResult becomes convertible to EquilibriumState",
    "quality acoustics sound-speed or global stability are added to phase-local result",
    "a failure result exposes state or an accepted live token",
    "parent-Q phase-payload generation runtime-epoch model schema options species energy or root binding drifts",
    "rollback or final reflash restores live-token validity from restored Q bytes",
    "restart hint is treated as an accepted live token or avoids fresh parent revalidation",
    "first recovery does not reproduce p T Z rho e h Y and root within declared residuals",
    "seed becomes branch authority or alternate-root fallback occurs",
    "vapor and supercritical branches auto-transition or liquid branch is admitted",
    "phase density basis silently becomes constant pressure or alpha authority changes",
    "phase-local result reaches Domain flux CFL transport checkpoint-equilibrium or generic solver admission",
    "final whole-mixture reflash is bypassed or unsupported topology is coerced",
    "the offline contract is represented as physical or production evidence",
]
ALLOWED_FILES = [
    "tests/reference/phase_caloric_feasibility_contract_v2.json",
    "tools/validate_phase_caloric_feasibility_contract_v2.py",
    "tests/test_phase_caloric_feasibility_contract_v2.py",
    "docs/phase_caloric_feasibility_contract_v2.md",
]


class ContractError(ValueError):
    """A structured fail-closed schema rejection."""


class Verifier:
    def __init__(self) -> None:
        self.checks = 0

    def exact(self, actual: Any, expected: Any, label: str) -> None:
        self.checks += 1
        if type(actual) is not type(expected) or actual != expected:
            raise ContractError(f"contract-mismatch:{label}")

    def obj(self, actual: Any, label: str, keys: set[str]) -> dict[str, Any]:
        self.checks += 1
        if not isinstance(actual, dict):
            raise ContractError(f"not-object:{label}")
        missing = keys - set(actual)
        unknown = set(actual) - keys
        if missing:
            raise ContractError(f"missing-keys:{label}:" + ",".join(sorted(missing)))
        if unknown:
            raise ContractError(f"unknown-keys:{label}:" + ",".join(sorted(unknown)))
        return actual

    def finite_positive(self, value: Any, label: str) -> None:
        self.checks += 1
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ContractError(f"invalid-number:{label}")
        if not math.isfinite(float(value)) or float(value) <= 0.0:
            raise ContractError(f"nonpositive-or-nonfinite:{label}")


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
        document = json.loads(
            raw.decode("utf-8"), object_pairs_hook=_reject_duplicate_pairs,
            parse_constant=_reject_nonfinite_json,
        )
    except ContractError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ContractError("contract-json-invalid") from exc
    if not isinstance(document, dict):
        raise ContractError("contract-not-object")
    return document, raw


def _verify(document: dict[str, Any], raw: bytes) -> int:
    v = Verifier()
    top = v.obj(document, "top", {
        "schema", "contract_version", "status", "scope", "approvals",
        "phase_caloric_schema", "parent_binding_schema", "root_token_schema",
        "first_recovery_and_residuals", "phase_local_validity", "chemistry_substep",
        "conservative_writeback_and_publication", "generic_solver_exclusion",
        "negative_cases_required", "evidence_boundary", "stop_conditions",
    })
    v.exact(top["schema"], SCHEMA, "schema")
    v.exact(top["contract_version"], 2, "contract-version")
    v.exact(top["status"], "blocked-offline-phase-caloric-feasibility-contract", "status")
    v.exact(top["scope"], "phase-local-vapor-or-supercritical-caloric-feasibility-no-runtime-implementation-no-two-phase-chemistry-success", "scope")
    approvals = v.obj(top["approvals"], "approvals", set(APPROVALS))
    for key, expected in APPROVALS.items():
        v.exact(approvals[key], expected, f"approval:{key}")
    if hashlib.sha256(raw).hexdigest() != EXPECTED_CONTRACT_SHA256:
        raise ContractError("contract-sha256-mismatch")

    phase = v.obj(top["phase_caloric_schema"], "phase-caloric-schema", {
        "id", "result_type", "result_is_equilibrium_state",
        "conversion_to_equilibrium_state", "forbidden_result_fields",
        "required_success_result_fields", "failure_result_policy",
    })
    v.exact(phase["id"], "rfam.phase-caloric/v2", "phase-schema-id")
    v.exact(phase["result_type"], "PhaseCaloricResult", "result-type")
    v.exact(phase["result_is_equilibrium_state"], False, "result-is-equilibrium-state")
    v.exact(phase["conversion_to_equilibrium_state"], "forbidden", "result-conversion")
    v.exact(phase["forbidden_result_fields"], FORBIDDEN_RESULT_FIELDS, "forbidden-result-fields")
    v.exact(phase["required_success_result_fields"], SUCCESS_RESULT_FIELDS, "success-result-fields")
    failure = v.obj(phase["failure_result_policy"], "failure-result-policy", {
        "result_kind", "forbidden_on_failure", "required_diagnostics",
    })
    v.exact(failure["result_kind"], "diagnostics-only-no-state-no-token", "failure-kind")
    v.exact(failure["forbidden_on_failure"], FAILURE_FORBIDDEN_FIELDS, "failure-forbidden-fields")
    v.exact(failure["required_diagnostics"], FAILURE_DIAGNOSTICS, "failure-diagnostics")

    binding = v.obj(top["parent_binding_schema"], "parent-binding", {
        "required_parent_fields", "generation_identity", "runtime_epoch",
        "required_phase_identities", "automatic_phase_transition",
        "automatic_vapor_to_supercritical", "automatic_vapor_to_liquid",
        "multi_root_fallback",
    })
    v.exact(binding["required_parent_fields"], PARENT_FIELDS, "parent-fields")
    v.exact(binding["required_phase_identities"], ["VAPOR", "SUPERCRITICAL"], "phase-identities")
    for key in ("automatic_phase_transition", "automatic_vapor_to_supercritical",
                "automatic_vapor_to_liquid", "multi_root_fallback"):
        v.exact(binding[key], "forbidden", f"branch-policy:{key}")
    generation = v.obj(binding["generation_identity"], "generation-identity", {
        "type", "advance_on", "restore_policy",
    })
    v.exact(generation["type"], "uint64-monotonic-accepted-generation", "generation-type")
    v.exact(generation["advance_on"], GENERATION_ADVANCE, "generation-advance")
    v.exact(generation["restore_policy"], "restored-q-bytes-must-not-restore-live-generation-or-token-validity", "generation-restore")
    epoch = v.obj(binding["runtime_epoch"], "runtime-epoch", {
        "type", "advance_on", "restore_policy",
    })
    v.exact(epoch["type"], "uint64-monotonic-runtime-epoch", "epoch-type")
    v.exact(epoch["advance_on"], RUNTIME_EPOCH_ADVANCE, "epoch-advance")
    v.exact(epoch["restore_policy"], "persisted-values-are-hints-not-live-epoch-authority", "epoch-restore")

    token = v.obj(top["root_token_schema"], "root-token", {
        "seed_role", "accepted_live_token_type", "token_issue_sequence_type",
        "token_binds", "previous_accepted_root_token_required_after_first",
        "root_continuity_required", "token_revalidation_against_parent_required",
        "invalidation_policy", "rollback_policy", "restart_hint",
    })
    v.exact(token["seed_role"], "numerical-hint-only-never-branch-authority", "seed-role")
    v.exact(token["accepted_live_token_type"], "opaque-parent-bound-runtime-epoch-root-continuity-token", "live-token-type")
    v.exact(token["token_issue_sequence_type"], "uint64-monotonic-per-runtime-epoch", "token-sequence-type")
    v.exact(token["token_binds"], TOKEN_BINDS, "token-binds")
    for key in ("previous_accepted_root_token_required_after_first",
                "root_continuity_required", "token_revalidation_against_parent_required"):
        v.exact(token[key], True, f"token-policy:{key}")
    v.exact(token["invalidation_policy"], "rollback-restart-final-reflash-or-parent-binding-drift-invalidates-live-token", "token-invalidation")
    v.exact(token["rollback_policy"], "advance-generation-and-invalidate-live-token-even-when-q-bytes-are-restored", "token-rollback")
    hint = v.obj(token["restart_hint"], "restart-hint", {
        "persistence", "never_accepted_as_live_token", "required_revalidation",
        "direct_reuse_across_restart",
    })
    v.exact(hint["persistence"], "opaque-hint-only", "hint-persistence")
    v.exact(hint["never_accepted_as_live_token"], True, "hint-live-token")
    v.exact(hint["required_revalidation"], [
        "complete-parent-replay", "fresh-runtime-epoch", "fresh-parent-generation",
        "fresh-live-token-issuance",
    ], "hint-revalidation")
    v.exact(hint["direct_reuse_across_restart"], "fail-closed", "hint-direct-reuse")

    recovery = v.obj(top["first_recovery_and_residuals"], "recovery-residuals", {
        "must_reproduce_parent_payload", "replay_required_fields",
        "relative_residual_max", "mass_fraction_l1_residual_max",
        "root_identity_must_match", "residual_schema",
    })
    v.exact(recovery["must_reproduce_parent_payload"], True, "parent-replay-required")
    v.exact(recovery["replay_required_fields"], REPLAY_FIELDS, "replay-fields")
    relative = v.obj(recovery["relative_residual_max"], "relative-residual-max", set(RELATIVE_RESIDUALS))
    for field in RELATIVE_RESIDUALS:
        v.finite_positive(relative[field], f"relative-residual:{field}")
    v.finite_positive(recovery["mass_fraction_l1_residual_max"], "mass-fraction-l1-limit")
    v.exact(recovery["root_identity_must_match"], True, "replay-root-identity")
    residual = v.obj(recovery["residual_schema"], "residual-schema", {
        "all_finite_required", "units", "required_residuals", "acceptance_max",
        "acceptance",
    })
    v.exact(residual["all_finite_required"], True, "residual-finite")
    v.exact(residual["units"], "dimensionless-normalized-residuals-except-mass-fraction-l1", "residual-units")
    v.exact(residual["required_residuals"], REQUIRED_RESIDUALS, "required-residuals")
    v.exact(residual["acceptance_max"], RESIDUAL_ACCEPTANCE_MAX, "residual-acceptance-max")
    v.exact(residual["acceptance"], "each-residual-must-meet-its-declared-finite-bound-before-token-issuance", "residual-acceptance")

    validity = v.obj(top["phase_local_validity"], "phase-local-validity", {
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
        v.exact(validity[key], True, f"validity:{key}")
    v.finite_positive(validity["composition_sum_absolute_tolerance"], "composition-tolerance")
    for key in ("root_loss", "spinodal_crossing", "near_critical_ambiguity"):
        v.exact(validity[key], "fail-closed", f"validity:{key}")
    v.exact(validity["metastable_vapor_statement"], "bounded-phase-local-transient-closure-only-not-global-equilibrium", "metastable-statement")

    substep = v.obj(top["chemistry_substep"], "chemistry-substep", {
        "basis", "constant_pressure_chemistry", "phase_density_relative_displacement_max",
        "specific_internal_energy_relative_displacement_max",
        "mass_fraction_l1_displacement_max", "envelope_exit",
        "initial_mid_rhs_final_recovery_required", "concurrent_reentrant_determinism_required",
    })
    v.exact(substep["basis"], "isochoric-fixed-phase-density-vapor-or-supercritical-chemistry", "substep-basis")
    v.exact(substep["constant_pressure_chemistry"], "forbidden", "constant-pressure")
    v.exact(substep["phase_density_relative_displacement_max"], 0.0, "fixed-density")
    v.finite_positive(substep["specific_internal_energy_relative_displacement_max"], "energy-displacement")
    v.finite_positive(substep["mass_fraction_l1_displacement_max"], "composition-displacement")
    v.exact(substep["envelope_exit"], "fail-closed-with-exact-rollback", "substep-exit")
    v.exact(substep["initial_mid_rhs_final_recovery_required"], True, "substep-stages")
    v.exact(substep["concurrent_reentrant_determinism_required"], True, "substep-determinism")

    publication = v.obj(top["conservative_writeback_and_publication"], "publication", {
        "mixture_volume_basis", "parent_vapor_volume_fraction_authority",
        "species_increment_formula", "energy_increment_formula",
        "vapor_volume_to_mixture_volume_rhoY_and_rhoE_required",
        "post_flash_radical_or_energy_fix", "only_publisher",
        "final_reflash_requirements", "three_phase_or_unsupported_topology",
    })
    v.exact(publication["mixture_volume_basis"], "committed-parent-mixture-volume", "mixture-basis")
    v.exact(publication["parent_vapor_volume_fraction_authority"], "committed-parent-phase-payload-only-not-recomputed-during-substep", "alpha-authority")
    v.exact(publication["species_increment_formula"], "delta-rhoY-kg-m3=alpha-parent-vapor-volume-fraction*rho-phase-kg-m3*delta-Y-phase", "species-formula")
    v.exact(publication["energy_increment_formula"], "delta-rhoE-j-m3=alpha-parent-vapor-volume-fraction*rho-phase-kg-m3*delta-e-phase-j-kg", "energy-formula")
    v.exact(publication["vapor_volume_to_mixture_volume_rhoY_and_rhoE_required"], True, "writeback-required")
    v.exact(publication["post_flash_radical_or_energy_fix"], "forbidden", "post-flash-fix")
    v.exact(publication["only_publisher"], "final-whole-mixture-EquilibriumThermoModel-recoverRhoE", "only-publisher")
    v.exact(publication["final_reflash_requirements"], FINAL_REFLASH_REQUIREMENTS, "final-reflash")
    v.exact(publication["three_phase_or_unsupported_topology"], "fail-closed", "topology-failure")

    exclusion = v.obj(top["generic_solver_exclusion"], "generic-solver-exclusion", {
        "phase_caloric_result_domain_publication",
        "phase_caloric_result_flux_cfl_transport_boundary",
        "phase_caloric_result_generic_solver_admission",
        "phase_caloric_result_checkpoint_as_equilibrium_state",
    })
    for key in exclusion:
        v.exact(exclusion[key], "forbidden", f"solver-exclusion:{key}")
    v.exact(top["negative_cases_required"], NEGATIVE_CASES, "negative-cases")

    evidence = v.obj(top["evidence_boundary"], "evidence-boundary", {
        "artifact_class", "physical_pre_validation", "runtime_implementation_present",
        "two_phase_chemistry_success_claimed", "source_files_allowed",
    })
    v.exact(evidence["artifact_class"], "OFFLINE_SCHEMA_AND_TEST_CONTRACT_ONLY", "artifact-class")
    for key in ("physical_pre_validation", "runtime_implementation_present",
                "two_phase_chemistry_success_claimed"):
        v.exact(evidence[key], False, f"evidence:{key}")
    v.exact(evidence["source_files_allowed"], ALLOWED_FILES, "allowed-files")
    v.exact(top["stop_conditions"], STOP_CONDITIONS, "stop-conditions")
    return v.checks


def validate(root: Path, contract_path: Path) -> int:
    if not root.resolve().is_dir():
        raise ContractError("root-not-directory")
    resolved = contract_path.resolve() if contract_path.is_absolute() else (root / contract_path).resolve()
    document, raw = _load(resolved)
    return _verify(document, raw)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--contract", type=Path,
                        default=Path("tests/reference/phase_caloric_feasibility_contract_v2.json"))
    args = parser.parse_args()
    try:
        checks = validate(args.root.resolve(), args.contract)
    except ContractError as exc:
        print(f"{FAIL_MARKER} {exc}", file=sys.stderr)
        return 2
    except Exception as exc:
        print(f"{FAIL_MARKER} unexpected:{type(exc).__name__}", file=sys.stderr)
        return 2
    print(f"{BLOCKED_MARKER} offline_schema=1 runtime_implementation=0 "
          f"two_phase_chemistry_enabled=0 phase_aware_chemistry_enabled=0 "
          f"production_ready=0 checks={checks}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
