#!/usr/bin/env python3
"""Validate the frozen same-mechanism R2 trajectory comparison contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path, PurePosixPath
from typing import Any


class ContractError(RuntimeError):
    pass


EXPECTED_ROOT_KEYS = {
    "schema",
    "contract_version",
    "status",
    "scope",
    "approvals",
    "artifacts",
    "species_order",
    "cases",
    "time_grid",
    "thermo_envelope",
    "energy_reference",
    "solver_profiles",
    "acceptance_tolerances",
    "expected_counts",
    "output_contract",
    "stop_conditions",
}

EXPECTED_APPROVALS = {
    "detailed_mechanism_selected": False,
    "experimental_data_selected": False,
    "r2_fitting_approved": False,
    "physical_validation_approved": False,
    "production_ready": False,
    "reacting_hot_start_approved": False,
}

EXPECTED_ARTIFACTS = [
    ("upstream_0d_contract", "tests/reference/r2_0d_calibration_contract_v1.json", "bbf01269ac1af96814446f462696ed62b5b57e80a777afc5a9d0e7fb25f14611"),
    ("initial_state_reference", "tests/reference/r2_0d_cantera_reference_v1.json", "11d1135a0d0beae0138a3b1dbdde9cae6ccdcee092c103d537311e72f497981d"),
    ("instantaneous_reference", "tests/reference/r2_instantaneous_reference_v1.json", "8bd3f2c05d70bb59af7668aa37c35bab3726c79d5ea8f3b650ba065980c1840a"),
    ("cantera_same_mechanism", "tests/reference/rfam_twostep_cantera_v1.yaml", "a8e88d5fd0114575afa258720f3b8c3dd6de558ebe00cc3ab82dec32cdb2f171"),
    ("rfam_runner_source", "tools/r2_0d_rfam_runner.cpp", "2f6acfe72015b4de91fe1ade806ea8f6d30bf8d9049fb7d89dc936e546eeabbb"),
    ("rfam_two_step_header", "src/combustion/FiniteRate/TwoStepCO.H", "9636dc4b2280af37d0231880a31672557556d5c9ee955bd4c46d373cad46d8c5"),
    ("rfam_two_step_source", "src/combustion/FiniteRate/TwoStepCO.C", "013d7af8f055a0ee834a095fa2946abd45214cf769ca5ccd43c3d83d16fd6465"),
    ("rfam_mass_action_thermo", "src/combustion/FiniteRate/MassActionThermo.H", "7cc012ea2b0bf9e6cee067b40a7590e3e5f3a27819b04fecca3edb8f2fe94d7b"),
    ("rfam_nasa7_thermo", "src/thermo/NASA7Thermo.H", "40ca86c9b52ee7332299ae8401680cc8ac84ad0dd035aab1183d4df9b5e5fafc"),
    ("rfam_species_registry", "src/thermo/SpeciesRegistry.H", "22b3c0e96d5bff8b6ffa6b69ecd7568a4f87422f5188af85de31c50136b27ba7"),
]

SPECIES = ["C3H8", "O2", "CO", "CO2", "H2O", "N2"]
CASE_IDS = [
    "train_phi0p8_diluted",
    "train_phi1p0",
    "train_phi1p2",
    "holdout_phi0p6_diluted",
    "holdout_phi1p1",
    "holdout_phi1p4",
]


def strict_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ContractError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def load_json(path: Path, label: str) -> Any:
    try:
        return json.loads(
            path.read_text(encoding="utf-8"), object_pairs_hook=strict_object
        )
    except ContractError:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContractError(f"json-load-failed:{label}:{error}") from error


def exact_object(value: Any, keys: set[str], label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ContractError(f"invalid-object:{label}")
    actual = set(value)
    if actual != keys:
        missing = ",".join(sorted(keys - actual))
        unknown = ",".join(sorted(actual - keys))
        raise ContractError(f"object-shape-mismatch:{label}:{missing}:{unknown}")
    return value


def exact(value: Any, expected: Any, label: str) -> None:
    if value != expected:
        raise ContractError(f"value-mismatch:{label}")


def finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"invalid-number:{label}")
    result = float(value)
    if not math.isfinite(result):
        raise ContractError(f"nonfinite-number:{label}")
    return result


def relative_path(root: Path, value: Any, label: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ContractError(f"invalid-path:{label}")
    path = PurePosixPath(value)
    if path.is_absolute() or "." in path.parts or ".." in path.parts:
        raise ContractError(f"path-escapes-root:{label}")
    candidate = (root / Path(*path.parts)).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as error:
        raise ContractError(f"path-escapes-root:{label}") from error
    return candidate


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise ContractError(f"artifact-open-failed:{path.name}:{error}") from error


def verify_artifacts(root: Path, artifacts: Any) -> dict[str, Path]:
    if not isinstance(artifacts, list) or len(artifacts) != len(EXPECTED_ARTIFACTS):
        raise ContractError("artifact-count-mismatch")
    paths: dict[str, Path] = {}
    for index, (expected_id, expected_path, expected_hash) in enumerate(
        EXPECTED_ARTIFACTS
    ):
        item = exact_object(
            artifacts[index], {"id", "path", "sha256"}, f"artifacts[{index}]"
        )
        exact(item["id"], expected_id, f"artifacts[{index}].id")
        exact(item["path"], expected_path, f"artifacts[{index}].path")
        exact(item["sha256"], expected_hash, f"artifacts[{index}].sha256")
        path = relative_path(root, item["path"], f"artifacts[{index}].path")
        if not path.is_file():
            raise ContractError(f"artifact-missing:{expected_id}")
        exact(sha256(path), expected_hash, f"artifact-bytes:{expected_id}")
        paths[expected_id] = path
    return paths


def expected_cases(reference: dict[str, Any]) -> list[dict[str, Any]]:
    if reference.get("major_species") != SPECIES:
        raise ContractError("reference-species-order-mismatch")
    cases = reference.get("cases")
    if not isinstance(cases, list) or [case.get("id") for case in cases] != CASE_IDS:
        raise ContractError("reference-case-order-mismatch")
    result = []
    for case in cases:
        initial = case.get("initial_major_species_mass_fractions")
        if not isinstance(initial, dict) or set(initial) != set(SPECIES):
            raise ContractError(f"reference-initial-composition-invalid:{case.get('id')}")
        result.append(
            {
                "id": case["id"],
                "partition": case["partition"],
                "temperature_k": case["temperature_k"],
                "pressure_pa": case["pressure_pa"],
                "mass_fractions": [initial[name] for name in SPECIES],
            }
        )
    return result


def validate(root: Path, contract_path: Path) -> int:
    contract = load_json(contract_path, "contract")
    exact_object(contract, EXPECTED_ROOT_KEYS, "contract")
    exact(contract["schema"], "rfam.r2-same-mechanism-trajectory-contract/v1", "schema")
    if type(contract["contract_version"]) is not int:
        raise ContractError("invalid-integer:contract_version")
    exact(contract["contract_version"], 1, "contract_version")
    exact(contract["status"], "blocked-implementation-equivalence-only", "status")
    exact(
        contract["scope"],
        "same-two-step-mechanism-ideal-gas-constant-pressure-no-fit-no-physical-validation",
        "scope",
    )
    exact(contract["approvals"], EXPECTED_APPROVALS, "approvals")

    artifacts = verify_artifacts(root, contract["artifacts"])
    reference = load_json(artifacts["initial_state_reference"], "initial-state-reference")
    exact(contract["species_order"], SPECIES, "species_order")
    exact(contract["cases"], expected_cases(reference), "cases")
    for index, case in enumerate(contract["cases"]):
        finite(case["temperature_k"], f"cases[{index}].temperature_k")
        finite(case["pressure_pa"], f"cases[{index}].pressure_pa")
        values = case["mass_fractions"]
        if any(finite(value, f"cases[{index}].Y") < 0.0 for value in values):
            raise ContractError(f"negative-composition:cases[{index}]")
        if abs(sum(values) - 1.0) > 2.0e-12:
            raise ContractError(f"composition-sum-mismatch:cases[{index}]")

    exact(
        contract["time_grid"],
        {
            "kind": "zero-plus-base10-logspace",
            "minimum_positive_time_s": 1.0e-14,
            "maximum_time_s": 1.0e-8,
            "log_points_including_endpoints": 121,
            "total_points_per_case": 122,
        },
        "time_grid",
    )
    exact(
        contract["thermo_envelope"],
        {
            "shared_minimum_temperature_k": 300.0,
            "comparison_maximum_temperature_k": 3400.0,
            "diagnostic_maximum_temperature_k": 3500.0,
            "eligibility_rule": "max-rfam-cantera-temperature-at-sample-le-3400K",
            "guard_rule": "3400K-lt-max-rfam-cantera-temperature-le-3500K-diagnostic-only",
            "outside_rule": "max-rfam-cantera-temperature-gt-3500K-thermo-extrapolation-excluded",
            "retain_all_samples_in_audit": True,
        },
        "thermo_envelope",
    )
    exact(
        contract["energy_reference"],
        {
            "reference_temperature_k": 298.15,
            "rfam_reference_policy": "nasa7-in-range-else-boundary-cp-linear-tail",
            "registry_formation_enthalpy_authority": "SpeciesRegistry",
            "reference_offset_definition": "h_ref_rfam_minus_hf_registry",
            "cantera_to_rfam_total_enthalpy": "h_total_rfam_equals_h_total_cantera_minus_reference_offset",
            "lower_linear_tail_species": ["C3H8", "N2"],
        },
        "energy_reference",
    )
    exact(
        contract["solver_profiles"],
        {
            "rfam": {
                "name": "strict-equivalence-v1",
                "relative_tolerance": 1.0e-12,
                "species_absolute_tolerance": 1.0e-16,
                "sensible_enthalpy_absolute_tolerance_j_per_kg": 1.0e-8,
                "maximum_internal_step_s": 1.0e-13,
            },
            "cantera": {
                "version": "3.1.0",
                "reactor": "IdealGasConstPressureReactor",
                "relative_tolerance": 1.0e-12,
                "absolute_tolerance": 1.0e-18,
                "maximum_time_step_s": 1.0e-13,
            },
        },
        "solver_profiles",
    )
    exact(
        contract["acceptance_tolerances"],
        {
            "trajectory_temperature_absolute_k": 2.0e-5,
            "trajectory_mass_fraction_absolute": 2.0e-9,
            "ignition_delay_absolute_s": 2.0e-15,
            "qdot_relative": 2.0e-8,
            "qdot_absolute_w_per_m3": 1.0e-3,
            "species_derivative_relative": 2.0e-8,
            "species_derivative_absolute_per_s": 1.0e-8,
            "temperature_derivative_relative": 2.0e-8,
            "temperature_derivative_absolute_k_per_s": 1.0e-6,
        },
        "acceptance_tolerances",
    )
    exact(
        contract["expected_counts"],
        {
            "cases": 6,
            "total_samples": 732,
            "eligible_samples": 646,
            "guard_samples": 0,
            "outside_samples": 86,
        },
        "expected_counts",
    )
    exact(
        contract["output_contract"],
        {
            "schema": "rfam.r2-same-mechanism-trajectory-evidence/v1",
            "required_case_order": CASE_IDS,
            "score_only_eligible_samples": True,
            "report_guard_and_outside_samples": True,
            "allow_sample_truncation": False,
            "fit_executed": False,
            "physical_validation_claimed": False,
            "production_ready": False,
        },
        "output_contract",
    )
    exact(
        contract["stop_conditions"],
        [
            "artifact path or sha256 mismatch",
            "case, species, grid, energy-reference or solver-profile drift",
            "any approval flag true",
            "any sample omitted from the audit",
            "guard or outside sample included in an acceptance metric",
            "any acceptance tolerance weakened",
            "attempt to treat implementation equivalence as physical R2 validation",
        ],
        "stop_conditions",
    )
    return len(EXPECTED_ARTIFACTS) + len(contract["cases"]) + 13


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    args = parser.parse_args()
    try:
        root = args.root.resolve()
        contract = args.contract.resolve()
        if not root.is_dir():
            raise ContractError("root-not-directory")
        checks = validate(root, contract)
    except ContractError as error:
        print(f"R2_SAME_MECHANISM_TRAJECTORY_CONTRACT_FAIL {error}", file=sys.stderr)
        return 2
    except Exception as error:
        print(
            "R2_SAME_MECHANISM_TRAJECTORY_CONTRACT_FAIL unexpected:"
            + type(error).__name__,
            file=sys.stderr,
        )
        return 2
    print(
        "R2_SAME_MECHANISM_TRAJECTORY_CONTRACT_BLOCKED "
        f"checks={checks} cases=6 samples=732 eligible=646 outside=86 "
        "fit_approved=0 physical_validation_approved=0 production_ready=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
