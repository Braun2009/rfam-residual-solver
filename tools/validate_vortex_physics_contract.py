#!/usr/bin/env python3
"""Validate the fail-closed RFAM/Vortex WP0 physics contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Any


SCHEMA = "rfam.vortex-physics-validation-contract/v1"
SNAPSHOT_SCHEMA = "rfam.external-json-snapshot/v1"
EXPECTED_ROOT_KEYS = {
    "schema",
    "contract_version",
    "status",
    "production_ready",
    "mapping_freeze_approved",
    "fluent_execution_approved",
    "state_ownership",
    "energy_basis",
    "artifacts",
    "external_origins",
    "target_envelope",
    "species",
    "r2",
    "kij",
    "hem",
    "condensed_carbon",
    "comparison",
}
EXPECTED_ARTIFACT_PATHS = {
    "hotstart_config_review": "docs/v20_hotstart_config_review.md",
    "injection_snapshot": "tests/reference/vortex_injection_3pintle_estimate_v1.json",
    "cantera_generator": "scripts/cantera_reference.py",
    "cantera_reference": "tests/reference/cantera_gri30_propane_oxygen_v1.txt",
    "cantera_mechanism": "tests/reference/cantera_gri30_v3_0.yaml",
    "two_step_source": "src/combustion/FiniteRate/TwoStepCO.C",
    "two_step_header": "src/combustion/FiniteRate/TwoStepCO.H",
    "arrhenius_rate_header": "src/combustion/FiniteRate/ArrheniusRate.H",
    "dynamic_chemistry_source": "src/combustion/FiniteRate/DynamicChemistry.H",
    "two_step_yaml": "tests/twostep_co.yaml",
    "species_registry": "src/thermo/SpeciesRegistry.H",
    "binary_interaction_source": "src/thermo/MixingRules/BinaryInteractionDB.C",
    "binary_interaction_header": "src/thermo/MixingRules/BinaryInteractionDB.H",
    "r2_instantaneous_generator": "scripts/generate_r2_instantaneous_reference.py",
    "r2_instantaneous_reference": "tests/reference/r2_instantaneous_reference_v1.json",
    "r2_cantera_mechanism": "tests/reference/rfam_twostep_cantera_v1.yaml",
    "r2_equivalence_validator": "tools/validate_r2_instantaneous_equivalence.py",
    "r2_probe_source": "test_r2_instantaneous_compare.cpp",
    "r2_provenance_gate": "tests/r2_instantaneous_provenance_gate.cmake",
    "r2_regeneration_gate": "tests/r2_reference_regeneration_gate.cmake",
    "r2_mass_action_thermo": "src/combustion/FiniteRate/MassActionThermo.H",
    "r2_nasa7_thermo": "src/thermo/NASA7Thermo.H",
    "r2_mechanism_reader": "src/combustion/FiniteRate/MechanismReader.H",
    "r2_thermo_provenance_header": "src/thermo/ThermoDataProvenance.H",
    "r2_thermo_provenance_source": "src/thermo/ThermoDataProvenance.C",
}
EXPECTED_ARTIFACTS = set(EXPECTED_ARTIFACT_PATHS)
EXPECTED_EXTERNAL_SOURCE_PATH = "F:/Vortex/config/injection_3pintle_estimate.json"
EXPECTED_ANCHORS = {
    "pressure_pa": [101325.0, 4000000.0],
    "temperature_k": [130.0, 450.0, 700.0, 2200.0, 2800.0],
    "mass_flow_kg_s": [0.0449333, 0.35, 0.9333],
}
EXPECTED_ENVELOPE_REQUIREMENTS = [
    "pressure bounds and uncertainty",
    "temperature bounds and uncertainty",
    "composition and equivalence-ratio bounds",
    "dilution and residence-time bounds",
    "single-phase, VLE, near-critical and solid-risk partition",
]
EXPECTED_R2_OBSERVABLES = [
    "forward and reverse rates at prescribed T/rho/Y",
    "reaction orders and rate units",
    "formation-energy qdot conversion",
    "ignition delay",
    "adiabatic temperature",
    "major-species histories",
    "heat-release history",
]
EXPECTED_KIJ_ZERO_FIELDS = [
    "source",
    "applicability envelope",
    "uncertainty",
    "nonlinear sensitivity evidence",
]
EXPECTED_KIJ_OUTPUTS = [
    "phase envelope",
    "vapor mass fraction",
    "phase compositions",
    "density",
    "internal energy",
    "vapor-side O2 and C3H8 concentrations",
]
EXPECTED_SOOT_PHYSICS = [
    "nucleation",
    "surface growth",
    "oxidation",
    "wall deposition",
    "dilute zero-slip applicability",
]
EXPECTED_EXCLUSIONS = [
    "quality",
    "latent_energy",
    "rhoSoot",
    "number_density",
    "wall_carbon",
]


class ContractError(ValueError):
    pass


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise ContractError(f"artifact-open-failed:{path}:{exc}") from exc
    return digest.hexdigest()


def _require_object(value: Any, context: str, keys: set[str]) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ContractError(f"not-object:{context}")
    actual = set(value)
    missing = keys - actual
    unknown = actual - keys
    if missing:
        raise ContractError(f"missing-keys:{context}:" + ",".join(sorted(missing)))
    if unknown:
        raise ContractError(f"unknown-keys:{context}:" + ",".join(sorted(unknown)))
    return value


def _require_string(mapping: dict[str, Any], key: str, context: str) -> str:
    value = mapping[key]
    if not isinstance(value, str) or not value.strip():
        raise ContractError(f"invalid-string:{context}.{key}")
    return value


def _require_bool(mapping: dict[str, Any], key: str, context: str) -> bool:
    value = mapping[key]
    if type(value) is not bool:
        raise ContractError(f"invalid-bool:{context}.{key}")
    return value


def _require_int(mapping: dict[str, Any], key: str, context: str) -> int:
    value = mapping[key]
    if type(value) is not int:
        raise ContractError(f"invalid-int:{context}.{key}")
    return value


def _require_finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"invalid-number:{label}")
    result = float(value)
    if not math.isfinite(result):
        raise ContractError(f"nonfinite-number:{label}")
    return result


def _require_exact_strings(value: Any, expected: list[str], context: str) -> None:
    if not isinstance(value, list):
        raise ContractError(f"not-list:{context}")
    for index, item in enumerate(value):
        if not isinstance(item, str) or not item:
            raise ContractError(f"invalid-string:{context}[{index}]")
    if value != expected:
        raise ContractError(f"list-contract-mismatch:{context}")


def _valid_sha256(value: Any, context: str, allow_zero: bool = False) -> str:
    if not isinstance(value, str):
        raise ContractError(f"invalid-sha256:{context}")
    digest = value.lower()
    if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
        raise ContractError(f"invalid-sha256:{context}")
    if not allow_zero and digest == "0" * 64:
        raise ContractError(f"zero-sha256-forbidden:{context}")
    return digest


def _resolve_under(root: Path, relative: str, label: str) -> Path:
    if not isinstance(relative, str) or not relative.strip():
        raise ContractError(f"invalid-path:{label}")
    candidate = (root / relative).resolve()
    root_resolved = root.resolve()
    if candidate != root_resolved and root_resolved not in candidate.parents:
        raise ContractError(f"path-escapes-root:{label}")
    return candidate


def _external_path(value: str) -> Path:
    windows = re.fullmatch(r"([A-Za-z]):[\\/](.*)", value)
    if windows and os.name != "nt":
        tail = windows.group(2).replace("\\", "/")
        return Path("/mnt") / windows.group(1).lower() / tail
    return Path(value)


def _validate_artifacts(root: Path, value: Any) -> dict[str, Path]:
    if not isinstance(value, list):
        raise ContractError("not-list:artifacts")
    paths: dict[str, Path] = {}
    relative_paths: set[str] = set()
    for index, entry_value in enumerate(value):
        entry = _require_object(
            entry_value, f"artifacts[{index}]", {"id", "path", "sha256"}
        )
        artifact_id = _require_string(entry, "id", f"artifacts[{index}]")
        if artifact_id in paths:
            raise ContractError(f"duplicate-artifact:{artifact_id}")
        relative = _require_string(entry, "path", f"artifact[{artifact_id}]")
        actual_path = _resolve_under(root, relative, f"artifact[{artifact_id}]")
        canonical = EXPECTED_ARTIFACT_PATHS.get(artifact_id)
        if canonical is not None and relative != canonical:
            raise ContractError(
                "artifact-path-contract-mismatch:"
                f"{artifact_id}:expected={canonical}:actual={relative}"
            )
        if relative in relative_paths:
            raise ContractError(f"duplicate-artifact-path:{relative}")
        relative_paths.add(relative)
        expected = _valid_sha256(entry["sha256"], f"artifact[{artifact_id}]")
        if not actual_path.is_file():
            raise ContractError(f"artifact-missing:{artifact_id}:{relative}")
        actual = _sha256(actual_path)
        if actual != expected:
            raise ContractError(
                f"artifact-sha256-mismatch:{artifact_id}:expected={expected}:actual={actual}"
            )
        paths[artifact_id] = actual_path
    actual_ids = set(paths)
    if actual_ids != EXPECTED_ARTIFACTS:
        missing = EXPECTED_ARTIFACTS - actual_ids
        unknown = actual_ids - EXPECTED_ARTIFACTS
        detail = []
        if missing:
            detail.append("missing=" + ",".join(sorted(missing)))
        if unknown:
            detail.append("unknown=" + ",".join(sorted(unknown)))
        raise ContractError("artifact-set-mismatch:" + ":".join(detail))
    return paths


def _validate_origins(value: Any, artifacts: dict[str, Path]) -> None:
    if not isinstance(value, list) or len(value) != 1:
        raise ContractError("external-origin-count-mismatch")
    origin = _require_object(
        value[0],
        "external_origins[0]",
        {"id", "source_path", "source_sha256", "snapshot_artifact", "snapshot_schema"},
    )
    if _require_string(origin, "id", "external_origins[0]") != "vortex_injection_source":
        raise ContractError("external-origin-id-mismatch")
    source_path = _require_string(origin, "source_path", "external_origins[0]")
    if source_path != EXPECTED_EXTERNAL_SOURCE_PATH:
        raise ContractError(
            "external-origin-source-path-contract-mismatch:"
            f"expected={EXPECTED_EXTERNAL_SOURCE_PATH}:actual={source_path}"
        )
    source_sha = _valid_sha256(origin["source_sha256"], "external_origins[0].source_sha256")
    snapshot_id = _require_string(origin, "snapshot_artifact", "external_origins[0]")
    if snapshot_id not in artifacts:
        raise ContractError(f"external-origin-snapshot-artifact-missing:{snapshot_id}")
    if _require_string(origin, "snapshot_schema", "external_origins[0]") != SNAPSHOT_SCHEMA:
        raise ContractError("external-origin-snapshot-schema-mismatch")
    try:
        snapshot = json.loads(artifacts[snapshot_id].read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(f"external-snapshot-invalid:{exc}") from exc
    snapshot = _require_object(
        snapshot, "external_snapshot", {"schema", "source_path", "source_sha256", "payload"}
    )
    if snapshot["schema"] != SNAPSHOT_SCHEMA:
        raise ContractError("external-snapshot-schema-mismatch")
    if snapshot["source_path"] != source_path:
        raise ContractError("external-snapshot-source-path-mismatch")
    if snapshot["source_sha256"] != source_sha:
        raise ContractError("external-snapshot-source-sha256-mismatch")
    if not isinstance(snapshot["payload"], dict):
        raise ContractError("external-snapshot-payload-not-object")

    live_source = _external_path(source_path)
    if live_source.is_file():
        if _sha256(live_source) != source_sha:
            raise ContractError("external-source-sha256-mismatch")
        try:
            live_payload = json.loads(live_source.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ContractError(f"external-source-json-invalid:{exc}") from exc
        if snapshot["payload"] != live_payload:
            raise ContractError("external-snapshot-payload-mismatch")


def _validate_contract(root: Path, contract_path: Path) -> list[str]:
    try:
        document = json.loads(contract_path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ContractError(f"contract-open-failed:{exc}") from exc
    except json.JSONDecodeError as exc:
        raise ContractError(f"contract-json-invalid:{exc.msg}") from exc
    document = _require_object(document, "root", EXPECTED_ROOT_KEYS)
    if document["schema"] != SCHEMA:
        raise ContractError("schema-mismatch")
    if _require_int(document, "contract_version", "root") != 1:
        raise ContractError("contract-version-mismatch")
    if _require_string(document, "status", "root") != "blocked":
        raise ContractError("status-must-remain-blocked")
    for key in ("production_ready", "mapping_freeze_approved", "fluent_execution_approved"):
        if _require_bool(document, key, "root"):
            raise ContractError(f"{key}-must-be-false")

    ownership = _require_object(
        document["state_ownership"],
        "state_ownership",
        {"a413_udrgm", "rfam_hem", "quality_authoritative_uds", "condensed_carbon_uds"},
    )
    expected_owners = {
        "a413_udrgm": "VALIDATED_SINGLE_PHASE_T_P_Y",
        "rfam_hem": "RFAM_RHO_E_Y_FLASH",
        "quality_authoritative_uds": False,
        "condensed_carbon_uds": "A408_DILUTE_ZERO_SLIP_LEDGER_ONLY",
    }
    if ownership != expected_owners:
        raise ContractError("state-ownership-contract-mismatch")

    energy = _require_object(
        document["energy_basis"],
        "energy_basis",
        {"flow_energy", "formation_reference", "temperature_reference", "formation_energy_in_flow"},
    )
    expected_energy = {
        "flow_energy": "sensible_internal_energy",
        "formation_reference": "species_standard_formation_enthalpy_at_298.15_K",
        "temperature_reference": "298.15 K",
        "formation_energy_in_flow": False,
    }
    if energy != expected_energy:
        raise ContractError("energy-basis-contract-mismatch")

    artifacts = _validate_artifacts(root, document["artifacts"])
    _validate_origins(document["external_origins"], artifacts)

    envelope = _require_object(
        document["target_envelope"],
        "target_envelope",
        {"status", "known_anchors", "required_before_fit"},
    )
    if envelope["status"] != "incomplete":
        raise ContractError("target-envelope-must-be-incomplete")
    anchors = _require_object(
        envelope["known_anchors"], "target_envelope.known_anchors", set(EXPECTED_ANCHORS)
    )
    for key, expected in EXPECTED_ANCHORS.items():
        values = anchors[key]
        if not isinstance(values, list) or not values:
            raise ContractError(f"target-envelope-anchor-list:{key}")
        parsed = [_require_finite(value, f"target_envelope.{key}[{index}]")
                  for index, value in enumerate(values)]
        if any(value <= 0.0 for value in parsed):
            raise ContractError(f"target-envelope-anchor-nonpositive:{key}")
        if parsed != sorted(set(parsed)):
            raise ContractError(f"target-envelope-anchor-order-or-duplicate:{key}")
        if parsed != expected:
            raise ContractError(f"target-envelope-anchor-contract-mismatch:{key}")
    _require_exact_strings(
        envelope["required_before_fit"],
        EXPECTED_ENVELOPE_REQUIREMENTS,
        "target_envelope.required_before_fit",
    )

    species = _require_object(document["species"], "species", {"registry_order"})
    _require_exact_strings(
        species["registry_order"],
        ["O2", "N2", "H2", "CH4", "C2H6", "C3H8", "H2O", "CO", "CO2"],
        "species.registry_order",
    )

    r2 = _require_object(
        document["r2"],
        "r2",
        {"status", "instantaneous_rate_equivalence", "physical_holdout_validation",
         "detailed_mechanism", "experimental_dataset", "required_observables"},
    )
    if r2["status"] != "not-validated":
        raise ContractError("r2-status-must-remain-not-validated")
    if r2["instantaneous_rate_equivalence"] != "candidate-pass-awaiting-independent-review":
        raise ContractError("r2-instantaneous-status-mismatch")
    if r2["physical_holdout_validation"] != "not-run":
        raise ContractError("r2-physical-holdout-must-remain-not-run")
    if r2["detailed_mechanism"] != "unselected":
        raise ContractError("r2-detailed-mechanism-must-remain-unselected")
    if r2["experimental_dataset"] != "unselected":
        raise ContractError("r2-experimental-dataset-must-remain-unselected")
    _require_exact_strings(r2["required_observables"], EXPECTED_R2_OBSERVABLES, "r2.required_observables")

    kij = _require_object(
        document["kij"], "kij", {"status", "missing_pair_policy", "explicit_zero", "sensitivity_outputs"}
    )
    if kij["status"] != "incomplete":
        raise ContractError("kij-status-must-remain-incomplete")
    if kij["missing_pair_policy"] != "fail-closed":
        raise ContractError("kij-missing-policy-mismatch")
    explicit_zero = _require_object(
        kij["explicit_zero"], "kij.explicit_zero", {"requires_evidence", "required_fields"}
    )
    if _require_bool(explicit_zero, "requires_evidence", "kij.explicit_zero") is not True:
        raise ContractError("validated-explicit-zero-evidence-required")
    _require_exact_strings(
        explicit_zero["required_fields"], EXPECTED_KIJ_ZERO_FIELDS, "kij.explicit_zero.required_fields"
    )
    _require_exact_strings(kij["sensitivity_outputs"], EXPECTED_KIJ_OUTPUTS, "kij.sensitivity_outputs")

    hem = _require_object(
        document["hem"],
        "hem",
        {"status", "production_capability", "box_a", "box_b", "box_c", "box_dependencies"},
    )
    if hem["status"] != "foundation-only":
        raise ContractError("hem-status-must-remain-foundation-only")
    if _require_bool(hem, "production_capability", "hem"):
        raise ContractError("hem-production-capability-must-be-false")
    for key in ("box_a", "box_b", "box_c"):
        if hem[key] != "not-run":
            raise ContractError(f"hem-{key}-must-remain-not-run")
    dependencies = _require_object(
        hem["box_dependencies"], "hem.box_dependencies", {"box_a", "box_b", "box_c"}
    )
    expected_dependencies = {
        "box_a": "pure O2 reference; no kij or R2 dependency",
        "box_b": "validated active binary interactions",
        "box_c": "R2 implementation and holdout gates plus all active binary interactions",
    }
    if dependencies != expected_dependencies:
        raise ContractError("hem-box-dependency-contract-mismatch")

    carbon = _require_object(
        document["condensed_carbon"],
        "condensed_carbon",
        {"status", "quantitative_soot_comparison", "physics_required"},
    )
    if carbon["status"] != "ledger-only":
        raise ContractError("condensed-carbon-status-mismatch")
    if _require_bool(carbon, "quantitative_soot_comparison", "condensed_carbon"):
        raise ContractError("quantitative-soot-comparison-must-be-false")
    _require_exact_strings(
        carbon["physics_required"], EXPECTED_SOOT_PHYSICS, "condensed_carbon.physics_required"
    )

    comparison = _require_object(
        document["comparison"], "comparison", {"quantitative_single_phase_only", "excluded_quantities"}
    )
    if _require_bool(comparison, "quantitative_single_phase_only", "comparison") is not True:
        raise ContractError("comparison-must-remain-single-phase-only")
    _require_exact_strings(
        comparison["excluded_quantities"], EXPECTED_EXCLUSIONS, "comparison.excluded_quantities"
    )

    return [
        "detailed mechanism and experimental holdout not selected",
        "target envelope bounds remain incomplete",
        "R2 physical holdout and kij validation have not run",
        "RFAM HEM boxes and phase-aware capabilities remain disabled",
    ]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    args = parser.parse_args(argv)
    root = args.root.resolve()
    contract = args.contract
    if not contract.is_absolute():
        contract = root / contract
    try:
        blockers = _validate_contract(root, contract.resolve())
    except ContractError as error:
        print(f"VORTEX_PHYSICS_CONTRACT_FAIL {error}", file=sys.stderr)
        return 2
    except Exception as error:  # fail closed without leaking an unstructured traceback
        print(
            "VORTEX_PHYSICS_CONTRACT_FAIL "
            f"internal-error:{type(error).__name__}:{error}",
            file=sys.stderr,
        )
        return 2
    print("VORTEX_PHYSICS_CONTRACT_BLOCKED")
    print("status=blocked production_ready=0 mapping_freeze_approved=0 fluent_execution_approved=0")
    for blocker in blockers:
        print(f"blocked_reason={blocker}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
