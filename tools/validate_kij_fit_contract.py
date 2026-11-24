#!/usr/bin/env python3
"""Validate the frozen O2/C3H8 PR kij pre-fit contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any


SCHEMA = "rfam.kij-fit-contract/v1"
BLOCKED_MARKER = "RFAM_KIJ_FIT_CONTRACT_BLOCKED"
FAIL_MARKER = "RFAM_KIJ_FIT_CONTRACT_FAIL"

EXPECTED_APPROVALS = {
    "data_rights_resolved": False,
    "raw_dataset_selected": False,
    "training_holdout_frozen": False,
    "fit_form_selected": False,
    "fit_executed": False,
    "kij_value_validated": False,
    "hem_enabled": False,
    "fluent_execution_approved": False,
    "hot_start_approved": False,
    "production_ready": False,
}

EXPECTED_UPSTREAM = [
    {
        "id": "wp3_sensitivity_contract",
        "path": "realFluidAllMach/tests/reference/kij_sensitivity_contract_v1.json",
        "sha256": "e1501157afdf1d2d51bb6fa0d03413b6da3f4d0cee4a0fbd9b8cfa9a39c99752",
    },
    {
        "id": "wp3_source_candidate_contract",
        "path": "realFluidAllMach/tests/reference/kij_o2_c3h8_source_candidate_v1.json",
        "sha256": "014eb6abd75f943a6a34c49d20b262fd61a3a16d3616babd8afd729a0c6bde5e",
    },
]

EXPECTED_PAIR_IDENTITY = {
    "components": [
        {
            "name": "O2",
            "formula": "O2",
            "inchi_key": "MYMOFIZGZYHOMD-UHFFFAOYSA-N",
        },
        {
            "name": "C3H8",
            "formula": "C3H8",
            "inchi_key": "ATUOYWHBWRKTHZ-UHFFFAOYSA-N",
        },
    ],
    "canonical_pair_key": "C3H8|O2",
    "canonical_order": "lexicographic-component-name",
    "eos": "Peng-Robinson",
    "mixing_rule": "van-der-Waals-one-fluid",
    "parameter": "binary-interaction-kij",
    "symmetry_requirement": (
        "kij(O2,C3H8) equals kij(C3H8,O2) exactly at every evaluated temperature"
    ),
    "missing_pair_policy": "fail-closed-when-nonlinear-sensitivity-is-material",
    "parameter_state": "not-fit-no-published-or-runtime-value",
}

EXPECTED_DATA_PROVENANCE = {
    "status": "no-authorized-fit-dataset",
    "source_candidate_role": "screening-only-not-training-or-holdout",
    "raw_artifact_repository_copy": False,
    "required_artifact_identity_fields": [
        "artifact_id",
        "provider",
        "citation_or_dataset_doi",
        "immutable_url_or_local_archive_uri",
        "sha256",
        "byte_count",
        "retrieved_on",
        "license_or_terms_url",
        "terms_sha256",
    ],
    "required_data_dictionary_fields": [
        "column_name",
        "physical_quantity",
        "unit",
        "basis",
        "phase_label",
        "uncertainty_kind",
        "uncertainty_confidence_or_coverage",
        "missing_value_reason",
    ],
    "no_unbound_reference_backend_values": True,
}

EXPECTED_RIGHTS = {
    "status": "unresolved-artifact-specific-numerical-use-and-redistribution",
    "fit_and_redistribution_permitted": False,
    "required_before_fit": [
        "artifact-specific terms allowing the intended numerical fitting use",
        (
            "artifact-specific terms allowing the intended repository redistribution "
            "or an explicit no-copy workflow"
        ),
        (
            "attribution and modification-notice requirements recorded in the "
            "selected-data manifest"
        ),
    ],
    "public_access_not_equivalent_to_permission": True,
}

EXPECTED_MEASUREMENT_SCHEMA = {
    "accepted_measurement_classes": ["binary-vle", "binary-pvt"],
    "composition_basis": "mole-fraction",
    "mass_density_unit": "kg/m3",
    "temperature_unit": "K",
    "pressure_unit": "Pa",
    "phase_label_enum": [
        "single-liquid",
        "single-vapor",
        "vapor-liquid-equilibrium",
        "vapor-liquid-liquid-equilibrium",
        "critical-or-near-critical",
        "solid-risk",
        "unknown",
    ],
    "required_common_fields": [
        "artifact_id",
        "dataset_id",
        "record_id",
        "split",
        "temperature_k",
        "pressure_pa",
        "measured_phase_label",
        "temperature_standard_uncertainty_k",
        "pressure_standard_uncertainty_pa",
    ],
    "required_vle_tie_line_fields": [
        "liquid_x_o2",
        "vapor_y_o2",
        "liquid_x_o2_standard_uncertainty",
        "vapor_y_o2_standard_uncertainty",
    ],
    "required_pvt_fields": [
        "overall_z_o2",
        "mixture_density_kg_m3",
        "overall_z_o2_standard_uncertainty",
        "mixture_density_standard_uncertainty_kg_m3",
    ],
    "phase_label_must_be_measured_or_source-qualified": True,
    "missing_uncertainty_policy": "forbid-empty-or-zero-substitution",
    "vle_fit_requires_both-liquid-and-vapor-composition": True,
}

EXPECTED_PARTITION = {
    "status": "unselected",
    "frozen_before_fit": True,
    "training_dataset_ids": [],
    "holdout_dataset_ids": [],
    "partition_unit": "independent-source-dataset-or-experimental-campaign",
    "disjointness_requirements": [
        "a record cannot appear in both training and holdout",
        "holdout data must not determine fit-form selection or parameter values",
        (
            "holdout source provenance must be independent from the primary training "
            "campaign unless a documented non-overlap justification is approved"
        ),
    ],
    "holdout_required_before_validation_claim": True,
}

EXPECTED_FIT_DEFINITION = {
    "status": "not-selected-not-executed",
    "candidate_functional_forms": [
        {"id": "constant", "expression": "kij(T)=c0", "parameter_count": 1},
        {
            "id": "reduced-temperature-linear",
            "expression": "kij(T)=c0+c1*(T_ref_K/T-1)",
            "parameter_count": 2,
        },
    ],
    "reference_temperature_definition": (
        "T_ref_K must be frozen with the selected form before fitting"
    ),
    "selected_form_id": None,
    "objective_and_weights": "unselected",
    "parameter_values": [],
    "fit_performed": False,
    "symmetry_enforcement": (
        "one canonical C3H8|O2 parameter record is queried symmetrically in both "
        "component orders"
    ),
    "historical_kij_0_055_policy": (
        "not-a-start-value-fit-target-fallback-or-validated-value"
    ),
}

EXPECTED_APPLICABILITY_ENVELOPE = {
    "status": "unselected",
    "required_before_fit": [
        "temperature range in K",
        "pressure range in Pa",
        "overall and phase-composition ranges on mole basis",
        "allowed phase labels",
        "near-critical exclusion or validated treatment",
        "solid-risk exclusion",
        "extrapolation policy",
    ],
    "outside_envelope_policy": "fail-closed-no-default-kij-or-single-root-fallback",
    "validated_envelope": None,
}

EXPECTED_UNCERTAINTY = {
    "status": "not-estimated",
    "required_before_validation_claim": [
        "measurement uncertainty model and coverage conversion",
        "fit-parameter covariance or bounded uncertainty set",
        "model-form uncertainty treatment",
        "uncertainty propagation through VLE/PVT observables",
        "holdout uncertainty-aware residuals",
    ],
    "zero_uncertainty_without_source_evidence": "forbidden",
    "parameter_uncertainty": None,
}

EXPECTED_NONLINEAR_SENSITIVITY = {
    "status": "not-run",
    "required_method": (
        "nonlinear sweep of the selected kij form and its declared uncertainty over "
        "the target T-p-z envelope"
    ),
    "required_observables": [
        "phase_region",
        "mass_quality",
        "mixture_density_kg_m3",
        "mixture_specific_internal_energy_j_kg",
        "liquid_x_o2",
        "vapor_y_o2",
        "vapor_o2_concentration_mol_m3",
        "vapor_c3h8_concentration_mol_m3",
    ],
    "required_outputs": [
        "state-by-state absolute-and-relative-spans",
        "phase-topology-change-records",
        "material-and-energy-closure-residuals",
        "sensitivity-classification-and-threshold-rationale",
    ],
    "material_missing_pair_action": "fail-closed-before-hem-or-production-use",
    "sensitivity_complete": False,
}

EXPECTED_ACCEPTANCE = {
    "status": "metrics-and-tolerances-unselected",
    "required_before_fit": [
        (
            "training residual definitions for pressure, liquid composition, vapor "
            "composition, density, and any selected PVT observable"
        ),
        "holdout residual definitions using the same physical bases and units",
        "uncertainty-aware weights and aggregation rule",
        "predeclared pass thresholds and phase-topology failure policy",
        "independent holdout report",
    ],
    "fit_success_is_not_validation": True,
    "validation_requires_independent_holdout": True,
}

EXPECTED_LIMITATIONS = [
    "this contract selects no external VLE or PVT dataset and stores no raw experimental data",
    "the current ThermoML candidate remains source-screening-only and is not a fit input",
    "Cantera is not a PR kij source, VLE/PVT fitting authority, or kij validation backend",
    "no value including historical kij=0.055 is fitted, published, recommended, or available to runtime",
    "no HEM, R2 physics, soot physics, Fluent, zero-advance, hot-start, or production capability is enabled",
]

EXPECTED_ROOT_KEYS = (
    "schema",
    "contract_version",
    "status",
    "scope",
    "approvals",
    "upstream",
    "pair_identity",
    "data_provenance",
    "rights",
    "measurement_schema",
    "partition",
    "fit_definition",
    "applicability_envelope",
    "uncertainty",
    "nonlinear_sensitivity",
    "acceptance",
    "limitations",
)


class KijFitContractError(ValueError):
    """A fail-closed pre-fit contract violation."""


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise KijFitContractError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def _decode_contract(path: Path) -> dict[str, Any]:
    try:
        payload = path.read_bytes()
    except OSError as exc:
        raise KijFitContractError(f"contract-read-failed:{exc}") from exc
    if b"\r" in payload or not payload.endswith(b"\n"):
        raise KijFitContractError("contract-not-canonical-lf")
    try:
        text = payload.decode("ascii")
        document = json.loads(text, object_pairs_hook=_reject_duplicate_pairs)
    except KijFitContractError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise KijFitContractError(f"contract-json-invalid:{exc}") from exc
    if not isinstance(document, dict):
        raise KijFitContractError("contract-root-not-object")
    canonical = json.dumps(document, indent=2, ensure_ascii=True) + "\n"
    if text != canonical:
        raise KijFitContractError("contract-not-canonical-json")
    return document


def _strict_equal(actual: Any, expected: Any, label: str) -> None:
    if type(actual) is not type(expected):
        raise KijFitContractError(f"contract-type-mismatch:{label}")
    if isinstance(expected, dict):
        actual_keys = list(actual)
        expected_keys = list(expected)
        if actual_keys != expected_keys:
            missing = [key for key in expected_keys if key not in actual]
            unknown = [key for key in actual_keys if key not in expected]
            if missing:
                raise KijFitContractError(
                    f"contract-missing-keys:{label}:" + ",".join(missing)
                )
            if unknown:
                raise KijFitContractError(
                    f"contract-unknown-keys:{label}:" + ",".join(unknown)
                )
            raise KijFitContractError(f"contract-key-order-mismatch:{label}")
        for key in expected_keys:
            _strict_equal(actual[key], expected[key], f"{label}.{key}")
        return
    if isinstance(expected, list):
        if len(actual) != len(expected):
            raise KijFitContractError(f"contract-list-length-mismatch:{label}")
        for index, expected_item in enumerate(expected):
            _strict_equal(actual[index], expected_item, f"{label}[{index}]")
        return
    if actual != expected:
        raise KijFitContractError(f"contract-value-mismatch:{label}")


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _workspace_path(root: Path, relative: str, label: str) -> Path:
    if not isinstance(relative, str) or not relative or "\\" in relative:
        raise KijFitContractError(f"invalid-path:{label}")
    if relative.startswith("/") or re.match(r"^[A-Za-z]:", relative):
        raise KijFitContractError(f"path-not-relative:{label}")
    resolved_root = root.resolve()
    candidate = (resolved_root / relative).resolve()
    if candidate != resolved_root and resolved_root not in candidate.parents:
        raise KijFitContractError(f"path-escapes-root:{label}")
    return candidate


def validate_contract(root: Path, contract_path: Path) -> dict[str, Any]:
    document = _decode_contract(contract_path)
    if list(document) != list(EXPECTED_ROOT_KEYS):
        missing = [key for key in EXPECTED_ROOT_KEYS if key not in document]
        unknown = [key for key in document if key not in EXPECTED_ROOT_KEYS]
        if missing:
            raise KijFitContractError("contract-missing-keys:root:" + ",".join(missing))
        if unknown:
            raise KijFitContractError("contract-unknown-keys:root:" + ",".join(unknown))
        raise KijFitContractError("contract-key-order-mismatch:root")

    _strict_equal(document["schema"], SCHEMA, "root.schema")
    _strict_equal(document["contract_version"], 1, "root.contract_version")
    _strict_equal(document["status"], "blocked-pre-fit-contract", "root.status")
    _strict_equal(
        document["scope"],
        "offline-o2-c3h8-pr-kij-fit-provenance-and-acceptance-schema-no-fit",
        "root.scope",
    )
    _strict_equal(document["approvals"], EXPECTED_APPROVALS, "root.approvals")
    _strict_equal(document["upstream"], EXPECTED_UPSTREAM, "root.upstream")
    _strict_equal(
        document["pair_identity"], EXPECTED_PAIR_IDENTITY, "root.pair_identity"
    )
    _strict_equal(
        document["data_provenance"], EXPECTED_DATA_PROVENANCE, "root.data_provenance"
    )
    _strict_equal(document["rights"], EXPECTED_RIGHTS, "root.rights")
    _strict_equal(
        document["measurement_schema"],
        EXPECTED_MEASUREMENT_SCHEMA,
        "root.measurement_schema",
    )
    _strict_equal(document["partition"], EXPECTED_PARTITION, "root.partition")
    _strict_equal(
        document["fit_definition"], EXPECTED_FIT_DEFINITION, "root.fit_definition"
    )
    _strict_equal(
        document["applicability_envelope"],
        EXPECTED_APPLICABILITY_ENVELOPE,
        "root.applicability_envelope",
    )
    _strict_equal(document["uncertainty"], EXPECTED_UNCERTAINTY, "root.uncertainty")
    _strict_equal(
        document["nonlinear_sensitivity"],
        EXPECTED_NONLINEAR_SENSITIVITY,
        "root.nonlinear_sensitivity",
    )
    _strict_equal(document["acceptance"], EXPECTED_ACCEPTANCE, "root.acceptance")
    _strict_equal(document["limitations"], EXPECTED_LIMITATIONS, "root.limitations")

    for item in EXPECTED_UPSTREAM:
        path = _workspace_path(root, item["path"], f"upstream:{item['id']}")
        try:
            payload = path.read_bytes()
        except OSError as exc:
            raise KijFitContractError(
                f"upstream-read-failed:{item['id']}:{exc}"
            ) from exc
        if _sha256(payload) != item["sha256"]:
            raise KijFitContractError(f"upstream-sha256-mismatch:{item['id']}")
    return document


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[2]
    default_contract = (
        Path(__file__).resolve().parents[1]
        / "tests/reference/kij_fit_contract_v1.json"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument("--contract", type=Path, default=default_contract)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        root = args.root.resolve()
        contract = args.contract
        if not contract.is_absolute():
            contract = (Path.cwd() / contract).resolve()
        validate_contract(root, contract)
    except Exception as exc:  # Unexpected input must remain structured and fail closed.
        reason = (
            str(exc)
            if isinstance(exc, KijFitContractError)
            else f"unexpected:{type(exc).__name__}:{exc}"
        )
        reason = reason.replace("\r", " ").replace("\n", " ")
        print(f"{FAIL_MARKER} reason={reason}")
        return 2

    print(
        f"{BLOCKED_MARKER} upstream=2 raw_dataset_selected=0 "
        "fit_executed=0 validated_kij=0 production_ready=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
