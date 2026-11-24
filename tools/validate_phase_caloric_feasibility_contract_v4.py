#!/usr/bin/env python3
"""Fail-closed validator for the blocked phase-caloric v4 design contract."""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
from pathlib import Path
import re
import struct
import sys
import types
from typing import Any, Dict, Mapping, Optional, Tuple


SCHEMA = "rfam.phase-caloric-feasibility-contract/v4"
BLOCKED_MARKER = "PHASE_CALORIC_FEASIBILITY_CONTRACT_V4_BLOCKED"
FAIL_MARKER = "PHASE_CALORIC_FEASIBILITY_CONTRACT_V4_FAIL"

CONTRACT_PATH = "tests/reference/phase_caloric_feasibility_contract_v4.json"
ORACLE_PATH = "tools/phase_caloric_lifecycle_oracle_v4.py"
VALIDATOR_PATH = "tools/validate_phase_caloric_feasibility_contract_v4.py"
TEST_PATH = "tests/test_phase_caloric_feasibility_contract_v4.py"
DOC_PATH = "docs/phase_caloric_feasibility_contract_v4.md"

EXPECTED_CONTRACT_SHA256 = (
    "035f4e3d003b01bf8f61c1e84fd0715fa2ebfecdf88ae954da7dc102e42ee33b"
)
EXPECTED_ORACLE_SHA256 = (
    "a0630407174dc18a53817fb5eb6e94dd8b68ef1e832eb16c8cf95ee97a434cc3"
)
EXPECTED_DOC_SHA256 = (
    "892ff0692cfc7650e4ae8676f3a746cb7241a5a1059dcada7a9e6bbdf1596a71"
)

SECTION_SHA256 = {
    "approvals": "9b2d31871c92f6faf338b3e54d20143b8f4b6934c3ae02b5e03f4770c5423c0f",
    "artifact_binding": "78417c05a7bc668dcfff9d85a30e502e2eab07d3dc9833fca53c0e667b5a1067",
    "canonical_serialization": "02ac89a7565e1ac58cdabf0a67100e6294b84362b9f6901fe6b8aee16e73e2ff",
    "typed_object_schemas": "e3219985d3d2396bf4d401bf928b581905785d268f23656ab3b3ba37ddec16e9",
    "lifecycle_ownership": "99a4f143d61b0fc0d6c9481366e11be75adbbcef4ed85e78d331bae112afe4e2",
    "event_order": "c1838072c050d85ed94da7abcedcd0bc84b4a31d308438a0e025400939e20fbe",
    "residual_formulas": "b36403ac046b833aa6c054ea1333c7443bfb655cb744bda1284ca1208d93da54",
    "applicability_envelope": "be63bff59c0ccfdd074f0037d02f738c1664d06226ef6fdce74a2a76b14881f4",
    "publication_boundary": "41c0a02dc289c4366634e1eec70980760b8ce9ee1f936578784b4ace315bdc06",
    "negative_cases_required": "e284c8e787d44669b13925f84bed70dd5b2a38b8386e4bc67ad737a71d4010e8",
    "evidence_boundary": "ff8f4cec26dff1e393a73c3eedb3da81b11bef92a21d96788be5ab2370025377",
    "stop_conditions": "394554b12eabeacd01f7c314c33edeab7098662e1382a91c7c9c5d18aa6cd763",
}

TOP_KEYS = {
    "schema",
    "contract_version",
    "status",
    "scope",
    "approvals",
    "artifact_binding",
    "canonical_serialization",
    "typed_object_schemas",
    "lifecycle_ownership",
    "event_order",
    "residual_formulas",
    "applicability_envelope",
    "publication_boundary",
    "negative_cases_required",
    "evidence_boundary",
    "stop_conditions",
}

APPROVAL_KEYS = {
    "phase_caloric_implementation_approved",
    "two_phase_chemistry_enabled",
    "phase_aware_chemistry_enabled",
    "equilibrium_result_admission_relaxed",
    "solver_consumer_admission",
    "checkpoint_runtime_integration",
    "restart_runtime_integration",
    "numerical_displacement_limits_validated",
    "physical_validation_approved",
    "fluent_execution_approved",
    "reacting_hot_start_approved",
    "production_ready",
}

SHA256_RE = re.compile(r"[0-9a-f]{64}")
OWNER_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.:/-]{0,127}")
F64BE_RE = re.compile(r"[0-9a-f]{16}")
UINT64_MAX = (1 << 64) - 1
FINGERPRINT_FIELDS = {
    "committed_parent_q_fingerprint": [
        "energy_convention",
        "ordered_species_count",
        "ordered_species_fingerprint",
        "parent_owner_id",
        "rho_c_kg_m3_f64be_ordered_scalars",
        "rho_e_j_m3_f64be",
        "rho_kg_m3_f64be",
        "rho_u_kg_m2_s_f64be_xyz",
        "rho_y_kg_m3_f64be_ordered_species",
        "rho_z_kg_m3_f64be_ordered_scalars",
        "schema",
    ],
    "parent_phase_payload_fingerprint": [
        "compressibility_z_f64be",
        "density_kg_m3_f64be",
        "energy_convention",
        "mass_fractions_f64be_ordered_species",
        "model_id",
        "ordered_species_count",
        "ordered_species_fingerprint",
        "parent_owner_id",
        "parent_phase_identity",
        "phase_caloric_options_fingerprint",
        "phase_caloric_schema_id",
        "pressure_pa_f64be",
        "root_identity",
        "schema",
        "specific_enthalpy_j_kg_f64be",
        "specific_internal_energy_j_kg_f64be",
        "temperature_k_f64be",
    ],
}


class ContractError(ValueError):
    """Structured fail-closed contract rejection."""


class Verifier:
    def __init__(self) -> None:
        self.checks = 0

    def check(self, condition: bool, label: str) -> None:
        self.checks += 1
        if not condition:
            raise ContractError(label)

    def exact(self, actual: Any, expected: Any, label: str) -> None:
        self.check(type(actual) is type(expected) and actual == expected, label)

    def obj(self, value: Any, keys: set[str], label: str) -> Dict[str, Any]:
        self.check(isinstance(value, dict), f"not-object:{label}")
        actual = set(value)
        self.check(not (keys - actual), f"missing-keys:{label}")
        self.check(not (actual - keys), f"unknown-keys:{label}")
        return value


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ContractError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def _reject_nonfinite_json(token: str) -> None:
    raise ContractError(f"nonfinite-json-number:{token}")


def _read(path: Path, label: str) -> bytes:
    try:
        return path.read_bytes()
    except OSError as exc:
        raise ContractError(f"artifact-read-failed:{label}") from exc


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _section_sha(value: Any) -> str:
    payload = json.dumps(
        value,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("ascii")
    return _sha256(payload)


def _canonical_fingerprint(domain: str, value: Mapping[str, Any]) -> str:
    payload = json.dumps(
        value,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("ascii")
    digest = hashlib.sha256()
    digest.update(domain.encode("ascii"))
    digest.update(b"\0")
    digest.update(payload)
    return digest.hexdigest()


def _decode(raw: bytes) -> Dict[str, Any]:
    try:
        value = json.loads(
            raw.decode("ascii"),
            object_pairs_hook=_reject_duplicate_pairs,
            parse_constant=_reject_nonfinite_json,
        )
    except ContractError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ContractError("contract-json-invalid") from exc
    if not isinstance(value, dict):
        raise ContractError("contract-not-object")
    return value


def _canonical_path(root: Path, relative: str, label: str) -> Path:
    if (
        not isinstance(relative, str)
        or not relative
        or "\\" in relative
        or relative.startswith("/")
        or re.match(r"^[A-Za-z]:", relative)
    ):
        raise ContractError(f"invalid-canonical-path:{label}")
    resolved_root = root.resolve()
    candidate = (resolved_root / relative).resolve()
    if candidate != resolved_root and resolved_root not in candidate.parents:
        raise ContractError(f"canonical-path-escapes-root:{label}")
    return candidate


def _valid_sha(value: Any, label: str) -> str:
    if (
        not isinstance(value, str)
        or SHA256_RE.fullmatch(value) is None
        or value == "0" * 64
    ):
        raise ContractError(f"invalid-sha256:{label}")
    return value


def _uint64(value: Any, label: str, *, positive: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ContractError(f"invalid-uint64:{label}")
    minimum = 1 if positive else 0
    if value < minimum or value > UINT64_MAX:
        raise ContractError(f"uint64-out-of-range:{label}")
    return value


def _finite_binary64(value: Any, label: str, *, positive: bool = False) -> float:
    if type(value) is not float or not math.isfinite(value):
        raise ContractError(f"invalid-binary64:{label}")
    if positive and value <= 0.0:
        raise ContractError(f"nonpositive-binary64:{label}")
    return value


def _f64be(value: Any, label: str) -> float:
    if not isinstance(value, str) or F64BE_RE.fullmatch(value) is None:
        raise ContractError(f"fingerprint-f64be:{label}")
    decoded = struct.unpack(">d", bytes.fromhex(value))[0]
    if not math.isfinite(decoded):
        raise ContractError(f"fingerprint-f64be-nonfinite:{label}")
    return decoded


def _nonblank_ascii(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ContractError(f"invalid-string:{label}")
    try:
        value.encode("ascii")
    except UnicodeEncodeError as exc:
        raise ContractError(f"non-ascii-string:{label}") from exc
    return value


def _load_verified_oracle(root: Path) -> types.ModuleType:
    path = _canonical_path(root, ORACLE_PATH, "lifecycle-oracle")
    payload = _read(path, "lifecycle-oracle")
    if _sha256(payload) != EXPECTED_ORACLE_SHA256:
        raise ContractError("artifact-sha256-mismatch:lifecycle-oracle")
    try:
        source = payload.decode("ascii")
    except UnicodeDecodeError as exc:
        raise ContractError("lifecycle-oracle-not-ascii") from exc
    name = "rfam_phase_caloric_lifecycle_oracle_v4_verified"
    module = types.ModuleType(name)
    module.__file__ = str(path)
    sys.modules[name] = module
    try:
        exec(compile(source, str(path), "exec"), module.__dict__)
    except Exception as exc:
        sys.modules.pop(name, None)
        raise ContractError(
            f"lifecycle-oracle-load-failed:{type(exc).__name__}"
        ) from exc
    return module


def _verify_meta_schema(document: Dict[str, Any], verifier: Verifier) -> None:
    typed = document["typed_object_schemas"]
    verifier.check(isinstance(typed, dict), "not-object:typed-object-schemas")
    verifier.exact(
        set(typed),
        {
            "common_object_policy",
            "request",
            "success",
            "failure",
            "live_root_token",
            "restart_hint",
            "residual_object",
            "residual_summary",
        },
        "typed-schema-names",
    )
    for name in ("request", "success"):
        schema = verifier.obj(
            typed[name], {"required_fields", "properties"}, f"schema:{name}"
        )
        verifier.check(
            isinstance(schema["required_fields"], list),
            f"schema-required-fields:{name}",
        )
        verifier.check(
            isinstance(schema["properties"], dict), f"schema-properties:{name}"
        )
        verifier.exact(
            schema["required_fields"],
            list(schema["properties"]),
            f"schema-required-property-order:{name}",
        )
    failure = verifier.obj(
        typed["failure"],
        {"required_fields", "properties", "forbidden_semantic_categories"},
        "schema:failure",
    )
    verifier.exact(
        failure["required_fields"],
        list(failure["properties"]),
        "schema-required-property-order:failure",
    )
    for name in ("live_root_token", "residual_summary"):
        schema = typed[name]
        verifier.check(isinstance(schema, dict), f"schema-object:{name}")
        verifier.exact(
            schema["required_fields"],
            list(schema["properties"]),
            f"schema-required-property-order:{name}",
        )
        verifier.exact(
            schema["additional_properties"], False, f"schema-closed:{name}"
        )
    restart_hint = verifier.obj(
        typed["restart_hint"],
        {
            "required_fields",
            "additional_properties",
            "accepted_as_live_token",
            "authority",
            "reuse_policy",
            "properties",
        },
        "schema:restart-hint",
    )
    verifier.exact(
        restart_hint["required_fields"],
        list(restart_hint["properties"]),
        "schema-required-property-order:restart_hint",
    )
    verifier.exact(
        restart_hint["additional_properties"],
        False,
        "schema-closed:restart_hint",
    )
    verifier.exact(
        restart_hint["accepted_as_live_token"],
        False,
        "restart-hint-live-token-admission",
    )
    verifier.exact(
        restart_hint["authority"],
        "non-authoritative-numerical-guidance-only",
        "restart-hint-authority",
    )
    verifier.exact(
        restart_hint["reuse_policy"],
        "reusable-after-each-exact-current-parent-replay-each-reuse-mints-fresh-sequence-one-live-token",
        "restart-hint-reuse-policy",
    )
    residual = verifier.obj(
        typed["residual_object"],
        {"required_fields", "additional_properties", "property_type"},
        "schema:residual-object",
    )
    verifier.exact(
        residual["required_fields"],
        document["residual_formulas"]["required_residuals"],
        "residual-schema-fields",
    )
    verifier.exact(
        residual["additional_properties"], False, "residual-schema-closed"
    )
    verifier.exact(
        residual["property_type"],
        "finite-nonnegative-binary64",
        "residual-schema-type",
    )
    verifier.check(
        "committed_parent_q_fingerprint"
        in typed["success"]["required_fields"],
        "success-missing-committed-parent-q-fingerprint",
    )
    verifier.check(
        "parent_phase_identity" in typed["success"]["required_fields"],
        "success-missing-parent-phase-identity",
    )
    for name in ("request", "success"):
        verifier.check(
            "ordered_species_count" in typed[name]["required_fields"],
            f"{name}-missing-ordered-species-count",
        )
        verifier.exact(
            typed[name]["properties"]["ordered_species_count"],
            {"type": "uint64-positive"},
            f"{name}-ordered-species-count-type",
        )
    verifier.check(
        "parent_q_fingerprint" not in typed["success"]["required_fields"],
        "ambiguous-parent-q-fingerprint-admitted",
    )


def _verify_document(
    document: Dict[str, Any],
    raw: bytes,
    root: Path,
    *,
    expected_contract_sha256: str,
) -> int:
    verifier = Verifier()
    top = verifier.obj(document, TOP_KEYS, "top")
    verifier.exact(top["schema"], SCHEMA, "contract-mismatch:schema")
    verifier.exact(top["contract_version"], 4, "contract-mismatch:version")
    verifier.exact(
        top["status"], "blocked-offline-design-only", "contract-mismatch:status"
    )
    verifier.exact(
        top["scope"],
        "typed-phase-caloric-identity-and-lifecycle-design-no-runtime-owner-no-thermodynamic-implementation",
        "contract-mismatch:scope",
    )
    if _sha256(raw) != _valid_sha(expected_contract_sha256, "contract-authority"):
        raise ContractError("contract-sha256-mismatch")

    approvals = verifier.obj(top["approvals"], APPROVAL_KEYS, "approvals")
    for name in sorted(APPROVAL_KEYS):
        verifier.exact(approvals[name], False, f"approval-not-false:{name}")

    for name, expected in SECTION_SHA256.items():
        verifier.check(name in top, f"missing-section:{name}")
        verifier.exact(
            _section_sha(top[name]), expected, f"section-sha256-mismatch:{name}"
        )

    binding = top["artifact_binding"]
    paths = binding["canonical_relative_paths"]
    verifier.exact(
        paths,
        {
            "contract": CONTRACT_PATH,
            "lifecycle_oracle": ORACLE_PATH,
            "validator": VALIDATOR_PATH,
            "focused_test": TEST_PATH,
            "documentation": DOC_PATH,
        },
        "artifact-canonical-paths",
    )
    verifier.exact(
        binding["lifecycle_oracle_sha256"],
        EXPECTED_ORACLE_SHA256,
        "artifact-contract-oracle-sha256",
    )
    verifier.exact(
        binding["documentation_sha256"],
        EXPECTED_DOC_SHA256,
        "artifact-contract-doc-sha256",
    )
    for name, relative in paths.items():
        artifact = _canonical_path(root, relative, f"artifact:{name}")
        verifier.check(artifact.is_file(), f"artifact-missing:{name}")
    verifier.exact(
        _sha256(_read(_canonical_path(root, ORACLE_PATH, "oracle"), "oracle")),
        EXPECTED_ORACLE_SHA256,
        "artifact-sha256-mismatch:lifecycle-oracle",
    )
    verifier.exact(
        _sha256(_read(_canonical_path(root, DOC_PATH, "docs"), "docs")),
        EXPECTED_DOC_SHA256,
        "artifact-sha256-mismatch:documentation",
    )

    _verify_meta_schema(top, verifier)
    canonical = top["canonical_serialization"]
    verifier.obj(
        canonical,
        {
            "algorithm",
            "encoding",
            "json_object_keys",
            "json_separators",
            "json_duplicate_keys",
            "json_nonfinite_numbers",
            "binary64_encoding",
            "integer_encoding",
            "string_encoding",
            "committed_parent_q_fingerprint",
            "parent_phase_payload_fingerprint",
        },
        "canonical-serialization",
    )
    for field, expected in (
        ("algorithm", "sha256-domain-separated-canonical-json-v1"),
        ("encoding", "ascii"),
        ("json_object_keys", "lexicographically-sorted"),
        ("json_separators", "comma-colon-no-whitespace"),
        ("json_duplicate_keys", "forbidden"),
        ("json_nonfinite_numbers", "forbidden"),
        ("integer_encoding", "unsigned-decimal-uint64-no-bool-alias"),
        ("string_encoding", "json-ascii-with-standard-escapes"),
    ):
        verifier.exact(
            canonical[field], expected, f"canonical-serialization-policy:{field}"
        )
    verifier.exact(
        canonical["binary64_encoding"],
        "f64be-16-lowercase-hex-preserve-sign-bit-finite-only",
        "canonical-binary64-policy",
    )
    verifier.exact(
        canonical["committed_parent_q_fingerprint"]["domain"],
        "rfam.phase-caloric.committed-parent-q/v4",
        "committed-parent-q-domain",
    )
    verifier.exact(
        canonical["parent_phase_payload_fingerprint"]["domain"],
        "rfam.phase-caloric.parent-phase-payload/v4",
        "parent-phase-payload-domain",
    )
    for name, fields in FINGERPRINT_FIELDS.items():
        spec = verifier.obj(
            canonical[name],
            {
                "domain",
                "required_fields",
                "field_set_policy",
                "input_object_insertion_order",
            },
            f"canonical-fingerprint-spec:{name}",
        )
        verifier.exact(
            spec["required_fields"],
            fields,
            f"canonical-fingerprint-required-fields:{name}",
        )
        verifier.exact(
            spec["required_fields"],
            sorted(spec["required_fields"]),
            f"canonical-fingerprint-field-list-order:{name}",
        )
        verifier.exact(
            spec["field_set_policy"],
            "exact-no-missing-or-extra",
            f"canonical-fingerprint-field-set-policy:{name}",
        )
        verifier.exact(
            spec["input_object_insertion_order"],
            "ignored-canonical-serialization-sorts-keys-lexicographically",
            f"canonical-fingerprint-input-order:{name}",
        )

    ownership = top["lifecycle_ownership"]
    verifier.exact(
        ownership["process_runtime_epoch"]["owner"],
        "process-runtime-owner",
        "runtime-epoch-owner",
    )
    verifier.exact(
        ownership["parent_generation"]["owner"],
        "committed-parent-state-owner-per-parent-owner-id",
        "parent-generation-owner",
    )
    verifier.exact(
        ownership["transaction_sequence"]["owner"],
        "one-phase-caloric-transaction-object",
        "transaction-sequence-owner",
    )
    verifier.exact(
        ownership["transaction_sequence"]["different_parents_have_independent_sequences"],
        True,
        "transaction-sequence-parent-independence",
    )

    limits = top["applicability_envelope"]
    verifier.exact(
        limits["numerical_limit_validation_status"],
        "UNVALIDATED_DESIGN_LIMITS_DO_NOT_AUTHORIZE_RUNTIME",
        "numerical-limit-validation-status",
    )
    verifier.exact(
        approvals["numerical_displacement_limits_validated"],
        False,
        "numerical-limits-approval",
    )
    verifier.exact(
        limits["runtime_emission_or_use"],
        "forbidden-until-separately-validated-and-approved",
        "numerical-limit-runtime-use",
    )

    evidence = top["evidence_boundary"]
    verifier.exact(
        evidence["source_files_allowed"],
        [CONTRACT_PATH, ORACLE_PATH, VALIDATOR_PATH, TEST_PATH, DOC_PATH],
        "evidence-source-files",
    )
    for field in (
        "runtime_phase_caloric_owner_present",
        "thermodynamic_recovery_present",
        "physical_validation_present",
        "two_phase_chemistry_success_claimed",
    ):
        verifier.exact(evidence[field], False, f"evidence-boundary:{field}")
    return verifier.checks


def _check_property_type(
    document: Dict[str, Any],
    type_spec: Dict[str, Any],
    value: Any,
    label: str,
    instance: Dict[str, Any],
) -> None:
    if not isinstance(type_spec, dict) or "type" not in type_spec:
        raise ContractError(f"invalid-property-type-spec:{label}")
    kind = type_spec["type"]
    if kind == "literal":
        if type(value) is not type(type_spec["value"]) or value != type_spec["value"]:
            raise ContractError(f"typed-value-mismatch:{label}")
    elif kind == "enum":
        if not isinstance(value, str) or value not in type_spec["values"]:
            raise ContractError(f"typed-enum-mismatch:{label}")
    elif kind == "nullable-enum":
        if value is not None and (
            not isinstance(value, str) or value not in type_spec["values"]
        ):
            raise ContractError(f"typed-nullable-enum-mismatch:{label}")
    elif kind == "parent-owner-id":
        if not isinstance(value, str) or OWNER_RE.fullmatch(value) is None:
            raise ContractError(f"typed-parent-owner-id:{label}")
    elif kind == "sha256-lower-hex-nonzero":
        _valid_sha(value, label)
    elif kind == "nullable-sha256-lower-hex-nonzero":
        if value is not None:
            _valid_sha(value, label)
    elif kind == "uint64-positive":
        _uint64(value, label, positive=True)
    elif kind == "nonblank-ascii-string":
        _nonblank_ascii(value, label)
    elif kind == "bool":
        if type(value) is not bool:
            raise ContractError(f"typed-bool:{label}")
    elif kind == "finite-binary64":
        _finite_binary64(value, label)
    elif kind == "finite-positive-binary64":
        _finite_binary64(value, label, positive=True)
    elif kind == "nullable-finite-positive-binary64":
        if value is not None:
            _finite_binary64(value, label, positive=True)
    elif kind == "nullable-finite-nonnegative-binary64":
        if value is not None and _finite_binary64(value, label) < 0.0:
            raise ContractError(f"typed-negative-binary64:{label}")
    elif kind == "mass-fraction-vector-ordered-species":
        if not isinstance(value, list) or not value:
            raise ContractError(f"typed-mass-fraction-vector:{label}")
        species_count = _uint64(
            instance.get("ordered_species_count"),
            f"{label}:ordered-species-count",
            positive=True,
        )
        if len(value) != species_count:
            raise ContractError(f"typed-ordered-species-count:{label}")
        total = 0.0
        for index, item in enumerate(value):
            fraction = _finite_binary64(item, f"{label}[{index}]")
            if fraction < 0.0 or fraction > 1.0:
                raise ContractError(f"typed-mass-fraction-range:{label}[{index}]")
            total += fraction
        if not math.isfinite(total) or abs(total - 1.0) > 1.0e-12:
            raise ContractError(f"typed-mass-fraction-sum:{label}")
    elif kind == "phase-compatible-root-identity":
        phase = instance.get("parent_phase_identity")
        expected = {
            "VAPOR": "PR_VAPOR_ROOT",
            "SUPERCRITICAL": "PR_UNIQUE_SUPERCRITICAL_ROOT",
        }.get(phase)
        if value != expected:
            raise ContractError(f"typed-phase-root-mismatch:{label}")
    elif kind == "live-root-token-v4":
        validate_typed_instance(document, "live_root_token", value)
    elif kind == "nullable-live-root-token-v4":
        if value is not None:
            validate_typed_instance(document, "live_root_token", value)
    elif kind == "residual-object-v4":
        validate_typed_instance(document, "residual_object", value)
    elif kind == "residual-summary-v4":
        validate_typed_instance(document, "residual_summary", value)
    else:
        raise ContractError(f"unknown-property-type:{label}:{kind}")


def _parent_binding_digest(instance: Dict[str, Any]) -> str:
    root_identity = {
        "VAPOR": "PR_VAPOR_ROOT",
        "SUPERCRITICAL": "PR_UNIQUE_SUPERCRITICAL_ROOT",
    }[instance["parent_phase_identity"]]
    return _canonical_fingerprint(
        "rfam.phase-caloric.parent-binding/v4",
        {
            "committed_parent_q_fingerprint": instance[
                "committed_parent_q_fingerprint"
            ],
            "energy_convention": instance["energy_convention"],
            "model_id": instance["model_id"],
            "ordered_species_count": instance["ordered_species_count"],
            "ordered_species_fingerprint": instance[
                "ordered_species_fingerprint"
            ],
            "parent_generation": instance["parent_generation"],
            "parent_owner_id": instance["parent_owner_id"],
            "parent_phase_identity": instance["parent_phase_identity"],
            "parent_phase_payload_fingerprint": instance[
                "parent_phase_payload_fingerprint"
            ],
            "phase_caloric_options_fingerprint": instance[
                "phase_caloric_options_fingerprint"
            ],
            "phase_caloric_schema_id": instance["phase_caloric_schema_id"],
            "process_runtime_epoch": instance["process_runtime_epoch"],
            "root_identity": root_identity,
        },
    )


def _validate_nested_live_token(
    instance: Dict[str, Any], token: Dict[str, Any], schema_name: str
) -> None:
    root_identity = {
        "VAPOR": "PR_VAPOR_ROOT",
        "SUPERCRITICAL": "PR_UNIQUE_SUPERCRITICAL_ROOT",
    }[instance["parent_phase_identity"]]
    expected = {
        "parent_owner_id": instance["parent_owner_id"],
        "transaction_id": instance["transaction_id"],
        "process_runtime_epoch": instance["process_runtime_epoch"],
        "parent_generation": instance["parent_generation"],
        "parent_binding_digest": _parent_binding_digest(instance),
        "root_identity": root_identity,
    }
    if schema_name == "success":
        expected["transaction_sequence"] = instance["transaction_sequence"]
    for field, value in expected.items():
        if type(token[field]) is not type(value) or token[field] != value:
            raise ContractError(f"typed-live-token-binding:{schema_name}:{field}")
    token_payload = {
        "parent_binding_digest": token["parent_binding_digest"],
        "parent_generation": token["parent_generation"],
        "parent_owner_id": token["parent_owner_id"],
        "process_runtime_epoch": token["process_runtime_epoch"],
        "root_identity": token["root_identity"],
        "transaction_id": token["transaction_id"],
        "transaction_sequence": token["transaction_sequence"],
    }
    if token["token_digest"] != _canonical_fingerprint(
        "rfam.phase-caloric.live-root-token/v4", token_payload
    ):
        raise ContractError(f"typed-live-token-digest:{schema_name}")


def _validate_restart_hint_digest(instance: Dict[str, Any]) -> None:
    payload = {
        "committed_parent_q_fingerprint": instance[
            "committed_parent_q_fingerprint"
        ],
        "parent_generation": instance["source_parent_generation"],
        "parent_owner_id": instance["parent_owner_id"],
        "parent_phase_identity": instance["parent_phase_identity"],
        "parent_phase_payload_fingerprint": instance[
            "parent_phase_payload_fingerprint"
        ],
        "root_identity": instance["root_identity"],
        "source_process_runtime_epoch": instance[
            "source_process_runtime_epoch"
        ],
        "source_token_digest": instance["source_token_digest"],
    }
    if instance["hint_digest"] != _canonical_fingerprint(
        "rfam.phase-caloric.restart-hint/v4", payload
    ):
        raise ContractError("typed-restart-hint-digest")


def validate_typed_instance(
    document: Dict[str, Any],
    schema_name: str,
    instance: Any,
) -> None:
    typed = document["typed_object_schemas"]
    if schema_name not in typed or not isinstance(instance, dict):
        raise ContractError(f"typed-instance-not-object:{schema_name}")
    schema = typed[schema_name]
    required = schema["required_fields"]
    if set(instance) != set(required):
        missing = sorted(set(required) - set(instance))
        unknown = sorted(set(instance) - set(required))
        if missing:
            raise ContractError(
                f"typed-instance-missing:{schema_name}:" + ",".join(missing)
            )
        raise ContractError(
            f"typed-instance-unknown:{schema_name}:" + ",".join(unknown)
        )
    if schema_name == "residual_object":
        bounds = document["residual_formulas"]["acceptance_max"]
        for name in required:
            value = _finite_binary64(instance[name], f"{schema_name}.{name}")
            if value < 0.0 or value > bounds[name]:
                raise ContractError(f"residual-out-of-range:{name}")
        return
    properties = schema["properties"]
    for name in required:
        _check_property_type(
            document,
            properties[name],
            instance[name],
            f"{schema_name}.{name}",
            instance,
        )
    if schema_name == "residual_summary":
        if instance["all_finite"]:
            if (
                instance["maximum_normalized_residual"] is None
                or instance["failed_residual"] is not None
            ):
                raise ContractError("residual-summary-consistency")
        elif instance["failed_residual"] is None:
            raise ContractError("residual-summary-failure-name-missing")
    elif schema_name == "request":
        token = instance["previous_accepted_root_token"]
        if token is not None:
            _validate_nested_live_token(instance, token, schema_name)
    elif schema_name == "success":
        _validate_nested_live_token(
            instance, instance["accepted_root_token"], schema_name
        )
    elif schema_name == "restart_hint":
        _validate_restart_hint_digest(instance)


def fingerprint_payload(
    document: Dict[str, Any],
    oracle: types.ModuleType,
    fingerprint_name: str,
    payload: Dict[str, Any],
) -> str:
    if fingerprint_name not in (
        "committed_parent_q_fingerprint",
        "parent_phase_payload_fingerprint",
    ):
        raise ContractError(f"unknown-fingerprint-name:{fingerprint_name}")
    spec = document["canonical_serialization"][fingerprint_name]
    required = spec["required_fields"]
    if not isinstance(payload, dict) or set(payload) != set(required):
        raise ContractError(f"fingerprint-field-coverage:{fingerprint_name}")
    species_count = _uint64(
        payload.get("ordered_species_count"),
        f"fingerprint-ordered-species-count:{fingerprint_name}",
        positive=True,
    )
    decoded: Dict[str, Any] = {}
    for name, value in payload.items():
        if name.endswith("_f64be"):
            decoded[name] = _f64be(value, f"{fingerprint_name}:{name}")
        elif "_f64be_" in name:
            if not isinstance(value, list):
                raise ContractError(f"fingerprint-f64be-array:{fingerprint_name}:{name}")
            if name == "rho_u_kg_m2_s_f64be_xyz" and len(value) != 3:
                raise ContractError("fingerprint-rhou-dimension")
            if name in (
                "rho_y_kg_m3_f64be_ordered_species",
                "mass_fractions_f64be_ordered_species",
            ):
                if not value:
                    raise ContractError(
                        f"fingerprint-f64be-array:{fingerprint_name}:{name}"
                    )
                if len(value) != species_count:
                    raise ContractError(
                        f"fingerprint-ordered-species-count:{fingerprint_name}:{name}"
                    )
            decoded[name] = [
                _f64be(item, f"{fingerprint_name}:{name}[{index}]")
                for index, item in enumerate(value)
            ]
    if payload.get("energy_convention") != (
        "SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT"
    ):
        raise ContractError(f"fingerprint-energy-convention:{fingerprint_name}")
    if not isinstance(payload.get("parent_owner_id"), str) or OWNER_RE.fullmatch(
        payload["parent_owner_id"]
    ) is None:
        raise ContractError(f"fingerprint-parent-owner:{fingerprint_name}")
    _valid_sha(
        payload.get("ordered_species_fingerprint"),
        f"fingerprint-ordered-species:{fingerprint_name}",
    )
    if fingerprint_name == "committed_parent_q_fingerprint":
        if payload.get("schema") != "rfam.committed-parent-q/v4":
            raise ContractError("fingerprint-q-schema")
        if decoded["rho_kg_m3_f64be"] <= 0.0:
            raise ContractError("fingerprint-q-density-nonpositive")
        if any(
            value < 0.0
            for value in decoded["rho_y_kg_m3_f64be_ordered_species"]
        ):
            raise ContractError("fingerprint-rhoy-negative")
    elif fingerprint_name == "parent_phase_payload_fingerprint":
        if payload.get("schema") != "rfam.parent-phase-payload/v4":
            raise ContractError("fingerprint-phase-schema")
        phase = payload["parent_phase_identity"]
        root = payload["root_identity"]
        expected = {
            "VAPOR": "PR_VAPOR_ROOT",
            "SUPERCRITICAL": "PR_UNIQUE_SUPERCRITICAL_ROOT",
        }.get(phase)
        if root != expected:
            raise ContractError("fingerprint-phase-root-mismatch")
        for name in (
            "pressure_pa_f64be",
            "temperature_k_f64be",
            "compressibility_z_f64be",
            "density_kg_m3_f64be",
        ):
            if decoded[name] <= 0.0:
                raise ContractError(f"fingerprint-phase-nonpositive:{name}")
        fractions = decoded["mass_fractions_f64be_ordered_species"]
        if any(value < 0.0 or value > 1.0 for value in fractions):
            raise ContractError("fingerprint-phase-mass-fraction-range")
        if abs(sum(fractions) - 1.0) > 1.0e-12:
            raise ContractError("fingerprint-phase-mass-fraction-sum")
        _nonblank_ascii(payload.get("model_id"), "fingerprint-model-id")
        if payload.get("phase_caloric_schema_id") != "rfam.phase-caloric/v4":
            raise ContractError("fingerprint-phase-caloric-schema")
        _valid_sha(
            payload.get("phase_caloric_options_fingerprint"),
            "fingerprint-phase-caloric-options",
        )
    else:
        raise ContractError(f"unknown-fingerprint-name:{fingerprint_name}")
    try:
        return oracle.canonical_fingerprint(spec["domain"], payload)
    except Exception as exc:
        raise ContractError(f"fingerprint-oracle-rejected:{fingerprint_name}") from exc


def validate(
    root: Path,
    contract_path: Optional[Path] = None,
    *,
    expected_contract_sha256: str = EXPECTED_CONTRACT_SHA256,
    enforce_runtime_path: bool = True,
    import_oracle: bool = True,
) -> Tuple[Dict[str, Any], int, Optional[types.ModuleType]]:
    resolved_root = root.resolve()
    if not resolved_root.is_dir():
        raise ContractError("root-not-directory")
    canonical_validator = _canonical_path(
        resolved_root, VALIDATOR_PATH, "validator-runtime"
    )
    if enforce_runtime_path and Path(__file__).resolve() != canonical_validator:
        raise ContractError("validator-runtime-path-mismatch")
    canonical_contract = _canonical_path(resolved_root, CONTRACT_PATH, "contract")
    selected = canonical_contract if contract_path is None else contract_path.resolve()
    if selected != canonical_contract:
        raise ContractError("contract-canonical-path-mismatch")
    raw = _read(selected, "contract")
    document = _decode(raw)
    checks = _verify_document(
        document,
        raw,
        resolved_root,
        expected_contract_sha256=expected_contract_sha256,
    )
    oracle = _load_verified_oracle(resolved_root) if import_oracle else None
    return document, checks, oracle


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--contract", type=Path)
    args = parser.parse_args()
    try:
        document, checks, _ = validate(args.root, args.contract)
    except ContractError as exc:
        print(f"{FAIL_MARKER} {exc}", file=sys.stderr)
        return 2
    except Exception as exc:
        print(f"{FAIL_MARKER} unexpected:{type(exc).__name__}", file=sys.stderr)
        return 2
    approvals = document["approvals"]
    print(
        f"{BLOCKED_MARKER} design_only=1 pure_lifecycle_oracle=1 "
        f"runtime_owner=0 thermo_recovery=0 "
        f"two_phase_chemistry={int(approvals['two_phase_chemistry_enabled'])} "
        f"phase_aware_chemistry={int(approvals['phase_aware_chemistry_enabled'])} "
        f"numerical_limits_validated={int(approvals['numerical_displacement_limits_validated'])} "
        f"production_ready={int(approvals['production_ready'])} checks={checks}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
