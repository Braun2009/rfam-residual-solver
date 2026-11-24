#!/usr/bin/env python3
"""Validate O2/C3H8 kij inputs and build a no-fit, no-runtime dry-run plan."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import sys
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any
from urllib.parse import urlparse


SCHEMA = "rfam.kij-fit-harness-manifest/v1"
AUDIT_SCHEMA = "rfam.kij-fit-harness-audit/v1"
PASS_MARKER = "RFAM_KIJ_FIT_HARNESS_INPUTS_VALIDATED"
TEST_PASS_MARKER = "RFAM_KIJ_FIT_HARNESS_SYNTHETIC_TEST_ONLY_PASS"
FAIL_MARKER = "RFAM_KIJ_FIT_HARNESS_FAIL"

SYNTHETIC_CLASS = "synthetic-nonphysical-test-only"
PHYSICAL_CLASS = "rights-cleared-experimental-vle-pvt"
SYNTHETIC_WATERMARK = (
    "SYNTHETIC NONPHYSICAL TEST FIXTURE - NEVER USE FOR FIT OR RUNTIME"
)
ROW_SYNTHETIC_WATERMARK = "SYNTHETIC_NONPHYSICAL_TEST_ONLY"

PAIR = {
    "components": ["C3H8", "O2"],
    "canonical_pair_key": "C3H8|O2",
    "canonical_order": "lexicographic-component-name",
    "eos": "Peng-Robinson",
    "mixing_rule": "van-der-Waals-one-fluid",
    "parameter": "symmetric-binary-interaction-kij",
}

UNITS = {
    "composition_basis": "mole-fraction",
    "temperature": "K",
    "pressure": "Pa",
    "mass_density": "kg/m3",
    "uncertainty_basis": "one-standard-uncertainty",
}

APPROVALS = {
    "fit_execution_approved": False,
    "kij_value_validated": False,
    "runtime_kij_emission_approved": False,
    "binary_interaction_db_write_approved": False,
    "hem_enabled": False,
    "fluent_execution_approved": False,
    "hot_start_approved": False,
    "production_ready": False,
}

ROOT_KEYS = (
    "schema",
    "contract_version",
    "status",
    "scope",
    "data_class",
    "watermark",
    "pair",
    "units",
    "rights",
    "artifacts",
    "partition",
    "envelope",
    "fit_plan",
    "approvals",
)
RIGHTS_KEYS = (
    "status",
    "numerical_fit_use_permitted",
    "repository_redistribution_permitted",
    "external_no_copy_workflow",
    "terms_artifact",
    "attribution",
    "modification_notice",
    "confirmed_by",
    "confirmed_on",
)
TERMS_KEYS = ("path", "sha256", "byte_count")
SYNTHETIC_ARTIFACT_KEYS = (
    "artifact_id",
    "path",
    "sha256",
    "byte_count",
    "media_type",
    "split",
    "dataset_id",
    "campaign_id",
)
PHYSICAL_ARTIFACT_KEYS = SYNTHETIC_ARTIFACT_KEYS + (
    "provider",
    "citation_or_dataset_doi",
    "immutable_source_uri",
    "retrieved_on",
    "license_or_terms_url",
    "license_or_terms_artifact",
    "canonical_source_authority_sha256",
    "canonical_campaign_id",
    "canonical_campaign_authority_sha256",
    "rights_disposition",
)
PARTITION_KEYS = (
    "training_dataset_ids",
    "holdout_dataset_ids",
    "record_identity",
    "authority_independence_policy",
    "same_authority_nonoverlap_justification",
    "frozen_before_fit",
)
ENVELOPE_KEYS = (
    "temperature_k",
    "pressure_pa",
    "overall_z_o2",
    "liquid_x_o2",
    "vapor_y_o2",
    "allowed_phase_labels",
    "near_critical_policy",
    "solid_risk_policy",
    "outside_envelope_policy",
)
RANGE_KEYS = ("minimum", "maximum")
FIT_PLAN_KEYS = (
    "mode",
    "functional_form_id",
    "parameter_names",
    "reference_temperature_k",
    "objective",
    "model_backend",
    "fit_execution_enabled",
    "emit_runtime_kij",
    "write_binary_interaction_db",
)
RIGHTS_DISPOSITION_KEYS = (
    "schema",
    "contract_version",
    "status",
    "artifact_id",
    "artifact_sha256",
    "provider",
    "citation_or_dataset_doi",
    "immutable_source_uri",
    "license_or_terms_url",
    "license_or_terms_sha256",
    "decision",
    "reviewed_by",
    "reviewed_on",
)
RIGHTS_DECISION_KEYS = (
    "numerical_fit_use_permitted",
    "repository_redistribution_permitted",
    "external_no_copy_workflow",
)
NONOVERLAP_REVIEW_KEYS = (
    "schema",
    "contract_version",
    "status",
    "training_dataset_ids",
    "holdout_dataset_ids",
    "overlapping_source_authority_sha256",
    "overlapping_campaign_authority_sha256",
    "basis",
    "reviewed_by",
    "reviewed_on",
)

CSV_FIELDS = (
    "watermark",
    "artifact_id",
    "dataset_id",
    "campaign_id",
    "record_id",
    "split",
    "measurement_class",
    "temperature_k",
    "pressure_pa",
    "measured_phase_label",
    "temperature_standard_uncertainty_k",
    "pressure_standard_uncertainty_pa",
    "liquid_x_o2",
    "vapor_y_o2",
    "liquid_x_o2_standard_uncertainty",
    "vapor_y_o2_standard_uncertainty",
    "overall_z_o2",
    "mixture_density_kg_m3",
    "overall_z_o2_standard_uncertainty",
    "mixture_density_standard_uncertainty_kg_m3",
)

ALLOWED_PHASE_LABELS = {
    "single-liquid",
    "single-vapor",
    "vapor-liquid-equilibrium",
}
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


class KijFitHarnessError(ValueError):
    """A structured fail-close harness validation error."""


@dataclass(frozen=True)
class Record:
    artifact_id: str
    dataset_id: str
    campaign_id: str
    record_id: str
    split: str
    measurement_class: str
    temperature_k: float
    pressure_pa: float
    measured_phase_label: str
    measurement_fingerprint: tuple[Any, ...]


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise KijFitHarnessError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _expect_keys(value: Any, expected: tuple[str, ...], label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise KijFitHarnessError(f"type-mismatch:{label}:expected-object")
    actual = tuple(value)
    if actual != expected:
        missing = [key for key in expected if key not in value]
        unknown = [key for key in actual if key not in expected]
        if missing:
            raise KijFitHarnessError(f"missing-keys:{label}:" + ",".join(missing))
        if unknown:
            raise KijFitHarnessError(f"unknown-keys:{label}:" + ",".join(unknown))
        raise KijFitHarnessError(f"key-order-mismatch:{label}")
    return value


def _exact(value: Any, expected: Any, label: str) -> None:
    if type(value) is not type(expected):
        raise KijFitHarnessError(f"value-mismatch:{label}")
    if isinstance(expected, dict):
        if tuple(value) != tuple(expected):
            raise KijFitHarnessError(f"value-mismatch:{label}")
        for key in expected:
            _exact(value[key], expected[key], f"{label}.{key}")
        return
    if isinstance(expected, list):
        if len(value) != len(expected):
            raise KijFitHarnessError(f"value-mismatch:{label}")
        for index, item in enumerate(expected):
            _exact(value[index], item, f"{label}[{index}]")
        return
    if value != expected:
        raise KijFitHarnessError(f"value-mismatch:{label}")


def _nonempty(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise KijFitHarnessError(f"invalid-nonempty-string:{label}")
    return value


def _positive_integer(value: Any, label: str) -> int:
    if type(value) is not int or value <= 0:
        raise KijFitHarnessError(f"invalid-positive-integer:{label}")
    return value


def _hex_digest(value: Any, label: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise KijFitHarnessError(f"invalid-sha256:{label}")
    return value


def _number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise KijFitHarnessError(f"invalid-number:{label}")
    result = float(value)
    if not math.isfinite(result):
        raise KijFitHarnessError(f"nonfinite-number:{label}")
    return result


def _csv_number(value: str, label: str, positive: bool = False) -> float:
    if not value or value != value.strip():
        raise KijFitHarnessError(f"missing-or-malformed-number:{label}")
    try:
        result = float(value)
    except ValueError as exc:
        raise KijFitHarnessError(f"invalid-number:{label}") from exc
    if not math.isfinite(result):
        raise KijFitHarnessError(f"nonfinite-number:{label}")
    if positive and result <= 0.0:
        raise KijFitHarnessError(f"nonpositive-uncertainty:{label}")
    return result


def _fraction(value: str, label: str) -> float:
    result = _csv_number(value, label)
    if result < 0.0 or result > 1.0:
        raise KijFitHarnessError(f"fraction-out-of-range:{label}")
    return result


def _workspace_path(root: Path, relative: Any, label: str) -> Path:
    if not isinstance(relative, str) or not relative or "\\" in relative:
        raise KijFitHarnessError(f"invalid-relative-path:{label}")
    if relative.startswith("/") or re.match(r"^[A-Za-z]:", relative):
        raise KijFitHarnessError(f"path-not-relative:{label}")
    resolved_root = root.resolve()
    candidate = (resolved_root / relative).resolve()
    if candidate != resolved_root and resolved_root not in candidate.parents:
        raise KijFitHarnessError(f"path-escapes-root:{label}")
    return candidate


def _read_bound_file(
    root: Path,
    authority: dict[str, Any],
    label: str,
    expected_suffix: str | None = None,
    authority_keys: tuple[str, ...] = TERMS_KEYS,
) -> tuple[Path, bytes]:
    _expect_keys(authority, authority_keys, label)
    path = _workspace_path(root, authority["path"], f"{label}.path")
    if expected_suffix is not None and path.suffix.lower() != expected_suffix:
        raise KijFitHarnessError(f"file-extension-mismatch:{label}")
    expected_hash = _hex_digest(authority["sha256"], f"{label}.sha256")
    expected_size = _positive_integer(authority["byte_count"], f"{label}.byte_count")
    try:
        payload = path.read_bytes()
    except OSError as exc:
        raise KijFitHarnessError(f"file-read-failed:{label}:{exc}") from exc
    if len(payload) != expected_size:
        raise KijFitHarnessError(f"byte-count-mismatch:{label}")
    if _sha256(payload) != expected_hash:
        raise KijFitHarnessError(f"sha256-mismatch:{label}")
    return path, payload


def _decode_bound_json(payload: bytes, label: str) -> dict[str, Any]:
    if b"\r" in payload or not payload.endswith(b"\n"):
        raise KijFitHarnessError(f"bound-json-not-canonical-lf:{label}")
    try:
        text = payload.decode("ascii")
        document = json.loads(text, object_pairs_hook=_reject_duplicate_pairs)
    except KijFitHarnessError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise KijFitHarnessError(f"bound-json-invalid:{label}:{exc}") from exc
    if not isinstance(document, dict):
        raise KijFitHarnessError(f"bound-json-root-not-object:{label}")
    if text != json.dumps(document, indent=2, ensure_ascii=True) + "\n":
        raise KijFitHarnessError(f"bound-json-not-canonical:{label}")
    return document


def _date(value: Any, label: str) -> str:
    result = _nonempty(value, label)
    if DATE_RE.fullmatch(result) is None:
        raise KijFitHarnessError(f"invalid-date:{label}")
    try:
        date.fromisoformat(result)
    except ValueError as exc:
        raise KijFitHarnessError(f"invalid-date:{label}") from exc
    return result


def _source_uri(value: Any, label: str) -> str:
    result = _nonempty(value, label)
    if any(character.isspace() for character in result):
        raise KijFitHarnessError(f"invalid-uri:{label}")
    parsed = urlparse(result)
    if parsed.scheme not in ("https", "http", "doi", "archive", "s3"):
        raise KijFitHarnessError(f"invalid-uri-scheme:{label}")
    if parsed.scheme in ("https", "http", "s3") and not parsed.netloc:
        raise KijFitHarnessError(f"invalid-uri-authority:{label}")
    if parsed.scheme in ("doi", "archive") and not parsed.path:
        raise KijFitHarnessError(f"invalid-uri-path:{label}")
    return result


def _canonical_authority_sha256(fields: dict[str, str]) -> str:
    payload = (json.dumps(fields, separators=(",", ":"), ensure_ascii=True) + "\n").encode("ascii")
    return _sha256(payload)


def _validate_rights_disposition(
    payload: bytes,
    artifact: dict[str, Any],
    license_terms_sha256: str,
    repository_copy: bool,
    external_no_copy: bool,
) -> None:
    label = f"rights-disposition:{artifact['artifact_id']}"
    document = _expect_keys(
        _decode_bound_json(payload, label), RIGHTS_DISPOSITION_KEYS, label
    )
    _exact(document["schema"], "rfam.kij-artifact-rights-disposition/v1", f"{label}.schema")
    _exact(document["contract_version"], 1, f"{label}.contract_version")
    _exact(document["status"], "approved-for-declared-numerical-fit-use", f"{label}.status")
    _exact(document["artifact_id"], artifact["artifact_id"], f"{label}.artifact_id")
    _exact(document["artifact_sha256"], artifact["sha256"], f"{label}.artifact_sha256")
    for key in (
        "provider",
        "citation_or_dataset_doi",
        "immutable_source_uri",
        "license_or_terms_url",
    ):
        _exact(document[key], artifact[key], f"{label}.{key}")
    _exact(document["license_or_terms_sha256"], license_terms_sha256, f"{label}.license_or_terms_sha256")
    decision = _expect_keys(document["decision"], RIGHTS_DECISION_KEYS, f"{label}.decision")
    _exact(decision["numerical_fit_use_permitted"], True, f"{label}.decision.numerical_fit_use_permitted")
    _exact(decision["repository_redistribution_permitted"], repository_copy, f"{label}.decision.repository_redistribution_permitted")
    _exact(decision["external_no_copy_workflow"], external_no_copy, f"{label}.decision.external_no_copy_workflow")
    _nonempty(document["reviewed_by"], f"{label}.reviewed_by")
    _date(document["reviewed_on"], f"{label}.reviewed_on")


def _validate_physical_artifact_authority(
    root: Path,
    artifact: dict[str, Any],
    repository_copy: bool,
    external_no_copy: bool,
    label: str,
) -> tuple[str, str]:
    provider = _nonempty(artifact["provider"], f"{label}.provider")
    citation = _nonempty(
        artifact["citation_or_dataset_doi"], f"{label}.citation_or_dataset_doi"
    )
    source_uri = _source_uri(
        artifact["immutable_source_uri"], f"{label}.immutable_source_uri"
    )
    _date(artifact["retrieved_on"], f"{label}.retrieved_on")
    license_url = _source_uri(
        artifact["license_or_terms_url"], f"{label}.license_or_terms_url"
    )
    source_sha256 = _canonical_authority_sha256(
        {
            "provider": provider,
            "citation_or_dataset_doi": citation,
            "immutable_source_uri": source_uri,
        }
    )
    _exact(
        artifact["canonical_source_authority_sha256"],
        source_sha256,
        f"{label}.canonical_source_authority_sha256",
    )
    campaign_id = _nonempty(
        artifact["canonical_campaign_id"], f"{label}.canonical_campaign_id"
    )
    campaign_sha256 = _canonical_authority_sha256(
        {
            "canonical_source_authority_sha256": source_sha256,
            "canonical_campaign_id": campaign_id,
        }
    )
    _exact(
        artifact["canonical_campaign_authority_sha256"],
        campaign_sha256,
        f"{label}.canonical_campaign_authority_sha256",
    )

    license_authority = _expect_keys(
        artifact["license_or_terms_artifact"], TERMS_KEYS, f"{label}.license_or_terms_artifact"
    )
    _, license_payload = _read_bound_file(
        root,
        license_authority,
        f"license-terms:{artifact['artifact_id']}",
        authority_keys=TERMS_KEYS,
    )
    del license_payload
    disposition_authority = _expect_keys(
        artifact["rights_disposition"], TERMS_KEYS, f"{label}.rights_disposition"
    )
    _, disposition_payload = _read_bound_file(
        root,
        disposition_authority,
        f"rights-disposition:{artifact['artifact_id']}",
        expected_suffix=".json",
        authority_keys=TERMS_KEYS,
    )
    _validate_rights_disposition(
        disposition_payload,
        artifact,
        license_authority["sha256"],
        repository_copy,
        external_no_copy,
    )
    if license_url != artifact["license_or_terms_url"]:
        raise KijFitHarnessError(f"license-url-normalization-mismatch:{label}")
    return source_sha256, campaign_sha256


def _decode_manifest(path: Path, expected_sha256: str) -> tuple[dict[str, Any], str]:
    expected_sha256 = _hex_digest(expected_sha256, "expected-manifest-sha256")
    try:
        payload = path.read_bytes()
    except OSError as exc:
        raise KijFitHarnessError(f"manifest-read-failed:{exc}") from exc
    actual_sha256 = _sha256(payload)
    if actual_sha256 != expected_sha256:
        raise KijFitHarnessError("manifest-sha256-mismatch")
    if b"\r" in payload or not payload.endswith(b"\n"):
        raise KijFitHarnessError("manifest-not-canonical-lf")
    try:
        text = payload.decode("ascii")
        document = json.loads(text, object_pairs_hook=_reject_duplicate_pairs)
    except KijFitHarnessError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise KijFitHarnessError(f"manifest-json-invalid:{exc}") from exc
    if not isinstance(document, dict):
        raise KijFitHarnessError("manifest-root-not-object")
    if text != json.dumps(document, indent=2, ensure_ascii=True) + "\n":
        raise KijFitHarnessError("manifest-not-canonical-json")
    return document, actual_sha256


def _validate_rights(
    root: Path, rights: Any, data_class: str, allow_synthetic: bool
) -> None:
    rights = _expect_keys(rights, RIGHTS_KEYS, "rights")
    _exact(rights["numerical_fit_use_permitted"], True, "rights.numerical_fit_use_permitted")
    repository_copy = rights["repository_redistribution_permitted"]
    external_no_copy = rights["external_no_copy_workflow"]
    if type(repository_copy) is not bool or type(external_no_copy) is not bool:
        raise KijFitHarnessError("invalid-rights-workflow-flags")
    if repository_copy == external_no_copy:
        raise KijFitHarnessError("rights-workflow-must-select-exactly-one-mode")
    terms = _expect_keys(rights["terms_artifact"], TERMS_KEYS, "rights.terms_artifact")
    _read_bound_file(root, terms, "terms")
    _nonempty(rights["attribution"], "rights.attribution")
    _nonempty(rights["modification_notice"], "rights.modification_notice")
    _nonempty(rights["confirmed_by"], "rights.confirmed_by")
    _date(rights["confirmed_on"], "rights.confirmed_on")

    if data_class == SYNTHETIC_CLASS:
        if not allow_synthetic:
            raise KijFitHarnessError("synthetic-fixture-requires-explicit-test-only-flag")
        _exact(rights["status"], "synthetic-test-fixture-author-rights", "rights.status")
        _exact(repository_copy, True, "rights.repository_redistribution_permitted")
        _exact(external_no_copy, False, "rights.external_no_copy_workflow")
    else:
        _exact(rights["status"], "resolved-for-declared-numerical-fit-use", "rights.status")


def _range(value: Any, label: str, *, positive: bool = False) -> tuple[float, float]:
    value = _expect_keys(value, RANGE_KEYS, label)
    minimum = _number(value["minimum"], f"{label}.minimum")
    maximum = _number(value["maximum"], f"{label}.maximum")
    if minimum >= maximum:
        raise KijFitHarnessError(f"invalid-range:{label}")
    if positive and minimum <= 0.0:
        raise KijFitHarnessError(f"nonpositive-range:{label}")
    return minimum, maximum


def _validate_envelope(value: Any) -> dict[str, tuple[float, float]]:
    envelope = _expect_keys(value, ENVELOPE_KEYS, "envelope")
    ranges = {
        "temperature_k": _range(envelope["temperature_k"], "envelope.temperature_k", positive=True),
        "pressure_pa": _range(envelope["pressure_pa"], "envelope.pressure_pa", positive=True),
        "overall_z_o2": _range(envelope["overall_z_o2"], "envelope.overall_z_o2"),
        "liquid_x_o2": _range(envelope["liquid_x_o2"], "envelope.liquid_x_o2"),
        "vapor_y_o2": _range(envelope["vapor_y_o2"], "envelope.vapor_y_o2"),
    }
    for label in ("overall_z_o2", "liquid_x_o2", "vapor_y_o2"):
        if ranges[label][0] < 0.0 or ranges[label][1] > 1.0:
            raise KijFitHarnessError(f"composition-envelope-out-of-range:{label}")
    phase_labels = envelope["allowed_phase_labels"]
    if not isinstance(phase_labels, list) or not phase_labels:
        raise KijFitHarnessError("invalid-phase-label-list")
    if any(not isinstance(item, str) or item not in ALLOWED_PHASE_LABELS for item in phase_labels):
        raise KijFitHarnessError("invalid-phase-label-enum")
    if phase_labels != sorted(set(phase_labels)):
        raise KijFitHarnessError("phase-labels-not-unique-sorted")
    _exact(envelope["near_critical_policy"], "exclude-and-fail-closed", "envelope.near_critical_policy")
    _exact(envelope["solid_risk_policy"], "exclude-and-fail-closed", "envelope.solid_risk_policy")
    _exact(envelope["outside_envelope_policy"], "fail-closed-no-extrapolation", "envelope.outside_envelope_policy")
    return ranges


def _parse_csv(
    payload: bytes,
    authority: dict[str, Any],
    data_class: str,
    ranges: dict[str, tuple[float, float]],
    allowed_phase_labels: list[str],
) -> list[Record]:
    label = f"artifact:{authority['artifact_id']}"
    if b"\r" in payload or not payload.endswith(b"\n"):
        raise KijFitHarnessError(f"csv-not-canonical-lf:{label}")
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as exc:
        raise KijFitHarnessError(f"csv-not-ascii:{label}") from exc
    try:
        rows = list(csv.reader(io.StringIO(text), strict=True))
    except csv.Error as exc:
        raise KijFitHarnessError(f"csv-invalid:{label}:{exc}") from exc
    if not rows or tuple(rows[0]) != CSV_FIELDS:
        raise KijFitHarnessError(f"csv-header-mismatch:{label}")
    if len(set(rows[0])) != len(rows[0]):
        raise KijFitHarnessError(f"csv-duplicate-header:{label}")
    if len(rows) < 2:
        raise KijFitHarnessError(f"csv-no-records:{label}")

    records: list[Record] = []
    for row_index, cells in enumerate(rows[1:], start=2):
        row_label = f"{label}:row-{row_index}"
        if len(cells) != len(CSV_FIELDS):
            raise KijFitHarnessError(f"csv-column-count-mismatch:{row_label}")
        row = dict(zip(CSV_FIELDS, cells))
        expected_watermark = (
            ROW_SYNTHETIC_WATERMARK if data_class == SYNTHETIC_CLASS else "RIGHTS_CLEARED_EXPERIMENTAL_DATA"
        )
        _exact(row["watermark"], expected_watermark, f"{row_label}.watermark")
        for key in ("artifact_id", "dataset_id", "campaign_id", "split"):
            _exact(row[key], authority[key], f"{row_label}.{key}")
        _nonempty(row["record_id"], f"{row_label}.record_id")
        measurement_class = row["measurement_class"]
        if measurement_class not in ("binary-vle", "binary-pvt"):
            raise KijFitHarnessError(f"invalid-measurement-class:{row_label}")
        temperature = _csv_number(row["temperature_k"], f"{row_label}.temperature_k")
        pressure = _csv_number(row["pressure_pa"], f"{row_label}.pressure_pa")
        _csv_number(row["temperature_standard_uncertainty_k"], f"{row_label}.temperature_standard_uncertainty_k", positive=True)
        _csv_number(row["pressure_standard_uncertainty_pa"], f"{row_label}.pressure_standard_uncertainty_pa", positive=True)
        phase_label = row["measured_phase_label"]
        if not phase_label:
            raise KijFitHarnessError(f"missing-phase-label:{row_label}")
        if phase_label not in allowed_phase_labels:
            raise KijFitHarnessError(f"phase-label-outside-envelope:{row_label}")
        if not (ranges["temperature_k"][0] <= temperature <= ranges["temperature_k"][1]):
            raise KijFitHarnessError(f"temperature-outside-envelope:{row_label}")
        if not (ranges["pressure_pa"][0] <= pressure <= ranges["pressure_pa"][1]):
            raise KijFitHarnessError(f"pressure-outside-envelope:{row_label}")

        values: tuple[Any, ...]
        if measurement_class == "binary-vle":
            if phase_label != "vapor-liquid-equilibrium":
                raise KijFitHarnessError(f"vle-phase-label-mismatch:{row_label}")
            x_o2 = _fraction(row["liquid_x_o2"], f"{row_label}.liquid_x_o2")
            y_o2 = _fraction(row["vapor_y_o2"], f"{row_label}.vapor_y_o2")
            _csv_number(row["liquid_x_o2_standard_uncertainty"], f"{row_label}.liquid_x_o2_standard_uncertainty", positive=True)
            _csv_number(row["vapor_y_o2_standard_uncertainty"], f"{row_label}.vapor_y_o2_standard_uncertainty", positive=True)
            if not (ranges["liquid_x_o2"][0] <= x_o2 <= ranges["liquid_x_o2"][1]):
                raise KijFitHarnessError(f"liquid-composition-outside-envelope:{row_label}")
            if not (ranges["vapor_y_o2"][0] <= y_o2 <= ranges["vapor_y_o2"][1]):
                raise KijFitHarnessError(f"vapor-composition-outside-envelope:{row_label}")
            for key in (
                "overall_z_o2",
                "mixture_density_kg_m3",
                "overall_z_o2_standard_uncertainty",
                "mixture_density_standard_uncertainty_kg_m3",
            ):
                if row[key]:
                    raise KijFitHarnessError(f"inactive-pvt-field-not-empty:{row_label}.{key}")
            values = (measurement_class, temperature, pressure, phase_label, x_o2, y_o2)
        else:
            if phase_label not in ("single-liquid", "single-vapor"):
                raise KijFitHarnessError(f"pvt-phase-label-mismatch:{row_label}")
            z_o2 = _fraction(row["overall_z_o2"], f"{row_label}.overall_z_o2")
            density = _csv_number(row["mixture_density_kg_m3"], f"{row_label}.mixture_density_kg_m3")
            if density <= 0.0:
                raise KijFitHarnessError(f"nonpositive-density:{row_label}")
            _csv_number(row["overall_z_o2_standard_uncertainty"], f"{row_label}.overall_z_o2_standard_uncertainty", positive=True)
            _csv_number(row["mixture_density_standard_uncertainty_kg_m3"], f"{row_label}.mixture_density_standard_uncertainty_kg_m3", positive=True)
            if not (ranges["overall_z_o2"][0] <= z_o2 <= ranges["overall_z_o2"][1]):
                raise KijFitHarnessError(f"overall-composition-outside-envelope:{row_label}")
            for key in (
                "liquid_x_o2",
                "vapor_y_o2",
                "liquid_x_o2_standard_uncertainty",
                "vapor_y_o2_standard_uncertainty",
            ):
                if row[key]:
                    raise KijFitHarnessError(f"inactive-vle-field-not-empty:{row_label}.{key}")
            values = (measurement_class, temperature, pressure, phase_label, z_o2, density)

        records.append(
            Record(
                artifact_id=authority["artifact_id"],
                dataset_id=authority["dataset_id"],
                campaign_id=authority["campaign_id"],
                record_id=row["record_id"],
                split=authority["split"],
                measurement_class=measurement_class,
                temperature_k=temperature,
                pressure_pa=pressure,
                measured_phase_label=phase_label,
                measurement_fingerprint=values,
            )
        )
    return records


def _validate_nonoverlap_justification(
    root: Path,
    authority: Any,
    training_ids: list[str],
    holdout_ids: list[str],
    overlapping_sources: list[str],
    overlapping_campaigns: list[str],
) -> None:
    if not overlapping_sources and not overlapping_campaigns:
        _exact(authority, None, "partition.same_authority_nonoverlap_justification")
        return
    descriptor = _expect_keys(
        authority,
        TERMS_KEYS,
        "partition.same_authority_nonoverlap_justification",
    )
    _, payload = _read_bound_file(
        root,
        descriptor,
        "partition-nonoverlap-justification",
        expected_suffix=".json",
        authority_keys=TERMS_KEYS,
    )
    label = "partition-nonoverlap-review"
    document = _expect_keys(
        _decode_bound_json(payload, label), NONOVERLAP_REVIEW_KEYS, label
    )
    _exact(document["schema"], "rfam.kij-partition-nonoverlap-review/v1", f"{label}.schema")
    _exact(document["contract_version"], 1, f"{label}.contract_version")
    _exact(document["status"], "approved-independent-nonoverlap", f"{label}.status")
    _exact(document["training_dataset_ids"], training_ids, f"{label}.training_dataset_ids")
    _exact(document["holdout_dataset_ids"], holdout_ids, f"{label}.holdout_dataset_ids")
    _exact(
        document["overlapping_source_authority_sha256"],
        overlapping_sources,
        f"{label}.overlapping_source_authority_sha256",
    )
    _exact(
        document["overlapping_campaign_authority_sha256"],
        overlapping_campaigns,
        f"{label}.overlapping_campaign_authority_sha256",
    )
    _nonempty(document["basis"], f"{label}.basis")
    _nonempty(document["reviewed_by"], f"{label}.reviewed_by")
    _date(document["reviewed_on"], f"{label}.reviewed_on")


def _validate_artifacts_and_partition(
    root: Path,
    artifacts_value: Any,
    partition_value: Any,
    data_class: str,
    ranges: dict[str, tuple[float, float]],
    allowed_phase_labels: list[str],
    repository_copy: bool,
    external_no_copy: bool,
) -> list[Record]:
    if not isinstance(artifacts_value, list) or len(artifacts_value) < 2:
        raise KijFitHarnessError("artifacts-require-training-and-holdout")
    artifacts: list[dict[str, Any]] = []
    source_authorities: dict[str, str] = {}
    campaign_authorities: dict[str, str] = {}
    artifact_keys = (
        PHYSICAL_ARTIFACT_KEYS
        if data_class == PHYSICAL_CLASS
        else SYNTHETIC_ARTIFACT_KEYS
    )
    for index, item in enumerate(artifacts_value):
        label = f"artifacts[{index}]"
        artifact = _expect_keys(item, artifact_keys, label)
        for key in ("artifact_id", "dataset_id", "campaign_id"):
            _nonempty(artifact[key], f"{label}.{key}")
        if artifact["split"] not in ("training", "holdout"):
            raise KijFitHarnessError(f"invalid-split:{label}")
        _exact(artifact["media_type"], "text/csv", f"{label}.media_type")
        _hex_digest(artifact["sha256"], f"{label}.sha256")
        _positive_integer(artifact["byte_count"], f"{label}.byte_count")
        if data_class == PHYSICAL_CLASS:
            _exact(
                artifact["campaign_id"],
                artifact["canonical_campaign_id"],
                f"{label}.campaign_id",
            )
            source_sha256, campaign_sha256 = _validate_physical_artifact_authority(
                root,
                artifact,
                repository_copy,
                external_no_copy,
                label,
            )
            source_authorities[artifact["artifact_id"]] = source_sha256
            campaign_authorities[artifact["artifact_id"]] = campaign_sha256
        artifacts.append(artifact)
    for key in ("artifact_id", "dataset_id"):
        values = [item[key] for item in artifacts]
        if len(values) != len(set(values)):
            raise KijFitHarnessError(f"duplicate-artifact-authority:{key}")

    partition = _expect_keys(partition_value, PARTITION_KEYS, "partition")
    _exact(partition["record_identity"], "dataset_id+record_id", "partition.record_identity")
    _exact(
        partition["authority_independence_policy"],
        "canonical-source-and-campaign-disjoint-or-reviewed-nonoverlap-justification",
        "partition.authority_independence_policy",
    )
    _exact(partition["frozen_before_fit"], True, "partition.frozen_before_fit")
    training_ids = partition["training_dataset_ids"]
    holdout_ids = partition["holdout_dataset_ids"]
    if not isinstance(training_ids, list) or not training_ids:
        raise KijFitHarnessError("missing-training-datasets")
    if not isinstance(holdout_ids, list) or not holdout_ids:
        raise KijFitHarnessError("missing-holdout-datasets")
    if any(not isinstance(item, str) or not item for item in training_ids + holdout_ids):
        raise KijFitHarnessError("invalid-partition-dataset-id")
    if training_ids != sorted(set(training_ids)) or holdout_ids != sorted(set(holdout_ids)):
        raise KijFitHarnessError("partition-dataset-ids-not-unique-sorted")
    if set(training_ids) & set(holdout_ids):
        raise KijFitHarnessError("training-holdout-dataset-overlap")
    actual_training = sorted(item["dataset_id"] for item in artifacts if item["split"] == "training")
    actual_holdout = sorted(item["dataset_id"] for item in artifacts if item["split"] == "holdout")
    if actual_training != training_ids:
        raise KijFitHarnessError("training-partition-authority-mismatch")
    if actual_holdout != holdout_ids:
        raise KijFitHarnessError("holdout-partition-authority-mismatch")
    if data_class == PHYSICAL_CLASS:
        training_sources = {
            source_authorities[item["artifact_id"]]
            for item in artifacts
            if item["split"] == "training"
        }
        holdout_sources = {
            source_authorities[item["artifact_id"]]
            for item in artifacts
            if item["split"] == "holdout"
        }
        training_campaigns = {
            campaign_authorities[item["artifact_id"]]
            for item in artifacts
            if item["split"] == "training"
        }
        holdout_campaigns = {
            campaign_authorities[item["artifact_id"]]
            for item in artifacts
            if item["split"] == "holdout"
        }
        _validate_nonoverlap_justification(
            root,
            partition["same_authority_nonoverlap_justification"],
            training_ids,
            holdout_ids,
            sorted(training_sources & holdout_sources),
            sorted(training_campaigns & holdout_campaigns),
        )
    else:
        _exact(
            partition["same_authority_nonoverlap_justification"],
            None,
            "partition.same_authority_nonoverlap_justification",
        )
        training_campaigns = {
            item["campaign_id"] for item in artifacts if item["split"] == "training"
        }
        holdout_campaigns = {
            item["campaign_id"] for item in artifacts if item["split"] == "holdout"
        }
        if training_campaigns & holdout_campaigns:
            raise KijFitHarnessError("training-holdout-campaign-overlap")

    records: list[Record] = []
    for artifact in artifacts:
        _, payload = _read_bound_file(
            root,
            artifact,
            f"artifact:{artifact['artifact_id']}",
            ".csv",
            authority_keys=artifact_keys,
        )
        records.extend(
            _parse_csv(payload, artifact, data_class, ranges, allowed_phase_labels)
        )
    identities: set[tuple[str, str]] = set()
    training_fingerprints: set[tuple[Any, ...]] = set()
    holdout_fingerprints: set[tuple[Any, ...]] = set()
    for record in records:
        identity = (record.dataset_id, record.record_id)
        if identity in identities:
            raise KijFitHarnessError("duplicate-record-identity")
        identities.add(identity)
        target = training_fingerprints if record.split == "training" else holdout_fingerprints
        target.add(record.measurement_fingerprint)
    if training_fingerprints & holdout_fingerprints:
        raise KijFitHarnessError("training-holdout-measurement-overlap")
    return records


def _validate_fit_plan(value: Any, records: list[Record], ranges: dict[str, tuple[float, float]]) -> None:
    plan = _expect_keys(value, FIT_PLAN_KEYS, "fit_plan")
    _exact(plan["mode"], "validate-inputs-and-build-plan-only", "fit_plan.mode")
    form = plan["functional_form_id"]
    if form not in ("constant", "reduced-temperature-linear"):
        raise KijFitHarnessError("unsupported-fit-form")
    expected_names = ["c0"] if form == "constant" else ["c0", "c1"]
    _exact(plan["parameter_names"], expected_names, "fit_plan.parameter_names")
    reference = plan["reference_temperature_k"]
    if form == "constant":
        _exact(reference, None, "fit_plan.reference_temperature_k")
    else:
        reference_value = _number(reference, "fit_plan.reference_temperature_k")
        if not ranges["temperature_k"][0] <= reference_value <= ranges["temperature_k"][1]:
            raise KijFitHarnessError("reference-temperature-outside-envelope")
        training_temperatures = {record.temperature_k for record in records if record.split == "training"}
        if len(training_temperatures) < 2:
            raise KijFitHarnessError("temperature-dependent-form-rank-deficient")
    _exact(plan["objective"], "uncertainty-aware-vle-pvt-residual-plan", "fit_plan.objective")
    _exact(plan["model_backend"], "not-connected-dry-run-scaffold", "fit_plan.model_backend")
    _exact(plan["fit_execution_enabled"], False, "fit_plan.fit_execution_enabled")
    _exact(plan["emit_runtime_kij"], False, "fit_plan.emit_runtime_kij")
    _exact(plan["write_binary_interaction_db"], False, "fit_plan.write_binary_interaction_db")


def validate_manifest(
    root: Path,
    manifest_path: Path,
    expected_manifest_sha256: str,
    allow_synthetic_test_fixture: bool = False,
) -> dict[str, Any]:
    document, manifest_sha256 = _decode_manifest(manifest_path, expected_manifest_sha256)
    document = _expect_keys(document, ROOT_KEYS, "root")
    _exact(document["schema"], SCHEMA, "schema")
    _exact(document["contract_version"], 1, "contract_version")
    _exact(document["status"], "input-candidate-dry-run-only", "status")
    _exact(document["scope"], "offline-o2-c3h8-input-validation-and-fit-plan-no-fit-no-runtime", "scope")
    data_class = document["data_class"]
    if data_class not in (PHYSICAL_CLASS, SYNTHETIC_CLASS):
        raise KijFitHarnessError("invalid-data-class")
    if data_class == SYNTHETIC_CLASS:
        _exact(document["watermark"], SYNTHETIC_WATERMARK, "watermark")
    else:
        _exact(document["watermark"], "RIGHTS-CLEARED EXPERIMENTAL INPUT CANDIDATE", "watermark")
    _exact(document["pair"], PAIR, "pair")
    _exact(document["units"], UNITS, "units")
    _validate_rights(root, document["rights"], data_class, allow_synthetic_test_fixture)
    ranges = _validate_envelope(document["envelope"])
    records = _validate_artifacts_and_partition(
        root,
        document["artifacts"],
        document["partition"],
        data_class,
        ranges,
        document["envelope"]["allowed_phase_labels"],
        document["rights"]["repository_redistribution_permitted"],
        document["rights"]["external_no_copy_workflow"],
    )
    _validate_fit_plan(document["fit_plan"], records, ranges)
    _exact(document["approvals"], APPROVALS, "approvals")

    training_count = sum(record.split == "training" for record in records)
    holdout_count = sum(record.split == "holdout" for record in records)
    if training_count == 0 or holdout_count == 0:
        raise KijFitHarnessError("empty-training-or-holdout-partition")
    return {
        "schema": AUDIT_SCHEMA,
        "status": "synthetic-test-only-inputs-validated" if data_class == SYNTHETIC_CLASS else "physical-inputs-structurally-validated",
        "manifest_sha256": manifest_sha256,
        "data_class": data_class,
        "pair_key": PAIR["canonical_pair_key"],
        "training_records": training_count,
        "holdout_records": holdout_count,
        "functional_form_id": document["fit_plan"]["functional_form_id"],
        "fit_executed": False,
        "parameter_values_emitted": False,
        "runtime_kij_emitted": False,
        "binary_interaction_db_modified": False,
        "physical_validation_claimed": False,
        "production_ready": False,
    }


def _write_audit(path: Path, audit: dict[str, Any]) -> None:
    protected_directory = Path(__file__).resolve().parents[1] / "src/thermo/MixingRules"
    protected_paths = {
        (protected_directory / "BinaryInteractionDB.C").resolve(),
        (protected_directory / "BinaryInteractionDB.H").resolve(),
    }
    resolved = path.resolve()
    if resolved in protected_paths:
        raise KijFitHarnessError("audit-output-targets-binary-interaction-db")
    if resolved.suffix.lower() != ".json":
        raise KijFitHarnessError("audit-output-must-be-json")
    if resolved.exists():
        raise KijFitHarnessError("audit-output-already-exists")
    payload = (json.dumps(audit, indent=2, ensure_ascii=True) + "\n").encode("ascii")
    resolved.parent.mkdir(parents=True, exist_ok=True)
    with resolved.open("xb") as stream:
        stream.write(payload)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--expected-manifest-sha256", required=True)
    parser.add_argument("--allow-synthetic-test-fixture", action="store_true")
    parser.add_argument("--audit-output", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        root = args.root.resolve()
        manifest = args.manifest
        if not manifest.is_absolute():
            manifest = (Path.cwd() / manifest).resolve()
        audit = validate_manifest(
            root,
            manifest,
            args.expected_manifest_sha256,
            args.allow_synthetic_test_fixture,
        )
        if args.audit_output is not None:
            output = args.audit_output
            if not output.is_absolute():
                output = (Path.cwd() / output).resolve()
            _write_audit(output, audit)
    except Exception as exc:  # All malformed or unexpected input fails closed.
        reason = str(exc) if isinstance(exc, KijFitHarnessError) else f"unexpected:{type(exc).__name__}:{exc}"
        print(f"{FAIL_MARKER} reason={reason.replace(chr(13), ' ').replace(chr(10), ' ')}")
        return 2

    marker = TEST_PASS_MARKER if audit["data_class"] == SYNTHETIC_CLASS else PASS_MARKER
    print(
        f"{marker} pair={audit['pair_key']} training={audit['training_records']} "
        f"holdout={audit['holdout_records']} fit_executed=0 runtime_kij_emitted=0 "
        "binary_interaction_db_modified=0 production_ready=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
