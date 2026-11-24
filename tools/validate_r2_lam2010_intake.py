#!/usr/bin/env python3
"""Validate the internal-only Lam 2010 propane IDT intake."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import stat
import sys
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple


SCHEMA = "rfam.r2-lam2010-internal-intake/v1"
BLOCKED_MARKER = "R2_LAM2010_INTAKE_BLOCKED"
MANIFEST_SHA256 = "9323e023667fb8c82ea40628a8cef991be9d45c956c17e1f5cce43c756ee3157"
SOURCE_PATH = "tmp/lam-2010-proci-supplement-mmc1.doc"
SOURCE_BYTE_COUNT = 115712
SOURCE_SHA256 = "154b690336f9bd2aa6a56c417386d620d5dc360e5e996a8b41900238bd17f977"
RIGHTS_EVIDENCE_PATH = "tmp/crossref-lam-2010.json"
RIGHTS_EVIDENCE_BYTE_COUNT = 5986
RIGHTS_EVIDENCE_SHA256 = "48cbb63189102c227f2da4818f811b14fe5440fdb585bb8b9c18eec97a88207a"
GENERATOR_PATH = "realFluidAllMach/tools/generate_r2_lam2010_intake.py"
GENERATOR_SHA256 = "a6e6cc0c41499d903c0313867a680738325d5f3151b5d9131e3ca91675a2dd72"
CSV_PATH = "realFluidAllMach/tests/reference/r2_lam2010_proci_intake_v1.csv"
CSV_BYTE_COUNT = 7122
CSV_SHA256 = "e87b042e44d7fefee4fdafbd3842bda284afa0475a2a5430b401d4785510574e"
CSV_FIELDS = [
    "table_id",
    "row_in_table",
    "facility_mode",
    "t5_k",
    "p5_atm",
    "dp5_dt_percent_per_ms",
    "dp5_dt_status",
    "idt_ms",
]

EXPECTED_ROOT_KEYS = [
    "schema",
    "contract_version",
    "status",
    "scope",
    "citation",
    "source",
    "extraction",
    "experiment",
    "uncertainty",
    "rights_evidence",
    "rights",
    "eligibility",
    "approvals",
    "limitations",
]

EXPECTED_CITATION = {
    "doi": "10.1016/j.proci.2010.06.131",
    "pii": "S1540748910002609",
    "title": "Shock tube ignition delay time measurements in propane/O2/argon mixtures at near-constant-volume conditions",
    "authors": ["K.-Y. Lam", "Z. Hong", "D. F. Davidson", "R. K. Hanson"],
    "journal": "Proceedings of the Combustion Institute",
    "volume": "33",
    "year": 2011,
}

EXPECTED_SOURCE = {
    "provider": "Elsevier",
    "artifact_role": "article-supplementary-data",
    "url": "https://ars.els-cdn.com/content/image/1-s2.0-S1540748910002609-mmc1.doc",
    "local_path": SOURCE_PATH,
    "media_type": "application/msword",
    "container": "OLE Compound Document Word 97-2003",
    "byte_count": SOURCE_BYTE_COUNT,
    "sha256": SOURCE_SHA256,
    "retrieved_on": "2026-07-28",
    "repository_copy": False,
}

EXPECTED_EXTRACTION = {
    "generator_path": GENERATOR_PATH,
    "generator_sha256": GENERATOR_SHA256,
    "runtime": "Microsoft Word COM read-only table API",
    "csv_path": CSV_PATH,
    "csv_byte_count": CSV_BYTE_COUNT,
    "csv_sha256": CSV_SHA256,
    "canonicalization": "ASCII RFC4180-compatible CSV LF final newline source lexical decimals",
    "fields": CSV_FIELDS,
    "regeneration_count": 2,
    "regenerations_byte_identical": True,
}

EXPECTED_FIELD_CONTRACT = {
    "t5_k": "post-reflected-shock temperature in K",
    "p5_atm": "post-reflected-shock pressure in atm",
    "idt_ms": "ignition delay time in milliseconds with event definition absent from supplement",
    "dp5_dt_percent_per_ms": "reported per-row dP5/dt in percent per millisecond when present",
    "dp5_dt_status": "reported not-available-cell or not-tabulated-caption-approximately-zero",
}

EXPECTED_TABLES = [
    {
        "table_id": "S1",
        "row_count": 17,
        "source_column_count": 4,
        "facility_mode": "conventional-pressure-rise",
        "caption_dp5_dt": "1.2-3.5%/ms",
        "t5_k_range": ["991", "1405"],
        "p5_atm_range": ["5.46", "7.44"],
        "idt_ms_range": ["0.083", "17.79"],
        "reported_dp5_dt_rows": 6,
        "reported_dp5_dt_range": ["1.2", "3.5"],
    },
    {
        "table_id": "S2",
        "row_count": 18,
        "source_column_count": 3,
        "facility_mode": "near-constant-volume",
        "caption_dp5_dt": "approximately-0%/ms",
        "t5_k_range": ["1021", "1381"],
        "p5_atm_range": ["5.26", "6.82"],
        "idt_ms_range": ["0.15", "21.23"],
        "reported_dp5_dt_rows": 0,
        "reported_dp5_dt_range": None,
    },
    {
        "table_id": "S3",
        "row_count": 19,
        "source_column_count": 4,
        "facility_mode": "conventional-pressure-rise",
        "caption_dp5_dt": "5.9-8.0%/ms",
        "t5_k_range": ["950", "1289"],
        "p5_atm_range": ["21.1", "29.5"],
        "idt_ms_range": ["0.20", "6.79"],
        "reported_dp5_dt_rows": 11,
        "reported_dp5_dt_range": ["5.9", "8.0"],
    },
    {
        "table_id": "S4",
        "row_count": 18,
        "source_column_count": 3,
        "facility_mode": "near-constant-volume",
        "caption_dp5_dt": "approximately-0%/ms",
        "t5_k_range": ["988", "1273"],
        "p5_atm_range": ["22.8", "31.1"],
        "idt_ms_range": ["0.26", "9.28"],
        "reported_dp5_dt_rows": 0,
        "reported_dp5_dt_range": None,
    },
    {
        "table_id": "S5",
        "row_count": 10,
        "source_column_count": 4,
        "facility_mode": "conventional-pressure-rise",
        "caption_dp5_dt": "6.2-8.5%/ms",
        "t5_k_range": ["980", "1244"],
        "p5_atm_range": ["54.6", "66.8"],
        "idt_ms_range": ["0.17", "3.26"],
        "reported_dp5_dt_rows": 6,
        "reported_dp5_dt_range": ["6.2", "8.5"],
    },
    {
        "table_id": "S6",
        "row_count": 13,
        "source_column_count": 3,
        "facility_mode": "near-constant-volume",
        "caption_dp5_dt": "approximately-0%/ms",
        "t5_k_range": ["994", "1247"],
        "p5_atm_range": ["53.7", "67.8"],
        "idt_ms_range": ["0.20", "4.75"],
        "reported_dp5_dt_rows": 0,
        "reported_dp5_dt_range": None,
    },
]

EXPECTED_UNCERTAINTY = {
    "t5_available": False,
    "p5_available": False,
    "idt_available": False,
    "composition_available": False,
    "facility_repeatability_available": False,
    "classification": "not-present-in-supplement-tables",
}

EXPECTED_RIGHTS_EVIDENCE = {
    "provider": "Crossref REST API",
    "artifact_role": "rights-metadata-screen-not-license",
    "url": "https://api.crossref.org/works/10.1016/j.proci.2010.06.131",
    "local_path": RIGHTS_EVIDENCE_PATH,
    "media_type": "application/json",
    "byte_count": RIGHTS_EVIDENCE_BYTE_COUNT,
    "sha256": RIGHTS_EVIDENCE_SHA256,
    "retrieved_on": "2026-07-28",
    "observed_license_content_version": "tdm",
    "observed_license_url": "https://www.elsevier.com/tdm/userlicense/1.0/",
    "authorization_effect": "none",
}

EXPECTED_RIGHTS = {
    "status": "unresolved",
    "publicly_reachable": True,
    "crossref_open_content_license_observed": False,
    "crossref_metadata_is_authorization": False,
    "crossref_license_url": "https://www.elsevier.com/tdm/userlicense/1.0/",
    "crossref_content_version": "tdm",
    "internal_screening_only": True,
    "fit_eligible": False,
    "redistribution": False,
    "commercial_or_production_use_authorized": False,
    "derived_csv_storage_scope": "internal-workspace-only-do-not-redistribute",
    "required_resolution": "rights-holder terms or permission covering fit retention derived outputs and intended production use",
}

EXPECTED_ELIGIBILITY = {
    "training_eligible": False,
    "holdout_eligible": False,
    "fit_performed": False,
    "candidate_role_if_resolved": "propane-o2-argon-idt-screen-candidate",
    "blocking_reasons": [
        "fit and retained-derived-output rights are unresolved",
        "the supplement does not state the composition basis Ar fraction or equivalence ratio",
        "the supplement does not state the ignition-delay event definition or timing origin",
        "measurement uncertainties and facility repeatability are absent",
        "training and independent holdout campaigns have not been frozen",
    ],
}

EXPECTED_APPROVALS = {
    "experimental_training_data_selected": False,
    "experimental_holdout_data_selected": False,
    "r2_physical_fit_approved": False,
    "r2_parameter_installation_approved": False,
    "hem_chemistry_enabled": False,
    "fluent_execution_approved": False,
    "production_ready": False,
}

EXPECTED_LIMITATIONS = [
    "composition basis Ar balance and equivalence ratio are not source-stated and no normalized composition is installed",
    "near-constant-volume captions state approximately zero dP5/dt but do not tabulate per-row zero values",
    "N.A. dP5/dt cells remain missing and are not imputed",
    "Crossref TDM metadata is screening evidence and not an open-content or production-use license",
    "the canonical CSV is an internal screening derivative and not a redistributable fit dataset",
    "no chemistry parameter capability Fluent action or production status is changed",
]

_INTEGER_RE = re.compile(r"^[1-9][0-9]*$")
_DECIMAL_RE = re.compile(r"^(?:0|[1-9][0-9]*)(?:\.[0-9]+)?$")


class IntakeError(ValueError):
    """A fail-closed Lam intake violation."""


def _reject_duplicate_pairs(pairs: Sequence[Tuple[str, Any]]) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise IntakeError("duplicate-json-key:%s" % key)
        result[key] = value
    return result


def _reject_nonfinite(token: str) -> None:
    raise IntakeError("nonfinite-json-number:%s" % token)


def _strict_equal(actual: Any, expected: Any, label: str) -> None:
    if type(actual) is not type(expected):
        raise IntakeError("contract-type-mismatch:%s" % label)
    if isinstance(expected, dict):
        actual_keys = list(actual)
        expected_keys = list(expected)
        if actual_keys != expected_keys:
            missing = [key for key in expected_keys if key not in actual]
            unknown = [key for key in actual_keys if key not in expected]
            if missing:
                raise IntakeError("contract-missing-keys:%s:%s" % (label, ",".join(missing)))
            if unknown:
                raise IntakeError("contract-unknown-keys:%s:%s" % (label, ",".join(unknown)))
            raise IntakeError("contract-key-order-mismatch:%s" % label)
        for key in expected_keys:
            _strict_equal(actual[key], expected[key], "%s.%s" % (label, key))
        return
    if isinstance(expected, list):
        if len(actual) != len(expected):
            raise IntakeError("contract-list-length-mismatch:%s" % label)
        for index, item in enumerate(expected):
            _strict_equal(actual[index], item, "%s[%d]" % (label, index))
        return
    if actual != expected:
        raise IntakeError("contract-value-mismatch:%s" % label)


def _decode_manifest(path: Path) -> Dict[str, Any]:
    try:
        payload = path.read_bytes()
    except OSError as exc:
        raise IntakeError("manifest-read-failed:%s" % exc) from exc
    if b"\r" in payload or not payload.endswith(b"\n"):
        raise IntakeError("manifest-not-lf-final-newline")
    try:
        text = payload.decode("ascii")
        document = json.loads(
            text,
            object_pairs_hook=_reject_duplicate_pairs,
            parse_constant=_reject_nonfinite,
        )
    except IntakeError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise IntakeError("manifest-json-invalid:%s" % exc) from exc
    if not isinstance(document, dict):
        raise IntakeError("manifest-root-not-object")
    if list(document) != EXPECTED_ROOT_KEYS:
        missing = [key for key in EXPECTED_ROOT_KEYS if key not in document]
        unknown = [key for key in document if key not in EXPECTED_ROOT_KEYS]
        if missing:
            raise IntakeError("contract-missing-keys:root:%s" % ",".join(missing))
        if unknown:
            raise IntakeError("contract-unknown-keys:root:%s" % ",".join(unknown))
        raise IntakeError("contract-key-order-mismatch:root")

    _strict_equal(document["schema"], SCHEMA, "root.schema")
    _strict_equal(document["contract_version"], 1, "root.contract_version")
    _strict_equal(document["status"], "blocked-internal-screening-only", "root.status")
    _strict_equal(
        document["scope"],
        "six-table-lam2010-propane-idt-offline-intake-no-fit-no-redistribution",
        "root.scope",
    )
    _strict_equal(document["citation"], EXPECTED_CITATION, "root.citation")
    _strict_equal(document["source"], EXPECTED_SOURCE, "root.source")
    _strict_equal(document["extraction"], EXPECTED_EXTRACTION, "root.extraction")

    experiment = document["experiment"]
    expected_experiment_keys = [
        "apparatus",
        "reported_initial_mixture",
        "reported_composition_basis",
        "reported_species_values",
        "normalized_composition_available",
        "equivalence_ratio_reported",
        "nominal_pressure_atm",
        "table_count",
        "row_count",
        "conventional_pressure_rise_rows",
        "near_constant_volume_rows",
        "field_contract",
        "tables",
        "overall_ranges",
    ]
    if not isinstance(experiment, dict):
        raise IntakeError("contract-type-mismatch:root.experiment")
    if list(experiment) != expected_experiment_keys:
        unknown = [key for key in experiment if key not in expected_experiment_keys]
        missing = [key for key in expected_experiment_keys if key not in experiment]
        if missing:
            raise IntakeError("contract-missing-keys:root.experiment:%s" % ",".join(missing))
        if unknown:
            raise IntakeError("contract-unknown-keys:root.experiment:%s" % ",".join(unknown))
        raise IntakeError("contract-key-order-mismatch:root.experiment")
    _strict_equal(experiment["apparatus"], "shock tube", "root.experiment.apparatus")
    _strict_equal(
        experiment["reported_initial_mixture"],
        "0.8% C3H8/ 8% O2/ Ar",
        "root.experiment.reported_initial_mixture",
    )
    _strict_equal(
        experiment["reported_composition_basis"],
        "not-stated-in-supplement",
        "root.experiment.reported_composition_basis",
    )
    _strict_equal(
        experiment["reported_species_values"],
        {"C3H8": "0.8%", "O2": "8%", "Ar": "balance-unspecified"},
        "root.experiment.reported_species_values",
    )
    _strict_equal(
        experiment["normalized_composition_available"],
        False,
        "root.experiment.normalized_composition_available",
    )
    _strict_equal(
        experiment["equivalence_ratio_reported"],
        False,
        "root.experiment.equivalence_ratio_reported",
    )
    _strict_equal(
        experiment["nominal_pressure_atm"],
        ["6", "24", "60"],
        "root.experiment.nominal_pressure_atm",
    )
    _strict_equal(experiment["table_count"], 6, "root.experiment.table_count")
    _strict_equal(experiment["row_count"], 95, "root.experiment.row_count")
    _strict_equal(
        experiment["conventional_pressure_rise_rows"],
        46,
        "root.experiment.conventional_pressure_rise_rows",
    )
    _strict_equal(
        experiment["near_constant_volume_rows"],
        49,
        "root.experiment.near_constant_volume_rows",
    )
    _strict_equal(experiment["field_contract"], EXPECTED_FIELD_CONTRACT, "root.experiment.field_contract")
    _strict_equal(experiment["tables"], EXPECTED_TABLES, "root.experiment.tables")
    _strict_equal(
        experiment["overall_ranges"],
        {
            "t5_k": ["950", "1405"],
            "p5_atm": ["5.26", "67.8"],
            "idt_ms": ["0.083", "21.23"],
            "reported_dp5_dt_percent_per_ms": ["1.2", "8.5"],
        },
        "root.experiment.overall_ranges",
    )
    _strict_equal(document["uncertainty"], EXPECTED_UNCERTAINTY, "root.uncertainty")
    _strict_equal(
        document["rights_evidence"],
        EXPECTED_RIGHTS_EVIDENCE,
        "root.rights_evidence",
    )
    _strict_equal(document["rights"], EXPECTED_RIGHTS, "root.rights")
    _strict_equal(document["eligibility"], EXPECTED_ELIGIBILITY, "root.eligibility")
    _strict_equal(document["approvals"], EXPECTED_APPROVALS, "root.approvals")
    _strict_equal(document["limitations"], EXPECTED_LIMITATIONS, "root.limitations")
    if hashlib.sha256(payload).hexdigest() != MANIFEST_SHA256:
        raise IntakeError("manifest-sha256-mismatch")
    return document


def _resolve_file(root: Path, relative: str, label: str) -> Path:
    if not relative or "\\" in relative or relative.startswith("/"):
        raise IntakeError("path-invalid:%s" % label)
    if re.match(r"^[A-Za-z]:", relative):
        raise IntakeError("path-not-relative:%s" % label)
    parts = relative.split("/")
    if any(part in ("", ".", "..") for part in parts):
        raise IntakeError("path-noncanonical:%s" % label)
    resolved_root = root.resolve()
    candidate = resolved_root.joinpath(*parts)
    resolved = candidate.resolve()
    if resolved != resolved_root and resolved_root not in resolved.parents:
        raise IntakeError("path-escapes-root:%s" % label)
    current = resolved_root
    for part in parts:
        current = current / part
        try:
            info = current.lstat()
        except OSError as exc:
            raise IntakeError("file-stat-failed:%s:%s" % (label, exc)) from exc
        if stat.S_ISLNK(info.st_mode):
            raise IntakeError("symlink-forbidden:%s" % label)
        attributes = getattr(info, "st_file_attributes", 0)
        if attributes & getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400):
            raise IntakeError("reparse-forbidden:%s" % label)
    if not candidate.is_file():
        raise IntakeError("not-regular-file:%s" % label)
    return candidate


def _verify_file(path: Path, size: int, expected_sha: str, label: str) -> bytes:
    try:
        payload = path.read_bytes()
    except OSError as exc:
        raise IntakeError("file-read-failed:%s:%s" % (label, exc)) from exc
    if len(payload) != size:
        raise IntakeError("file-byte-count-mismatch:%s" % label)
    if hashlib.sha256(payload).hexdigest() != expected_sha:
        raise IntakeError("file-sha256-mismatch:%s" % label)
    return payload


def _validate_rights_evidence(payload: bytes) -> None:
    try:
        document = json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=_reject_duplicate_pairs,
            parse_constant=_reject_nonfinite,
        )
    except IntakeError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise IntakeError("rights-evidence-json-invalid:%s" % exc) from exc
    if not isinstance(document, dict):
        raise IntakeError("rights-evidence-root-not-object")
    _strict_equal(document.get("status"), "ok", "rights_evidence.status")
    message = document.get("message")
    if not isinstance(message, dict):
        raise IntakeError("rights-evidence-message-not-object")
    _strict_equal(
        message.get("DOI"),
        EXPECTED_CITATION["doi"],
        "rights_evidence.message.DOI",
    )
    _strict_equal(
        message.get("title"),
        [EXPECTED_CITATION["title"]],
        "rights_evidence.message.title",
    )
    _strict_equal(
        message.get("alternative-id"),
        [EXPECTED_CITATION["pii"]],
        "rights_evidence.message.alternative-id",
    )
    licenses = message.get("license")
    if not isinstance(licenses, list) or len(licenses) != 1:
        raise IntakeError("rights-evidence-license-count-mismatch")
    license_entry = licenses[0]
    if not isinstance(license_entry, dict):
        raise IntakeError("rights-evidence-license-not-object")
    _strict_equal(
        license_entry.get("content-version"),
        "tdm",
        "rights_evidence.message.license[0].content-version",
    )
    _strict_equal(
        license_entry.get("URL"),
        EXPECTED_RIGHTS["crossref_license_url"],
        "rights_evidence.message.license[0].URL",
    )


def _decimal(text: str, integer: bool, label: str) -> Decimal:
    pattern = _INTEGER_RE if integer else _DECIMAL_RE
    if not pattern.fullmatch(text):
        raise IntakeError("csv-number-invalid:%s" % label)
    try:
        value = Decimal(text)
    except InvalidOperation as exc:
        raise IntakeError("csv-number-invalid:%s" % label) from exc
    if not value.is_finite() or value <= 0:
        raise IntakeError("csv-number-nonpositive:%s" % label)
    return value


def _parse_csv(payload: bytes, tables: Sequence[Dict[str, Any]]) -> List[Dict[str, str]]:
    if b"\r" in payload or not payload.endswith(b"\n") or b"\x00" in payload:
        raise IntakeError("csv-not-canonical-lf")
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as exc:
        raise IntakeError("csv-not-ascii") from exc
    try:
        parsed = list(csv.reader(io.StringIO(text, newline=""), strict=True))
    except csv.Error as exc:
        raise IntakeError("csv-parse-failed:%s" % exc) from exc
    if not parsed or parsed[0] != CSV_FIELDS:
        raise IntakeError("csv-header-mismatch")
    if len(parsed) != 96:
        raise IntakeError("csv-row-count-mismatch:%d" % (len(parsed) - 1))

    rows: List[Dict[str, str]] = []
    for line_number, values in enumerate(parsed[1:], start=2):
        if len(values) != len(CSV_FIELDS):
            raise IntakeError("csv-column-count-mismatch:%d" % line_number)
        rows.append(dict(zip(CSV_FIELDS, values)))
    regenerated = io.StringIO(newline="")
    writer = csv.DictWriter(regenerated, fieldnames=CSV_FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    if regenerated.getvalue() != text:
        raise IntakeError("csv-not-canonical-writer-output")

    cursor = 0
    mode_counts = {"conventional-pressure-rise": 0, "near-constant-volume": 0}
    all_t5: List[Decimal] = []
    all_p5: List[Decimal] = []
    all_idt: List[Decimal] = []
    all_dpdt: List[Decimal] = []
    for table in tables:
        count = table["row_count"]
        selected = rows[cursor : cursor + count]
        cursor += count
        if len(selected) != count:
            raise IntakeError("csv-table-row-count-mismatch:%s" % table["table_id"])
        t5_values: List[Decimal] = []
        p5_values: List[Decimal] = []
        idt_values: List[Decimal] = []
        dpdt_values: List[Decimal] = []
        for index, row in enumerate(selected, start=1):
            label = "%s-R%d" % (table["table_id"], index)
            if row["table_id"] != table["table_id"]:
                raise IntakeError("csv-table-order-mismatch:%s" % label)
            if row["row_in_table"] != str(index):
                raise IntakeError("csv-row-index-mismatch:%s" % label)
            if row["facility_mode"] != table["facility_mode"]:
                raise IntakeError("csv-facility-mode-mismatch:%s" % label)
            t5 = _decimal(row["t5_k"], True, label + ".t5_k")
            p5 = _decimal(row["p5_atm"], False, label + ".p5_atm")
            idt = _decimal(row["idt_ms"], False, label + ".idt_ms")
            status = row["dp5_dt_status"]
            dpdt_text = row["dp5_dt_percent_per_ms"]
            if table["facility_mode"] == "near-constant-volume":
                if status != "not-tabulated-caption-approximately-zero" or dpdt_text:
                    raise IntakeError("csv-near-volume-dpdt-mismatch:%s" % label)
            elif status == "not-available-cell":
                if dpdt_text:
                    raise IntakeError("csv-missing-dpdt-has-value:%s" % label)
            elif status == "reported":
                dpdt_values.append(_decimal(dpdt_text, False, label + ".dp5_dt"))
            else:
                raise IntakeError("csv-dpdt-status-invalid:%s" % label)
            t5_values.append(t5)
            p5_values.append(p5)
            idt_values.append(idt)
            mode_counts[row["facility_mode"]] += 1

        expected_ranges = {
            "t5_k": (t5_values, table["t5_k_range"]),
            "p5_atm": (p5_values, table["p5_atm_range"]),
            "idt_ms": (idt_values, table["idt_ms_range"]),
        }
        for field, (values, expected) in expected_ranges.items():
            if [min(values), max(values)] != [Decimal(expected[0]), Decimal(expected[1])]:
                raise IntakeError("csv-table-range-mismatch:%s.%s" % (table["table_id"], field))
        if len(dpdt_values) != table["reported_dp5_dt_rows"]:
            raise IntakeError("csv-table-dpdt-count-mismatch:%s" % table["table_id"])
        expected_dpdt = table["reported_dp5_dt_range"]
        if expected_dpdt is None:
            if dpdt_values:
                raise IntakeError("csv-table-unexpected-dpdt:%s" % table["table_id"])
        elif [min(dpdt_values), max(dpdt_values)] != [
            Decimal(expected_dpdt[0]),
            Decimal(expected_dpdt[1]),
        ]:
            raise IntakeError("csv-table-dpdt-range-mismatch:%s" % table["table_id"])
        all_t5.extend(t5_values)
        all_p5.extend(p5_values)
        all_idt.extend(idt_values)
        all_dpdt.extend(dpdt_values)

    if cursor != len(rows):
        raise IntakeError("csv-unassigned-rows")
    if mode_counts != {"conventional-pressure-rise": 46, "near-constant-volume": 49}:
        raise IntakeError("csv-facility-mode-count-mismatch")
    overall = [
        [min(all_t5), max(all_t5)],
        [min(all_p5), max(all_p5)],
        [min(all_idt), max(all_idt)],
        [min(all_dpdt), max(all_dpdt)],
    ]
    expected_overall = [
        [Decimal("950"), Decimal("1405")],
        [Decimal("5.26"), Decimal("67.8")],
        [Decimal("0.083"), Decimal("21.23")],
        [Decimal("1.2"), Decimal("8.5")],
    ]
    if overall != expected_overall:
        raise IntakeError("csv-overall-range-mismatch")
    return rows


def validate(root: Path, manifest_path: Path) -> Dict[str, Any]:
    manifest = _decode_manifest(manifest_path)
    source = _resolve_file(root, SOURCE_PATH, "source")
    _verify_file(source, SOURCE_BYTE_COUNT, SOURCE_SHA256, "source")
    rights_evidence = _resolve_file(root, RIGHTS_EVIDENCE_PATH, "rights-evidence")
    rights_payload = _verify_file(
        rights_evidence,
        RIGHTS_EVIDENCE_BYTE_COUNT,
        RIGHTS_EVIDENCE_SHA256,
        "rights-evidence",
    )
    _validate_rights_evidence(rights_payload)
    generator = _resolve_file(root, GENERATOR_PATH, "generator")
    _verify_file(generator, generator.stat().st_size, GENERATOR_SHA256, "generator")
    csv_path = _resolve_file(root, CSV_PATH, "csv")
    payload = _verify_file(csv_path, CSV_BYTE_COUNT, CSV_SHA256, "csv")
    rows = _parse_csv(payload, manifest["experiment"]["tables"])
    return {"manifest": manifest, "rows": rows}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    try:
        args = parser.parse_args(argv)
        result = validate(args.root, args.manifest)
        print(
            "%s tables=6 rows=%d internal_screening_only=1 fit_eligible=0 redistribution=0"
            % (BLOCKED_MARKER, len(result["rows"]))
        )
        return 0
    except IntakeError as exc:
        print(
            json.dumps(
                {
                    "schema": "rfam.r2-lam2010-intake-check/v1",
                    "status": "FAIL",
                    "reason": str(exc),
                },
                sort_keys=True,
                separators=(",", ":"),
            ),
            file=sys.stderr,
        )
        return 2
    except Exception as exc:  # pragma: no cover - fail-closed boundary
        print(
            json.dumps(
                {
                    "schema": "rfam.r2-lam2010-intake-check/v1",
                    "status": "FAIL",
                    "reason": "unexpected-%s" % type(exc).__name__,
                },
                sort_keys=True,
                separators=(",", ":"),
            ),
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
