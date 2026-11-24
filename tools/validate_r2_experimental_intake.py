#!/usr/bin/env python3
"""Validate and extract the bounded SciExpeM R2 experimental-data screen."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
import sys
import tempfile
import xml.etree.ElementTree as ET
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple


SCHEMA = "rfam.r2-experimental-intake/v1"
RESULT_SCHEMA = "rfam.r2-experimental-intake-result/v1"
BLOCKED_MARKER = "R2_EXPERIMENTAL_INTAKE_BLOCKED"
FAIL_MARKER = "R2_EXPERIMENTAL_INTAKE_FAIL"

PUBLICATION_DOI = "10.1016/j.combustflame.2013.05.012"
PUBLICATION_TITLE = (
    "Further study on the ignition delay times of "
    "propane-hydrogen-oxygen-argon mixtures: Effect of equivalence ratio"
)


def _source(record_id: int, byte_count: int, sha256: str, md5: str) -> Dict[str, Any]:
    doi = "10.5281/zenodo.%d" % record_id
    filename = "10.5281-zenodo.%d.xml" % record_id
    return {
        "record_id": record_id,
        "record_doi": doi,
        "record_api_url": "https://zenodo.org/api/records/%d" % record_id,
        "source_url": (
            "https://zenodo.org/api/records/%d/files/"
            "10.5281/zenodo.%d.xml/content" % (record_id, record_id)
        ),
        "repository_path": (
            "realFluidAllMach/tests/reference/r2_experimental_sciexpem_v1/"
            + filename
        ),
        "media_type": "application/xml",
        "license_id": "cc-by-4.0",
        "retrieved_on": "2026-07-28",
        "byte_count": byte_count,
        "sha256": sha256,
        "zenodo_md5": md5,
    }


EXPECTED_SOURCES = [
    _source(
        4639716,
        2328,
        "efba3accb488551674d242b48996ece4aacbeda8788b26c9ca6269e2da187a0c",
        "f6ada662e1ae2348678b05e8b84baccd",
    ),
    _source(
        4639727,
        2327,
        "14b1b1728de3ac1a6341c8b684b07dc6bcb2292e842597b59f1db46f8402d23b",
        "0f5fe91186cb91806c34400645d11c7e",
    ),
    _source(
        4639733,
        2249,
        "06120c83868b48cd4fc5e828864e0b7b0e8bbda53d674c91bc4c3283cac64f58",
        "3b1d6e25a1c8285fdf0301f136064c96",
    ),
    _source(
        4639741,
        2092,
        "c751adafa3090a8775c0e39c0f8e97be4c4627617fb57cf61cb53ce8163bfcbd",
        "9fa7b6c05939afe99e118ce4d3b19aed",
    ),
    _source(
        4639751,
        2171,
        "25a28bd2d0c566c382e6c84ef088bc446c4e33801d9ed121da735ba58ae42209",
        "e61f6213cadf2e482036f0ba3e5ae227",
    ),
    _source(
        4639759,
        2331,
        "71cf410dcf41ab32e1f342a5d8989e37d1967cf17ec6f3e37d0a775896e45bde",
        "ea26dd07ba6bb21e7bd8d3d7fac86157",
    ),
]

EXPECTED_CAMPAIGN = {
    "dataset_family": "SciExpeM-ReSpecTh",
    "publication_doi": PUBLICATION_DOI,
    "publication_title": PUBLICATION_TITLE,
    "dataset_author": "Politecnico di Milano, Italy",
    "experiment_type": "ignition delay measurement",
    "apparatus": "shock tube",
    "record_count": 6,
    "independent_campaign_count": 1,
}

EXPECTED_LICENSE = {
    "provider": "Zenodo",
    "record_metadata_observed_on": "2026-07-28",
    "record_license_id": "cc-by-4.0",
    "spdx_id": "CC-BY-4.0",
    "license_url": "https://creativecommons.org/licenses/by/4.0/",
    "attribution_required": True,
    "legal_review_performed": False,
    "disposition": "record-metadata-declares-cc-by-4.0-retain-attribution",
}

EXPECTED_STRUCTURE = {
    "respec_th_version": "2.3",
    "pressure_units": "atm",
    "pressure_source_type": "reported",
    "composition_units": "mole fraction",
    "composition_source_type": "reported",
    "temperature_units": "K",
    "temperature_source_type": "digitized",
    "ignition_delay_units": "us",
    "ignition_delay_source_type": "digitized",
    "ignition_target": "OH",
    "ignition_method": "d/dt max",
    "species_set": ["AR", "H2", "O2"],
    "pressure_values_atm": ["1.2000000000", "10.0000000000"],
    "curve_count": 6,
    "point_count": 48,
    "c3h8_curve_count": 0,
}

EXPECTED_UNCERTAINTY = {
    "pressure": False,
    "temperature": False,
    "composition": False,
    "ignition_delay": False,
    "classification": "none-present-in-six-xml-artifacts",
}

EXPECTED_EXTRACTION = {
    "result_schema": RESULT_SCHEMA,
    "canonicalization": "ASCII JSON indent=2 LF final newline lexical source decimals",
    "expected_byte_count": 14874,
    "expected_sha256": "79b74c88c78222246c489b3e3e66609131616bae4aebd3e975c201b53a527ab8",
}

EXPECTED_ELIGIBILITY = {
    "source_role": "open-data-discovery-negative-for-propane",
    "training_eligible": False,
    "holdout_eligible": False,
    "fit_performed": False,
    "blocking_reasons": [
        "none of the six initial compositions contains C3H8",
        (
            "temperature and ignition-delay values are digitized rather than "
            "reported raw measurements"
        ),
        (
            "pressure temperature composition and ignition-delay uncertainties "
            "are absent"
        ),
        (
            "all six records cite one publication and cannot provide an "
            "independent campaign holdout"
        ),
    ],
    "required_next": [
        (
            "continue the official-source search for C3H8-containing "
            "ignition-delay records"
        ),
        "obtain a separate rights-cleared experimental campaign for holdout",
        (
            "freeze facility corrections uncertainties partitions metrics and "
            "tolerances before physical fitting"
        ),
    ],
}

EXPECTED_APPROVALS = {
    "experimental_training_data_selected": False,
    "experimental_holdout_data_selected": False,
    "r2_physical_fit_approved": False,
    "r2_parameter_publication_approved": False,
    "hem_chemistry_enabled": False,
    "fluent_execution_approved": False,
    "production_ready": False,
}

EXPECTED_LIMITATIONS = [
    (
        "the six-record list is a bounded source screen and is not a complete "
        "search of SciExpeM or Zenodo"
    ),
    (
        "the publication title mentioning propane does not establish propane "
        "in any individual data artifact"
    ),
    (
        "CC BY 4.0 record metadata does not repair missing scientific "
        "observables or uncertainty fields"
    ),
    (
        "no chemistry parameter model capability Fluent action or production "
        "status is changed"
    ),
]

EXPECTED_MANIFEST = {
    "schema": SCHEMA,
    "contract_version": 1,
    "status": "blocked-no-propane-curves",
    "scope": "offline-six-record-sciexpem-source-intake-no-fit-no-parameter-change",
    "campaign": EXPECTED_CAMPAIGN,
    "license": EXPECTED_LICENSE,
    "sources": EXPECTED_SOURCES,
    "observed_structure": EXPECTED_STRUCTURE,
    "uncertainty_availability": EXPECTED_UNCERTAINTY,
    "extraction": EXPECTED_EXTRACTION,
    "eligibility": EXPECTED_ELIGIBILITY,
    "approvals": EXPECTED_APPROVALS,
    "limitations": EXPECTED_LIMITATIONS,
}

_DECIMAL_RE = re.compile(r"^(?:0|[1-9][0-9]*)(?:\.[0-9]+)?$")


class IntakeError(ValueError):
    """A fail-closed manifest, source, or extraction violation."""


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
                raise IntakeError(
                    "contract-missing-keys:%s:%s" % (label, ",".join(missing))
                )
            if unknown:
                raise IntakeError(
                    "contract-unknown-keys:%s:%s" % (label, ",".join(unknown))
                )
            raise IntakeError("contract-key-order-mismatch:%s" % label)
        for key in expected_keys:
            _strict_equal(actual[key], expected[key], "%s.%s" % (label, key))
        return
    if isinstance(expected, list):
        if len(actual) != len(expected):
            raise IntakeError("contract-list-length-mismatch:%s" % label)
        for index, expected_item in enumerate(expected):
            _strict_equal(actual[index], expected_item, "%s[%d]" % (label, index))
        return
    if actual != expected:
        raise IntakeError("contract-value-mismatch:%s" % label)


def _decode_manifest(path: Path) -> Dict[str, Any]:
    try:
        payload = path.read_bytes()
    except OSError as exc:
        raise IntakeError("manifest-read-failed:%s" % exc) from exc
    if b"\r" in payload or not payload.endswith(b"\n"):
        raise IntakeError("manifest-not-canonical-lf")
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
    canonical = json.dumps(document, indent=2, ensure_ascii=True) + "\n"
    if text != canonical:
        raise IntakeError("manifest-not-canonical-json")
    _strict_equal(document, EXPECTED_MANIFEST, "root")
    return document


def _resolve_source(root: Path, relative: str, label: str) -> Path:
    if not relative or "\\" in relative or relative.startswith("/"):
        raise IntakeError("source-path-invalid:%s" % label)
    if re.match(r"^[A-Za-z]:", relative):
        raise IntakeError("source-path-not-relative:%s" % label)
    parts = relative.split("/")
    if any(part in ("", ".", "..") for part in parts):
        raise IntakeError("source-path-noncanonical:%s" % label)
    resolved_root = root.resolve()
    candidate = resolved_root.joinpath(*parts)
    resolved_candidate = candidate.resolve()
    if resolved_candidate != resolved_root and resolved_root not in resolved_candidate.parents:
        raise IntakeError("source-path-escapes-root:%s" % label)

    current = resolved_root
    for part in parts:
        current = current / part
        try:
            info = current.lstat()
        except OSError as exc:
            raise IntakeError("source-stat-failed:%s:%s" % (label, exc)) from exc
        if stat.S_ISLNK(info.st_mode):
            raise IntakeError("source-symlink-forbidden:%s" % label)
        attributes = getattr(info, "st_file_attributes", 0)
        reparse_flag = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)
        if attributes & reparse_flag:
            raise IntakeError("source-reparse-forbidden:%s" % label)
    if not candidate.is_file():
        raise IntakeError("source-not-regular-file:%s" % label)
    return candidate


def _local_name(element: ET.Element) -> str:
    if not isinstance(element.tag, str):
        raise IntakeError("source-non-element-node")
    if "}" in element.tag:
        return element.tag.rsplit("}", 1)[1]
    return element.tag


def _expect_attrs(element: ET.Element, expected: Dict[str, str], label: str) -> None:
    if element.attrib != expected:
        raise IntakeError("source-attribute-contract-mismatch:%s" % label)


def _expect_children(element: ET.Element, expected: Sequence[str], label: str) -> None:
    actual = [_local_name(child) for child in list(element)]
    if actual != list(expected):
        raise IntakeError("source-child-contract-mismatch:%s" % label)


def _text(element: ET.Element, label: str) -> str:
    if element.text is None or not element.text.strip():
        raise IntakeError("source-missing-text:%s" % label)
    return element.text.strip()


def _leaf_text(element: ET.Element, label: str) -> str:
    _expect_attrs(element, {}, label)
    _expect_children(element, [], label)
    return _text(element, label)


def _positive_decimal(text: str, label: str) -> Decimal:
    if not _DECIMAL_RE.fullmatch(text):
        raise IntakeError("source-noncanonical-number:%s" % label)
    try:
        value = Decimal(text)
    except InvalidOperation as exc:
        raise IntakeError("source-invalid-number:%s" % label) from exc
    if not value.is_finite() or value <= 0:
        raise IntakeError("source-nonpositive-number:%s" % label)
    return value


def _parse_dataset(payload: bytes, authority: Dict[str, Any]) -> Dict[str, Any]:
    if b"\r" in payload or not payload.endswith(b"\n"):
        raise IntakeError("source-not-canonical-lf:%d" % authority["record_id"])
    upper = payload.upper()
    if b"<!DOCTYPE" in upper or b"<!ENTITY" in upper:
        raise IntakeError("source-forbidden-xml-declaration:%d" % authority["record_id"])
    try:
        root = ET.fromstring(payload)
    except ET.ParseError as exc:
        raise IntakeError(
            "source-xml-invalid:%d:%s" % (authority["record_id"], exc)
        ) from exc

    label = "record-%d" % authority["record_id"]
    if _local_name(root) != "experiment":
        raise IntakeError("source-root-mismatch:%s" % label)
    _expect_attrs(root, {}, label)
    _expect_children(
        root,
        [
            "fileAuthor",
            "fileDOI",
            "ReSpecThVersion",
            "bibliographyLink",
            "experimentType",
            "apparatus",
            "commonProperties",
            "dataGroup",
            "ignitionType",
        ],
        label,
    )
    children = list(root)

    if _leaf_text(children[0], label + ".fileAuthor") != EXPECTED_CAMPAIGN["dataset_author"]:
        raise IntakeError("source-author-mismatch:%s" % label)
    if _leaf_text(children[1], label + ".fileDOI") != authority["record_doi"]:
        raise IntakeError("source-doi-mismatch:%s" % label)

    version = children[2]
    _expect_attrs(version, {}, label + ".ReSpecThVersion")
    _expect_children(version, ["major", "minor"], label + ".ReSpecThVersion")
    if [_leaf_text(item, label + ".version") for item in list(version)] != ["2", "3"]:
        raise IntakeError("source-version-mismatch:%s" % label)

    bibliography = children[3]
    _expect_attrs(bibliography, {}, label + ".bibliographyLink")
    _expect_children(
        bibliography, ["description", "referenceDOI"], label + ".bibliographyLink"
    )
    _leaf_text(list(bibliography)[0], label + ".description")
    if _leaf_text(list(bibliography)[1], label + ".referenceDOI") != PUBLICATION_DOI:
        raise IntakeError("source-reference-doi-mismatch:%s" % label)

    if _leaf_text(children[4], label + ".experimentType") != "ignition delay measurement":
        raise IntakeError("source-experiment-type-mismatch:%s" % label)
    apparatus = children[5]
    _expect_attrs(apparatus, {}, label + ".apparatus")
    _expect_children(apparatus, ["kind"], label + ".apparatus")
    if _leaf_text(list(apparatus)[0], label + ".apparatus.kind") != "shock tube":
        raise IntakeError("source-apparatus-mismatch:%s" % label)

    common = children[6]
    _expect_attrs(common, {}, label + ".commonProperties")
    _expect_children(common, ["property", "property"], label + ".commonProperties")
    pressure_node, composition_node = list(common)
    _expect_attrs(
        pressure_node,
        {"name": "pressure", "sourcetype": "reported", "units": "atm"},
        label + ".pressure",
    )
    _expect_children(pressure_node, ["value"], label + ".pressure")
    pressure_text = _leaf_text(list(pressure_node)[0], label + ".pressure.value")
    _positive_decimal(pressure_text, label + ".pressure.value")

    _expect_attrs(
        composition_node,
        {"name": "initial composition", "sourcetype": "reported"},
        label + ".composition",
    )
    composition: List[Dict[str, str]] = []
    seen_species = set()
    total = Decimal("0")
    for index, component in enumerate(list(composition_node)):
        component_label = "%s.composition[%d]" % (label, index)
        if _local_name(component) != "component":
            raise IntakeError("source-composition-child-mismatch:%s" % component_label)
        _expect_attrs(component, {}, component_label)
        _expect_children(component, ["speciesLink", "amount"], component_label)
        species_node, amount_node = list(component)
        if set(species_node.attrib) != {"preferredKey"} or list(species_node):
            raise IntakeError("source-species-contract-mismatch:%s" % component_label)
        species = species_node.attrib["preferredKey"]
        if not re.fullmatch(r"[A-Z][A-Z0-9]*", species) or species in seen_species:
            raise IntakeError("source-species-invalid:%s" % component_label)
        seen_species.add(species)
        _expect_attrs(amount_node, {"units": "mole fraction"}, component_label + ".amount")
        _expect_children(amount_node, [], component_label + ".amount")
        amount_text = _text(amount_node, component_label + ".amount")
        amount = _positive_decimal(amount_text, component_label + ".amount")
        total += amount
        composition.append({"species": species, "mole_fraction": amount_text})
    if len(composition) != 3 or total != Decimal("1.0000000000"):
        raise IntakeError("source-composition-closure-mismatch:%s" % label)

    data_group = children[7]
    _expect_attrs(data_group, {"id": "dg1", "label": "experimental_data"}, label + ".dataGroup")
    data_children = list(data_group)
    if len(data_children) < 3 or [_local_name(item) for item in data_children[:2]] != [
        "property",
        "property",
    ]:
        raise IntakeError("source-data-group-contract-mismatch:%s" % label)
    if any(_local_name(item) != "dataPoint" for item in data_children[2:]):
        raise IntakeError("source-data-point-contract-mismatch:%s" % label)

    property_by_id: Dict[str, str] = {}
    source_types: Dict[str, str] = {}
    for index, property_node in enumerate(data_children[:2]):
        property_label = "%s.dataProperty[%d]" % (label, index)
        attrs = property_node.attrib
        if set(attrs) != {"name", "id", "sourcetype", "units"}:
            raise IntakeError("source-data-property-attributes:%s" % property_label)
        name = attrs["name"]
        expected_units = {"temperature": "K", "ignition delay": "us"}.get(name)
        if expected_units is None or attrs["units"] != expected_units:
            raise IntakeError("source-data-property-units:%s" % property_label)
        if attrs["sourcetype"] != "digitized" or attrs["id"] in property_by_id:
            raise IntakeError("source-data-property-provenance:%s" % property_label)
        if list(property_node):
            raise IntakeError("source-data-property-not-empty:%s" % property_label)
        property_by_id[attrs["id"]] = name
        source_types[name] = attrs["sourcetype"]
    if set(property_by_id.values()) != {"temperature", "ignition delay"}:
        raise IntakeError("source-observable-set-mismatch:%s" % label)

    property_order = [item.attrib["id"] for item in data_children[:2]]
    points: List[Dict[str, str]] = []
    for index, point_node in enumerate(data_children[2:]):
        point_label = "%s.dataPoint[%d]" % (label, index)
        _expect_attrs(point_node, {}, point_label)
        _expect_children(point_node, property_order, point_label)
        values: Dict[str, str] = {}
        for value_node in list(point_node):
            value_label = point_label + "." + _local_name(value_node)
            value_text = _leaf_text(value_node, value_label)
            _positive_decimal(value_text, value_label)
            values[property_by_id[_local_name(value_node)]] = value_text
        points.append(
            {
                "temperature_k": values["temperature"],
                "ignition_delay_us": values["ignition delay"],
            }
        )

    ignition = children[8]
    _expect_attrs(ignition, {"target": "OH", "type": "d/dt max"}, label + ".ignitionType")
    _expect_children(ignition, [], label + ".ignitionType")

    for element in root.iter():
        if "uncert" in _local_name(element).lower():
            raise IntakeError("source-unexpected-uncertainty-element:%s" % label)
        if any("uncert" in key.lower() for key in element.attrib):
            raise IntakeError("source-unexpected-uncertainty-attribute:%s" % label)

    contains_c3h8 = "C3H8" in seen_species
    return {
        "record_id": authority["record_id"],
        "record_doi": authority["record_doi"],
        "pressure": {
            "value_atm": pressure_text,
            "source_type": "reported",
        },
        "initial_composition": composition,
        "composition_units": "mole fraction",
        "composition_source_type": "reported",
        "observables": {
            "temperature": {"units": "K", "source_type": source_types["temperature"]},
            "ignition_delay": {"units": "us", "source_type": source_types["ignition delay"]},
        },
        "ignition_definition": {"target": "OH", "method": "d/dt max"},
        "uncertainty_available": {
            "pressure": False,
            "temperature": False,
            "composition": False,
            "ignition_delay": False,
        },
        "contains_c3h8": contains_c3h8,
        "point_count": len(points),
        "points": points,
    }


def validate_and_extract(root: Path, manifest_path: Path) -> Dict[str, Any]:
    manifest = _decode_manifest(manifest_path)
    datasets: List[Dict[str, Any]] = []
    for authority in EXPECTED_SOURCES:
        label = str(authority["record_id"])
        source_path = _resolve_source(root, authority["repository_path"], label)
        try:
            payload = source_path.read_bytes()
        except OSError as exc:
            raise IntakeError("source-read-failed:%s:%s" % (label, exc)) from exc
        if len(payload) != authority["byte_count"]:
            raise IntakeError("source-byte-count-mismatch:%s" % label)
        if hashlib.sha256(payload).hexdigest() != authority["sha256"]:
            raise IntakeError("source-sha256-mismatch:%s" % label)
        if hashlib.md5(payload).hexdigest() != authority["zenodo_md5"]:
            raise IntakeError("source-zenodo-md5-mismatch:%s" % label)
        datasets.append(_parse_dataset(payload, authority))

    species = sorted(
        {
            item["species"]
            for dataset in datasets
            for item in dataset["initial_composition"]
        }
    )
    pressures: List[str] = []
    for dataset in datasets:
        pressure = dataset["pressure"]["value_atm"]
        if pressure not in pressures:
            pressures.append(pressure)
    summary = {
        "record_count": len(datasets),
        "curve_count": len(datasets),
        "point_count": sum(dataset["point_count"] for dataset in datasets),
        "species_set": species,
        "pressure_values_atm": pressures,
        "c3h8_curve_count": sum(bool(dataset["contains_c3h8"]) for dataset in datasets),
        "independent_campaign_count": 1,
    }
    observed = manifest["observed_structure"]
    comparisons = {
        "curve_count": summary["curve_count"],
        "point_count": summary["point_count"],
        "species_set": summary["species_set"],
        "pressure_values_atm": summary["pressure_values_atm"],
        "c3h8_curve_count": summary["c3h8_curve_count"],
    }
    for key, actual in comparisons.items():
        if actual != observed[key]:
            raise IntakeError("extraction-summary-mismatch:%s" % key)
    if any(dataset["contains_c3h8"] for dataset in datasets):
        raise IntakeError("eligibility-disposition-stale:c3h8-curve-observed")

    return {
        "schema": RESULT_SCHEMA,
        "status": "BLOCKED",
        "scope": manifest["scope"],
        "campaign": manifest["campaign"],
        "license": manifest["license"],
        "datasets": datasets,
        "summary": summary,
        "eligibility": manifest["eligibility"],
        "approvals": manifest["approvals"],
    }


def _canonical_result(result: Dict[str, Any]) -> bytes:
    return (json.dumps(result, indent=2, ensure_ascii=True) + "\n").encode("ascii")


def _write_atomic(path: Path, payload: bytes, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise IntakeError("output-exists-use-overwrite")
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise IntakeError("output-parent-create-failed:%s" % exc) from exc
    temp_name: Optional[str] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", prefix=path.name + ".", suffix=".tmp", dir=str(path.parent), delete=False
        ) as handle:
            temp_name = handle.name
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_name, str(path))
        temp_name = None
    except OSError as exc:
        raise IntakeError("output-write-failed:%s" % exc) from exc
    finally:
        if temp_name is not None:
            try:
                os.unlink(temp_name)
            except OSError:
                pass


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate and extract the bounded SciExpeM R2 intake screen."
    )
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        args = _parser().parse_args(argv)
        if args.overwrite and args.output is None:
            raise IntakeError("overwrite-requires-output")
        result = validate_and_extract(args.root, args.manifest)
        payload = _canonical_result(result)
        extraction = EXPECTED_EXTRACTION
        if len(payload) != extraction["expected_byte_count"]:
            raise IntakeError("extraction-byte-count-mismatch")
        if hashlib.sha256(payload).hexdigest() != extraction["expected_sha256"]:
            raise IntakeError("extraction-sha256-mismatch")
        if args.output is not None:
            _write_atomic(args.output, payload, args.overwrite)
        summary = result["summary"]
        print(
            "%s records=%d curves=%d points=%d c3h8_curves=%d "
            "training_eligible=0 holdout_eligible=0"
            % (
                BLOCKED_MARKER,
                summary["record_count"],
                summary["curve_count"],
                summary["point_count"],
                summary["c3h8_curve_count"],
            )
        )
        return 0
    except IntakeError as exc:
        print(
            json.dumps(
                {"schema": "rfam.r2-experimental-intake-check/v1", "status": "FAIL", "reason": str(exc)},
                sort_keys=True,
                separators=(",", ":"),
            ),
            file=sys.stderr,
        )
        return 2
    except Exception as exc:  # pragma: no cover - fail-closed safety boundary
        print(
            json.dumps(
                {
                    "schema": "rfam.r2-experimental-intake-check/v1",
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
