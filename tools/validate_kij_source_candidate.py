#!/usr/bin/env python3
"""Validate and optionally extract the O2/C3H8 ThermoML source candidate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import sys
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any


SCHEMA = "rfam.kij-source-candidate/v1"
CSV_SCHEMA = "rfam.kij-source-candidate-csv/v1"
BLOCKED_MARKER = "RFAM_KIJ_SOURCE_CANDIDATE_BLOCKED"
FAIL_MARKER = "RFAM_KIJ_SOURCE_CANDIDATE_FAIL"

SOURCE_URL = "https://trc.nist.gov/ThermoML/10.1021/je100574y.xml"
SOURCE_SHA256 = "e96a620f208243982427805dd2b4194a6a24a54b1e7764382b61e09134decc9a"
SOURCE_BYTE_COUNT = 26823
SOURCE_LAST_MODIFIED_HTTP = "Tue, 21 May 2024 23:20:23 GMT"
CANDIDATE_SHA256 = "5e6c646fda6e03c84db4aa30554c927b7e27fe3e4cd28dae7a2a524b1071ccf5"
THERMOML_NAMESPACE = "http://www.iupac.org/namespaces/ThermoML"
NS = {"t": THERMOML_NAMESPACE}

EXPECTED_APPROVALS = {
    "source_selected": False,
    "fit_approved": False,
    "kij_value_validated": False,
    "production_ready": False,
    "hem_enabled": False,
    "fluent_execution_approved": False,
    "hot_start_approved": False,
}

EXPECTED_UPSTREAM = [
    {
        "id": "wp3_sensitivity_contract",
        "path": "realFluidAllMach/tests/reference/kij_sensitivity_contract_v1.json",
        "sha256": "e1501157afdf1d2d51bb6fa0d03413b6da3f4d0cee4a0fbd9b8cfa9a39c99752",
    },
    {
        "id": "wp3_sensitivity_owner_report",
        "path": "tmp/a408-wp3-o2-c3h8-kij-sensitivity-screen-replacement-v2-20260722.md",
        "sha256": "142d01fb4cf057885de56c521fb8e504dfa1ae48b68daa81641eb8a85795e851",
    },
    {
        "id": "wp3_sensitivity_independent_review",
        "path": "tmp/a408-wp3-o2-c3h8-kij-sensitivity-screen-replacement-v2-independent-review-20260722.md",
        "sha256": "f02c945863f51f69ee6d00f81deefed5afa02b789a372c9283bc67723ddbc1e3",
    },
]

EXPECTED_PAIR = {
    "species": ["O2", "C3H8"],
    "source_composition_basis": "liquid-phase mole fraction of O2",
    "fit_model_target": "RFAM Peng-Robinson symmetric kij",
    "fit_functional_form_selected": False,
}

EXPECTED_CITATION = {
    "doi": "10.1021/je100574y",
    "title": (
        "Equilibrium Data for the Oxygen + Propane Binary System at "
        "Temperatures of (110.22, 120.13, 130.58, and 139.95) K"
    ),
    "authors": [
        "Houssin-Agbomson, D.[Deborah]",
        "Coquelet, C.[Christophe]",
        "Arpentinier, P.",
        "Delcorso, F.[Fabrice]",
        "Richon, D.[Dominique]",
    ],
    "journal": "J. Chem. Eng. Data",
    "year": 2010,
    "volume": "55",
    "pages": "4412-4415",
}

EXPECTED_SOURCE_ARTIFACT = {
    "provider": "NIST Thermodynamics Research Center",
    "archive": "ThermoML Archive",
    "url": SOURCE_URL,
    "media_type": "application/xml",
    "thermoml_version": [2, 0],
    "retrieved_on": "2026-07-22",
    "last_modified_utc": "2024-05-21T23:20:23Z",
    "sha256": SOURCE_SHA256,
    "byte_count": SOURCE_BYTE_COUNT,
    "repository_copy": False,
}

EXPECTED_RIGHTS = {
    "classification": "unresolved-artifact-specific-redistribution-and-fit-use",
    "source_xml_license_marker": "none-observed",
    "source_xml_copyright_marker": "none-observed",
    "evidence": [
        {
            "id": "nist-general-use",
            "url": "https://www.nist.gov/copyrights-disclaimers",
            "accessed_on": "2026-07-22",
            "observation": (
                "NIST says unmarked public information may be distributed or "
                "copied, but product-specific terms can apply."
            ),
        },
        {
            "id": "nist-srd-license-policy",
            "url": "https://www.nist.gov/open/license",
            "accessed_on": "2026-07-22",
            "observation": (
                "NIST states that Standard Reference Data compilations can be "
                "copyrighted and licensed."
            ),
        },
        {
            "id": "trc-public-availability",
            "url": "https://trc.nist.gov/",
            "accessed_on": "2026-07-22",
            "observation": (
                "TRC states that journal-submitted data are publicly available "
                "through the ThermoML Archive."
            ),
        },
        {
            "id": "crossref-license-metadata",
            "url": "https://api.crossref.org/works/10.1021%2Fje100574y",
            "accessed_on": "2026-07-22",
            "observation": "The Crossref record contains no license entry.",
        },
        {
            "id": "openalex-access-metadata",
            "url": "https://api.openalex.org/works/doi:10.1021/je100574y",
            "accessed_on": "2026-07-22",
            "observation": (
                "OpenAlex classifies the article as closed access and reports "
                "no best open-access location."
            ),
        },
        {
            "id": "hal-access-metadata",
            "url": (
                "https://api.archives-ouvertes.fr/search/"
                "?q=halId_s%3Ahal-00573772"
            ),
            "accessed_on": "2026-07-22",
            "observation": (
                "HAL record hal-00573772 reports openAccess=false and exposes "
                "metadata without an article file."
            ),
        },
    ],
    "contract_disposition": (
        "block-fit-and-repository-redistribution-pending-resolution"
    ),
    "required_before_selection": [
        "determine whether this ThermoML artifact is governed as SRD or non-SRD data",
        "record explicit terms permitting the intended numerical use and redistribution",
        "retain attribution and any required modification notice in the selected-data manifest",
    ],
}

EXPECTED_EXPERIMENT = {
    "method": "static-analytic method coupled to gas chromatograph analysis",
    "method_source": "ThermoML citation abstract",
    "components": [
        {
            "name": "oxygen",
            "formula": "O2",
            "inchi": "InChI=1S/O2/c1-2",
            "inchi_key": "MYMOFIZGZYHOMD-UHFFFAOYSA-N",
            "supplier_purity_mol_percent": 99.9995,
        },
        {
            "name": "propane",
            "formula": "C3H8",
            "inchi": "InChI=1S/C3H8/c1-3-2/h3H2,1-2H3",
            "inchi_key": "ATUOYWHBWRKTHZ-UHFFFAOYSA-N",
            "supplier_purity_mol_percent": 99.95,
        },
    ],
    "temperature_k": [110.22, 120.13, 130.58, 139.95],
    "datasets": [
        {
            "dataset_id": 1,
            "liquid_phase_label": "Liquid mixture 2",
            "branch": "oxygen-rich",
            "row_count": 8,
            "liquid_x_o2_range": [0.9115, 0.9776],
            "pressure_kpa_range": [535.9, 2611.6],
            "expanded_pressure_uncertainty_kpa_range": [3.8, 16.7],
        },
        {
            "dataset_id": 2,
            "liquid_phase_label": "Liquid mixture 1",
            "branch": "propane-rich",
            "row_count": 20,
            "liquid_x_o2_range": [0.0266, 0.419],
            "pressure_kpa_range": [246.6, 2435.3],
            "expanded_pressure_uncertainty_kpa_range": [10.8, 156.7],
        },
    ],
    "total_row_count": 28,
    "reported_variables": ["temperature_k", "liquid_x_o2"],
    "reported_property": "vapor or sublimation pressure, kPa",
    "pressure_uncertainty": {
        "kind": "expanded",
        "confidence_percent": 95,
        "evaluator": "Data file compiler",
        "method": "Propagation of evaluated standard uncertainties",
        "present_for_every_row": True,
    },
    "missing_measurement_fields": [
        "temperature_uncertainty_k",
        "liquid_x_o2_uncertainty",
        "vapor_y_o2",
    ],
    "abstract_reports_vlle_visual_confirmation": True,
    "xml_contains_explicit_vlle_tie_lines": False,
}

CANDIDATE_FIELDS = (
    "dataset_id",
    "liquid_phase",
    "branch",
    "temperature_k",
    "temperature_uncertainty_k",
    "liquid_x_o2",
    "liquid_x_o2_uncertainty",
    "vapor_y_o2",
    "pressure_kpa",
    "expanded_pressure_uncertainty_kpa",
    "pressure_uncertainty_confidence_percent",
)

EXPECTED_CANDIDATE_EXTRACTION = {
    "schema": CSV_SCHEMA,
    "fields": list(CANDIDATE_FIELDS),
    "row_count": 28,
    "canonicalization": (
        "ASCII CSV with LF endings, source lexical numbers, missing fields empty"
    ),
    "sha256": CANDIDATE_SHA256,
    "repository_copy": False,
    "fit_input": False,
}

EXPECTED_ELIGIBILITY = {
    "source_role": "candidate-pressure-liquid-composition-screening-only",
    "training_eligible": False,
    "holdout_eligible": False,
    "fit_performed": False,
    "validated_kij": False,
    "blocking_reasons": [
        "artifact-specific numerical-use and redistribution terms are unresolved",
        "temperature measurement uncertainty is absent from the ThermoML artifact",
        "liquid-composition measurement uncertainty is absent from the ThermoML artifact",
        (
            "vapor-phase composition is absent, so dew and phase-composition "
            "validation are unavailable"
        ),
        "no independent holdout dataset is selected",
    ],
    "required_next": [
        (
            "obtain an authorized article or author/source record for "
            "measurement accuracies and terms"
        ),
        "identify an independent O2-C3H8 VLE or PVT holdout source",
        (
            "freeze constant-versus-temperature-dependent kij forms, "
            "partition, metrics, and tolerances before fitting"
        ),
    ],
}

EXPECTED_LIMITATIONS = [
    "public URL accessibility is not treated as an artifact-specific license",
    "the 28 P-x rows do not contain vapor composition or density data",
    "the candidate CSV is a temporary source audit artifact and is not a fit dataset",
    "no value including historical kij=0.055 is validated or recommended",
    (
        "no HEM, R2 physics, soot, Fluent, zero-advance, hot-start, or "
        "production capability is enabled"
    ),
]

EXPECTED_ROOT_KEYS = (
    "schema",
    "contract_version",
    "status",
    "scope",
    "approvals",
    "upstream",
    "pair",
    "citation",
    "source_artifact",
    "rights",
    "experiment",
    "candidate_extraction",
    "eligibility",
    "limitations",
)


class SourceCandidateError(ValueError):
    """A fail-closed manifest, source, or extraction violation."""


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise SourceCandidateError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def _decode_manifest(path: Path) -> dict[str, Any]:
    try:
        payload = path.read_bytes()
    except OSError as exc:
        raise SourceCandidateError(f"manifest-read-failed:{exc}") from exc
    if b"\r" in payload or not payload.endswith(b"\n"):
        raise SourceCandidateError("manifest-not-canonical-lf")
    try:
        text = payload.decode("ascii")
        document = json.loads(text, object_pairs_hook=_reject_duplicate_pairs)
    except SourceCandidateError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SourceCandidateError(f"manifest-json-invalid:{exc}") from exc
    if not isinstance(document, dict):
        raise SourceCandidateError("manifest-root-not-object")
    canonical = json.dumps(document, indent=2, ensure_ascii=True) + "\n"
    if text != canonical:
        raise SourceCandidateError("manifest-not-canonical-json")
    return document


def _strict_equal(actual: Any, expected: Any, label: str) -> None:
    if type(actual) is not type(expected):
        raise SourceCandidateError(f"contract-type-mismatch:{label}")
    if isinstance(expected, dict):
        actual_keys = list(actual)
        expected_keys = list(expected)
        if actual_keys != expected_keys:
            missing = [key for key in expected_keys if key not in actual]
            unknown = [key for key in actual_keys if key not in expected]
            if missing:
                raise SourceCandidateError(
                    f"contract-missing-keys:{label}:" + ",".join(missing)
                )
            if unknown:
                raise SourceCandidateError(
                    f"contract-unknown-keys:{label}:" + ",".join(unknown)
                )
            raise SourceCandidateError(f"contract-key-order-mismatch:{label}")
        for key in expected_keys:
            _strict_equal(actual[key], expected[key], f"{label}.{key}")
        return
    if isinstance(expected, list):
        if len(actual) != len(expected):
            raise SourceCandidateError(f"contract-list-length-mismatch:{label}")
        for index, expected_item in enumerate(expected):
            _strict_equal(actual[index], expected_item, f"{label}[{index}]")
        return
    if actual != expected:
        raise SourceCandidateError(f"contract-value-mismatch:{label}")


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _workspace_path(root: Path, relative: str, label: str) -> Path:
    if not isinstance(relative, str) or not relative or "\\" in relative:
        raise SourceCandidateError(f"invalid-path:{label}")
    if relative.startswith("/") or re.match(r"^[A-Za-z]:", relative):
        raise SourceCandidateError(f"path-not-relative:{label}")
    resolved_root = root.resolve()
    candidate = (resolved_root / relative).resolve()
    if candidate != resolved_root and resolved_root not in candidate.parents:
        raise SourceCandidateError(f"path-escapes-root:{label}")
    return candidate


def validate_manifest(root: Path, manifest_path: Path) -> dict[str, Any]:
    document = _decode_manifest(manifest_path)
    if list(document) != list(EXPECTED_ROOT_KEYS):
        missing = [key for key in EXPECTED_ROOT_KEYS if key not in document]
        unknown = [key for key in document if key not in EXPECTED_ROOT_KEYS]
        if missing:
            raise SourceCandidateError("contract-missing-keys:root:" + ",".join(missing))
        if unknown:
            raise SourceCandidateError("contract-unknown-keys:root:" + ",".join(unknown))
        raise SourceCandidateError("contract-key-order-mismatch:root")

    _strict_equal(document["schema"], SCHEMA, "root.schema")
    _strict_equal(document["contract_version"], 1, "root.contract_version")
    _strict_equal(document["status"], "blocked-source-candidate", "root.status")
    _strict_equal(
        document["scope"],
        "offline-o2-c3h8-vle-source-screen-no-fit",
        "root.scope",
    )
    _strict_equal(document["approvals"], EXPECTED_APPROVALS, "root.approvals")
    _strict_equal(document["upstream"], EXPECTED_UPSTREAM, "root.upstream")
    _strict_equal(document["pair"], EXPECTED_PAIR, "root.pair")
    _strict_equal(document["citation"], EXPECTED_CITATION, "root.citation")
    _strict_equal(
        document["source_artifact"], EXPECTED_SOURCE_ARTIFACT, "root.source_artifact"
    )
    _strict_equal(document["rights"], EXPECTED_RIGHTS, "root.rights")
    _strict_equal(document["experiment"], EXPECTED_EXPERIMENT, "root.experiment")
    _strict_equal(
        document["candidate_extraction"],
        EXPECTED_CANDIDATE_EXTRACTION,
        "root.candidate_extraction",
    )
    _strict_equal(document["eligibility"], EXPECTED_ELIGIBILITY, "root.eligibility")
    _strict_equal(document["limitations"], EXPECTED_LIMITATIONS, "root.limitations")

    for item in EXPECTED_UPSTREAM:
        path = _workspace_path(root, item["path"], f"upstream:{item['id']}")
        try:
            payload = path.read_bytes()
        except OSError as exc:
            raise SourceCandidateError(
                f"upstream-read-failed:{item['id']}:{exc}"
            ) from exc
        if _sha256(payload) != item["sha256"]:
            raise SourceCandidateError(f"upstream-sha256-mismatch:{item['id']}")
    return document


def _local_name(element: ET.Element) -> str:
    return element.tag.rsplit("}", 1)[-1]


def _text(parent: ET.Element, path: str, label: str) -> str:
    element = parent.find(path, NS)
    if element is None or element.text is None or not element.text.strip():
        raise SourceCandidateError(f"source-missing-text:{label}")
    return element.text.strip()


def _children(element: ET.Element) -> list[str]:
    return [_local_name(child) for child in element]


def _expect_children(element: ET.Element, expected: list[str], label: str) -> None:
    if _children(element) != expected:
        raise SourceCandidateError(f"source-child-contract-mismatch:{label}")


def _finite(text: str, label: str) -> float:
    try:
        value = float(text)
    except ValueError as exc:
        raise SourceCandidateError(f"source-invalid-number:{label}") from exc
    if not math.isfinite(value):
        raise SourceCandidateError(f"source-nonfinite-number:{label}")
    return value


def _validate_citation(root: ET.Element) -> None:
    citation = root.find("t:Citation", NS)
    if citation is None:
        raise SourceCandidateError("source-citation-missing")
    observed = {
        "doi": _text(citation, "t:sDOI", "citation.doi"),
        "title": _text(citation, "t:sTitle", "citation.title"),
        "authors": [
            (element.text or "").strip()
            for element in citation.findall("t:sAuthor", NS)
        ],
        "journal": _text(citation, "t:sPubName", "citation.journal"),
        "year": int(_text(citation, "t:yrPubYr", "citation.year")),
        "volume": _text(citation, "t:sVol", "citation.volume"),
        "pages": _text(citation, "t:sPage", "citation.pages"),
    }
    _strict_equal(observed, EXPECTED_CITATION, "source.citation")
    abstract = _text(citation, "t:sAbstract", "citation.abstract")
    required_fragments = (
        "static-analytic",
        "gas chromatograph",
        "vapor-liquid-liquid equilibrium",
        "confirmed by visual observations",
    )
    if any(fragment not in abstract for fragment in required_fragments):
        raise SourceCandidateError("source-abstract-method-contract-mismatch")


def _validate_compounds(root: ET.Element) -> None:
    expected = [
        {
            "number": "1",
            "name": "oxygen",
            "formula": "O2",
            "inchi": "InChI=1S/O2/c1-2",
            "key": "MYMOFIZGZYHOMD-UHFFFAOYSA-N",
            "purity": "99.9995",
        },
        {
            "number": "2",
            "name": "propane",
            "formula": "C3H8",
            "inchi": "InChI=1S/C3H8/c1-3-2/h3H2,1-2H3",
            "key": "ATUOYWHBWRKTHZ-UHFFFAOYSA-N",
            "purity": "99.95",
        },
    ]
    compounds = root.findall("t:Compound", NS)
    if len(compounds) != len(expected):
        raise SourceCandidateError("source-compound-count-mismatch")
    for index, (compound, item) in enumerate(zip(compounds, expected)):
        observed = {
            "number": _text(compound, "t:RegNum/t:nOrgNum", f"compound[{index}].number"),
            "name": _text(compound, "t:sCommonName", f"compound[{index}].name"),
            "formula": _text(compound, "t:sFormulaMolec", f"compound[{index}].formula"),
            "inchi": _text(compound, "t:sStandardInChI", f"compound[{index}].inchi"),
            "key": _text(compound, "t:sStandardInChIKey", f"compound[{index}].key"),
            "purity": _text(
                compound, "t:Sample/t:purity/t:nPurityMol", f"compound[{index}].purity"
            ),
        }
        _strict_equal(observed, item, f"source.compound[{index}]")


def _validate_dataset_definition(
    data: ET.Element, dataset_id: str, liquid_phase: str, label: str
) -> str:
    components = [
        _text(component, "t:RegNum/t:nOrgNum", f"{label}.component")
        for component in data.findall("t:Component", NS)
    ]
    _strict_equal(components, ["2", "1"], f"source.{label}.components")
    _strict_equal(
        _text(data, "t:nPureOrMixtureDataNumber", f"{label}.id"),
        dataset_id,
        f"source.{label}.id",
    )

    prop = data.find("t:Property", NS)
    if prop is None:
        raise SourceCandidateError(f"source-property-missing:{label}")
    property_contract = {
        "number": _text(prop, "t:nPropNumber", f"{label}.property.number"),
        "name": _text(
            prop,
            "t:Property-MethodID/t:PropertyGroup/t:VaporPBoilingTAzeotropTandP/t:ePropName",
            f"{label}.property.name",
        ),
        "method": _text(
            prop,
            "t:Property-MethodID/t:PropertyGroup/t:VaporPBoilingTAzeotropTandP/t:sMethodName",
            f"{label}.property.method",
        ),
        "phase": _text(prop, "t:PropPhaseID/t:ePropPhase", f"{label}.property.phase"),
        "presentation": _text(prop, "t:ePresentation", f"{label}.property.presentation"),
        "evaluator": _text(
            prop, "t:CombinedUncertainty/t:sCombUncertEvaluator", f"{label}.uncertainty.evaluator"
        ),
        "uncertainty_method": _text(
            prop,
            "t:CombinedUncertainty/t:eCombUncertEvalMethod",
            f"{label}.uncertainty.method",
        ),
        "confidence": _text(
            prop,
            "t:CombinedUncertainty/t:nCombUncertLevOfConfid",
            f"{label}.uncertainty.confidence",
        ),
    }
    expected_property = {
        "number": "1",
        "name": "Vapor or sublimation pressure, kPa",
        "method": "CCELL:UFactor:2",
        "phase": liquid_phase,
        "presentation": "Direct value, X",
        "evaluator": "Data file compiler",
        "uncertainty_method": "Propagation of evaluated standard uncertainties",
        "confidence": "95",
    }
    _strict_equal(property_contract, expected_property, f"source.{label}.property")

    phases = [
        _text(phase, "t:ePhase", f"{label}.phase")
        for phase in data.findall("t:PhaseID", NS)
    ]
    _strict_equal(phases, [liquid_phase, "Gas"], f"source.{label}.phases")
    variables = data.findall("t:Variable", NS)
    if len(variables) != 2:
        raise SourceCandidateError(f"source-variable-count-mismatch:{label}")
    first, second = variables
    first_contract = {
        "number": _text(first, "t:nVarNumber", f"{label}.variable[0].number"),
        "type": _text(
            first, "t:VariableID/t:VariableType/t:eTemperature", f"{label}.variable[0].type"
        ),
        "phase": _text(first, "t:VarPhaseID/t:eVarPhase", f"{label}.variable[0].phase"),
    }
    second_contract = {
        "number": _text(second, "t:nVarNumber", f"{label}.variable[1].number"),
        "type": _text(
            second,
            "t:VariableID/t:VariableType/t:eComponentComposition",
            f"{label}.variable[1].type",
        ),
        "component": _text(
            second, "t:VariableID/t:RegNum/t:nOrgNum", f"{label}.variable[1].component"
        ),
        "phase": _text(second, "t:VarPhaseID/t:eVarPhase", f"{label}.variable[1].phase"),
    }
    _strict_equal(
        first_contract,
        {"number": "1", "type": "Temperature, K", "phase": liquid_phase},
        f"source.{label}.variable[0]",
    )
    _strict_equal(
        second_contract,
        {
            "number": "2",
            "type": "Mole fraction",
            "component": "1",
            "phase": liquid_phase,
        },
        f"source.{label}.variable[1]",
    )
    return property_contract["confidence"]


def _candidate_csv(root: ET.Element) -> bytes:
    datasets = root.findall("t:PureOrMixtureData", NS)
    if len(datasets) != 2:
        raise SourceCandidateError("source-dataset-count-mismatch")
    definitions = (
        ("1", "Liquid mixture 2", "oxygen-rich", 8),
        ("2", "Liquid mixture 1", "propane-rich", 20),
    )
    stream = io.StringIO(newline="")
    stream.write(f"# {CSV_SCHEMA}\n")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(CANDIDATE_FIELDS)
    row_count = 0
    for dataset_index, (data, definition) in enumerate(zip(datasets, definitions)):
        dataset_id, liquid_phase, branch, expected_rows = definition
        label = f"dataset[{dataset_index}]"
        confidence = _validate_dataset_definition(data, dataset_id, liquid_phase, label)
        rows = data.findall("t:NumValues", NS)
        if len(rows) != expected_rows:
            raise SourceCandidateError(f"source-row-count-mismatch:{label}")
        for row_index, row in enumerate(rows):
            row_label = f"{label}.row[{row_index}]"
            variables = row.findall("t:VariableValue", NS)
            if len(variables) != 2:
                raise SourceCandidateError(f"source-row-variable-count-mismatch:{row_label}")
            values: dict[str, str] = {}
            for variable_index, variable in enumerate(variables):
                _expect_children(
                    variable,
                    ["nVarNumber", "nVarValue", "nVarDigits"],
                    f"{row_label}.variable[{variable_index}]",
                )
                number = _text(variable, "t:nVarNumber", f"{row_label}.variable.number")
                if number in values:
                    raise SourceCandidateError(f"source-duplicate-variable:{row_label}:{number}")
                values[number] = _text(
                    variable, "t:nVarValue", f"{row_label}.variable.value"
                )
            if list(values) != ["1", "2"]:
                raise SourceCandidateError(f"source-row-variable-order-mismatch:{row_label}")

            properties = row.findall("t:PropertyValue", NS)
            if len(properties) != 1:
                raise SourceCandidateError(f"source-row-property-count-mismatch:{row_label}")
            prop = properties[0]
            _expect_children(
                prop,
                ["nPropNumber", "nPropValue", "nPropDigits", "CombinedUncertainty"],
                f"{row_label}.property",
            )
            if _text(prop, "t:nPropNumber", f"{row_label}.property.number") != "1":
                raise SourceCandidateError(f"source-row-property-number-mismatch:{row_label}")
            uncertainty = prop.find("t:CombinedUncertainty", NS)
            if uncertainty is None:
                raise SourceCandidateError(f"source-row-uncertainty-missing:{row_label}")
            _expect_children(
                uncertainty,
                ["nCombUncertAssessNum", "nCombExpandUncertValue"],
                f"{row_label}.uncertainty",
            )
            if _text(
                uncertainty, "t:nCombUncertAssessNum", f"{row_label}.uncertainty.number"
            ) != "1":
                raise SourceCandidateError(f"source-row-uncertainty-number-mismatch:{row_label}")
            temperature = values["1"]
            composition = values["2"]
            pressure = _text(prop, "t:nPropValue", f"{row_label}.pressure")
            pressure_uncertainty = _text(
                uncertainty,
                "t:nCombExpandUncertValue",
                f"{row_label}.pressure_uncertainty",
            )
            for name, text in (
                ("temperature", temperature),
                ("composition", composition),
                ("pressure", pressure),
                ("pressure_uncertainty", pressure_uncertainty),
            ):
                _finite(text, f"{row_label}.{name}")
            writer.writerow(
                (
                    dataset_id,
                    liquid_phase,
                    branch,
                    temperature,
                    "",
                    composition,
                    "",
                    "",
                    pressure,
                    pressure_uncertainty,
                    confidence,
                )
            )
            row_count += 1
    if row_count != 28:
        raise SourceCandidateError("source-total-row-count-mismatch")
    stream.write(f"# {BLOCKED_MARKER} rows=28 fit_approved=0\n")
    try:
        payload = stream.getvalue().encode("ascii")
    except UnicodeEncodeError as exc:
        raise SourceCandidateError("candidate-csv-not-ascii") from exc
    if _sha256(payload) != CANDIDATE_SHA256:
        raise SourceCandidateError("candidate-csv-sha256-mismatch")
    return payload


def validate_source(payload: bytes) -> bytes:
    if len(payload) != SOURCE_BYTE_COUNT:
        raise SourceCandidateError("source-byte-count-mismatch")
    if _sha256(payload) != SOURCE_SHA256:
        raise SourceCandidateError("source-sha256-mismatch")
    lowered = payload.lower()
    if b"<!doctype" in lowered or b"<!entity" in lowered:
        raise SourceCandidateError("source-xml-doctype-or-entity-forbidden")
    if any(marker in lowered for marker in (b"<license", b"<copyright", b"<rights")):
        raise SourceCandidateError("source-rights-marker-contract-mismatch")
    try:
        root = ET.fromstring(payload)
    except ET.ParseError as exc:
        raise SourceCandidateError(f"source-xml-invalid:{exc}") from exc
    if root.tag != f"{{{THERMOML_NAMESPACE}}}DataReport":
        raise SourceCandidateError("source-root-mismatch")
    version = [
        int(_text(root, "t:Version/t:nVersionMajor", "version.major")),
        int(_text(root, "t:Version/t:nVersionMinor", "version.minor")),
    ]
    _strict_equal(version, [2, 0], "source.version")
    _validate_citation(root)
    _validate_compounds(root)
    return _candidate_csv(root)


def fetch_source() -> bytes:
    request = urllib.request.Request(
        SOURCE_URL,
        headers={
            "Accept": "application/xml",
            "User-Agent": "Vortex-WP3-source-audit/1.0",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            if response.geturl() != SOURCE_URL:
                raise SourceCandidateError("source-final-url-mismatch")
            media_type = response.headers.get_content_type()
            if media_type != "application/xml":
                raise SourceCandidateError("source-media-type-mismatch")
            if response.headers.get("Last-Modified") != SOURCE_LAST_MODIFIED_HTTP:
                raise SourceCandidateError("source-last-modified-mismatch")
            payload = response.read(SOURCE_BYTE_COUNT + 1)
    except SourceCandidateError:
        raise
    except (OSError, urllib.error.URLError) as exc:
        raise SourceCandidateError(f"source-fetch-failed:{exc}") from exc
    return payload


def _write_candidate(path: Path, payload: bytes, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise SourceCandidateError("candidate-output-exists")
    if not path.parent.is_dir():
        raise SourceCandidateError("candidate-output-parent-missing")
    try:
        path.write_bytes(payload)
    except OSError as exc:
        raise SourceCandidateError(f"candidate-output-write-failed:{exc}") from exc


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[2]
    default_manifest = (
        Path(__file__).resolve().parents[1]
        / "tests/reference/kij_o2_c3h8_source_candidate_v1.json"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument("--manifest", type=Path, default=default_manifest)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--source-xml", type=Path)
    source.add_argument("--fetch-source", action="store_true")
    parser.add_argument("--candidate-output", type=Path)
    parser.add_argument("--overwrite-candidate-output", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        root = args.root.resolve()
        manifest = args.manifest
        if not manifest.is_absolute():
            manifest = (Path.cwd() / manifest).resolve()
        validate_manifest(root, manifest)
        if args.candidate_output is not None and not (
            args.source_xml is not None or args.fetch_source
        ):
            raise SourceCandidateError("candidate-output-requires-source-verification")
        if args.overwrite_candidate_output and args.candidate_output is None:
            raise SourceCandidateError("overwrite-requires-candidate-output")

        candidate: bytes | None = None
        if args.source_xml is not None:
            try:
                source_payload = args.source_xml.read_bytes()
            except OSError as exc:
                raise SourceCandidateError(f"source-read-failed:{exc}") from exc
            candidate = validate_source(source_payload)
        elif args.fetch_source:
            candidate = validate_source(fetch_source())

        if args.candidate_output is not None:
            assert candidate is not None
            _write_candidate(
                args.candidate_output.resolve(),
                candidate,
                args.overwrite_candidate_output,
            )
    except Exception as exc:  # Keep malformed or unexpected input fail-closed.
        reason = (
            str(exc)
            if isinstance(exc, SourceCandidateError)
            else f"unexpected:{type(exc).__name__}:{exc}"
        )
        reason = reason.replace("\r", " ").replace("\n", " ")
        print(f"{FAIL_MARKER} reason={reason}")
        return 2

    live = 1 if candidate is not None else 0
    print(
        f"{BLOCKED_MARKER} upstream=3 rows=28 live_source_verified={live} "
        f"candidate_sha256={CANDIDATE_SHA256} source_selected=0 "
        "fit_approved=0 production_ready=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
