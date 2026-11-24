#!/usr/bin/env python3
"""Validate the provenance-only WP2 detailed-mechanism candidate.

The only successful state is a blocked, internally consistent provenance
snapshot. This tool does not fit R2 parameters or authorize a physical run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path
from typing import Any


SCHEMA = "rfam.r2-physics-candidate/v1"
EXPECTED_ROOT_KEYS = {
    "schema", "contract_version", "status", "scope", "approvals", "upstream",
    "artifacts", "mechanism_candidate", "conversion", "runtime_probe",
    "experimental_evidence", "validation_protocol", "cross_scope", "limitations",
}

WORKSPACE_ARTIFACTS = {
    "ucsd_mechanism": {
        "path": "Vortex/tools/openfoam/sandiego20161214_mechCK.txt",
        "url": "https://web.eng.ucsd.edu/mae/groups/combustion/sdmech/sandiego20161214/sandiego20161214_mechCK.txt",
        "sha256": "728aa95e372aa71c59528850ac35511bf4c5f457954b9803c6ef823996163c73",
        "kind": "official-chemkin-mechanism",
    },
    "ucsd_thermo": {
        "path": "Vortex/tools/openfoam/sandiego20160815_therm.txt",
        "url": "https://web.eng.ucsd.edu/mae/groups/combustion/sdmech/sandiego20161214/sandiego20160815_therm.txt",
        "sha256": "ba9518b7007d007353ca6ac927f881b05d6b302bbb56bc32c9e1ee80d3333d92",
        "kind": "official-chemkin-thermo",
    },
    "ucsd_transport": {
        "path": "Vortex/tools/openfoam/sandiego20160815_trans.txt",
        "url": "https://web.eng.ucsd.edu/mae/groups/combustion/sdmech/sandiego20161214/sandiego20160815_trans.txt",
        "sha256": "dc8134f2c62cd54c7999d2428ac1de22743e2e331caef2da8286763bedd5ce94",
        "kind": "official-chemkin-transport",
    },
    "converted_yaml": {
        "path": "Vortex/tools/openfoam/sandiego2016_c3h8.yaml",
        "url": "",
        "sha256": "fbdbf6ce0bdbe85448229395f42096b59eb08828e92416d0f026bb6ae5040eaa",
        "kind": "cantera-converted-mechanism",
    },
    "cantera_wheel": {
        "path": "tmp/cantera-3.2.0-cp312-cp312-win_amd64.whl",
        "url": "https://files.pythonhosted.org/packages/b6/f7/9bd49339a73a6f45849ff860fdbf29004ea63a5ce693a89e81d75204076b/cantera-3.2.0-cp312-cp312-win_amd64.whl",
        "sha256": "258b95a4e29f9afd50788eab284e88440174302123c8b3a4a3a493904b40e238",
        "kind": "pypi-runtime-wheel",
    },
    "regenerated_yaml": {
        "path": "tmp/sandiego2016-regenerated-ct320.yaml",
        "url": "",
        "sha256": "dfc4abe7ff40e56279c0a143351cf98021dacb7e35f57727304a6de203b18908",
        "kind": "regeneration-output",
    },
}

UPSTREAM_ARTIFACTS = {
    "wp0_contract": {
        "path": "realFluidAllMach/tests/reference/vortex_physics_validation_contract_v1.json",
        "sha256": "34b15eae9c89d3c44638748e3283fed59429912f16faa33a224f70f78710439a",
    },
    "wp0_owner_report": {
        "path": "tmp/a408-wp0-wp2-replacement-v4-20260722.md",
        "sha256": "05ddf8a971cc79ca76196505ead2be37258cf7797a5f404eb14ec72fd5067420",
    },
    "wp0_independent_review": {
        "path": "tmp/a408-wp0-wp2-replacement-v4-independent-review-20260722.md",
        "sha256": "b05980810f2db2eaddfa48839479e9cf5cf2abd80940f700e9f1db986541d626",
    },
}

RUNTIME_ARTIFACTS = {
    "cantera320_python": {
        "path": "C:/Users/guhai/AppData/Local/Programs/Python/Python312/python.exe",
        "sha256": "4d6f5f81a4bca11191c4c7c6b43632694d0a4ce74e068619d8fdc161d469859a",
    },
    "cantera320_ck2yaml": {
        "path": "F:/tmp/cantera-3.2.0-python312/cantera/ck2yaml.py",
        "sha256": "30f95b1a0ecd7b6e7dc7f102c53fcec5a246e48472054ca3912140f31085017b",
    },
    "cantera320_extension": {
        "path": "F:/tmp/cantera-3.2.0-python312/cantera/_cantera.cp312-win_amd64.pyd",
        "sha256": "00f00749cc90b76a7266656fd02e6a7d40be82e7d7e1389472ad285123064e9d",
    },
    "cantera320_metadata": {
        "path": "F:/tmp/cantera-3.2.0-python312/cantera-3.2.0.dist-info/METADATA",
        "sha256": "f483c140b917b5e758f57ac8232ece2aeb7756325801395428470a64467824c7",
    },
    "cantera320_record": {
        "path": "F:/tmp/cantera-3.2.0-python312/cantera-3.2.0.dist-info/RECORD",
        "sha256": "8c9ccfd6d02ce4e423ea37267618b6475eb963b58b9bfd96fb39a7f7a8f7d630",
    },
}

EXPECTED_SOURCE_PAGE = "https://web.eng.ucsd.edu/mae/groups/combustion/mechanism.html"
EXPECTED_UPDATE_LOG_URL = "https://web.eng.ucsd.edu/mae/groups/combustion/sdmech/sandiego20161214/20161214_log.txt"
EXPECTED_UPDATE_LOG_SHA = "c15a1deda4a1c6938ba9a34ab3f05975bba883299ed69b4a66fd22cf4b3a1672"
EXPECTED_CANONICAL_YAML_SHA = "9885d44ff8ebcafa5fb52b261f994f46aecc85920bb2e058def71c859ba7e80d"
EXPECTED_RUNTIME_SPECIES_SHA = "a20ecbeb194b8fad036878944dc03519fdedf5d60822faf726010eb1ae7f9d6f"
EXPECTED_RUNTIME_EQUATION_SHA = "15f53906e58986bbdccce32864928315124af067b34a30a55cee714069fc94bb"
EXPECTED_RUNTIME_WARNING_SHA = "dfa92908878b2534ed72741806ac4b2ad6584fe6722aa06a97e8a29f0b61c876"

EXPECTED_WARNING_RECORDS = [
    {"species": "C2H4OOH", "property": "cp/R", "Tmid_K": 1000.0},
    {"species": "C2H4OOH", "property": "h/RT", "Tmid_K": 1000.0},
    {"species": "OC2H3OOH", "property": "cp/R", "Tmid_K": 1000.0},
    {"species": "OC2H3OOH", "property": "h/RT", "Tmid_K": 1000.0},
    {"species": "CH2CH2OH", "property": "cp/R", "Tmid_K": 1000.0},
    {"species": "CH2CH2OH", "property": "h/RT", "Tmid_K": 1000.0},
    {"species": "CH3CHOH", "property": "cp/R", "Tmid_K": 1000.0},
    {"species": "CH3CHOH", "property": "h/RT", "Tmid_K": 1000.0},
    {"species": "CH3CH2O", "property": "cp/R", "Tmid_K": 1000.0},
    {"species": "CH3CH2O", "property": "h/RT", "Tmid_K": 1000.0},
    {"species": "C3H6OOH", "property": "cp/R", "Tmid_K": 1000.0},
    {"species": "C3H6OOH", "property": "h/RT", "Tmid_K": 1000.0},
    {"species": "OC3H5OOH", "property": "cp/R", "Tmid_K": 1000.0},
    {"species": "OC3H5OOH", "property": "h/RT", "Tmid_K": 1000.0},
]


class CandidateError(ValueError):
    """A contract violation that must be reported without a traceback."""


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise CandidateError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def _decode_json(text: str, label: str) -> Any:
    try:
        return json.loads(text, object_pairs_hook=_reject_duplicate_pairs)
    except CandidateError:
        raise
    except json.JSONDecodeError as exc:
        raise CandidateError(f"{label}-json-invalid:{exc.msg}") from exc


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise CandidateError(f"artifact-open-failed:{path}:{exc}") from exc
    return digest.hexdigest()


def _valid_sha(value: Any, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise CandidateError(f"invalid-sha256:{label}")
    value = value.lower()
    if any(c not in "0123456789abcdef" for c in value) or value == "0" * 64:
        raise CandidateError(f"invalid-sha256:{label}")
    return value


def _object(value: Any, label: str, keys: set[str]) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise CandidateError(f"not-object:{label}")
    actual = set(value)
    missing = keys - actual
    unknown = actual - keys
    if missing:
        raise CandidateError(f"missing-keys:{label}:" + ",".join(sorted(missing)))
    if unknown:
        raise CandidateError(f"unknown-keys:{label}:" + ",".join(sorted(unknown)))
    return value


def _string(mapping: dict[str, Any], key: str, label: str) -> str:
    value = mapping[key]
    if not isinstance(value, str) or not value:
        raise CandidateError(f"invalid-string:{label}.{key}")
    return value


def _bool(mapping: dict[str, Any], key: str, label: str) -> bool:
    value = mapping[key]
    if type(value) is not bool:
        raise CandidateError(f"invalid-bool:{label}.{key}")
    return value


def _int(mapping: dict[str, Any], key: str, label: str) -> int:
    value = mapping[key]
    if type(value) is not int:
        raise CandidateError(f"invalid-int:{label}.{key}")
    return value


def _exact(actual: Any, expected: Any, label: str) -> None:
    if actual != expected or type(actual) is not type(expected):
        raise CandidateError(f"contract-mismatch:{label}")


def _workspace_path(workspace_root: Path, relative: str, label: str) -> Path:
    if not isinstance(relative, str) or not relative or "\\" in relative:
        raise CandidateError(f"invalid-workspace-path:{label}")
    if relative.startswith("/") or re.match(r"^[A-Za-z]:", relative):
        raise CandidateError(f"non-relative-workspace-path:{label}")
    candidate = (workspace_root / relative).resolve()
    root = workspace_root.resolve()
    if candidate != root and root not in candidate.parents:
        raise CandidateError(f"path-escapes-workspace:{label}")
    return candidate


def _absolute_path(value: str, label: str) -> Path:
    if not isinstance(value, str) or not re.match(r"^[A-Za-z]:/", value):
        raise CandidateError(f"non-canonical-absolute-path:{label}")
    return Path(value)


def _normalize_yaml(path: Path) -> bytes:
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise CandidateError(f"yaml-read-failed:{path}:{exc}") from exc
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    text = re.sub(r"(?m)^date:.*(?:\n|$)", "", text)
    if not text.endswith("\n"):
        text += "\n"
    return text.encode("utf-8")


def _normalized_yaml_sha(path: Path) -> str:
    return hashlib.sha256(_normalize_yaml(path)).hexdigest()


def _validate_artifacts(workspace_root: Path, value: Any) -> dict[str, Path]:
    if not isinstance(value, list) or len(value) != len(WORKSPACE_ARTIFACTS):
        raise CandidateError("artifact-list-contract-mismatch")
    result: dict[str, Path] = {}
    seen_paths: set[str] = set()
    for index, raw in enumerate(value):
        entry = _object(raw, f"artifacts[{index}]", {"id", "path", "url", "sha256", "kind"})
        artifact_id = _string(entry, "id", f"artifacts[{index}]")
        if artifact_id in result or artifact_id not in WORKSPACE_ARTIFACTS:
            raise CandidateError(f"unknown-or-duplicate-artifact:{artifact_id}")
        if artifact_id != list(WORKSPACE_ARTIFACTS)[index]:
            raise CandidateError(f"artifact-order-contract-mismatch:{artifact_id}")
        expected = WORKSPACE_ARTIFACTS[artifact_id]
        path = _string(entry, "path", f"artifact[{artifact_id}]")
        if path != expected["path"]:
            raise CandidateError(f"artifact-path-contract-mismatch:{artifact_id}")
        if path in seen_paths:
            raise CandidateError(f"duplicate-artifact-path:{path}")
        seen_paths.add(path)
        _exact(entry["url"], expected["url"], f"artifact[{artifact_id}].url")
        _exact(entry["kind"], expected["kind"], f"artifact[{artifact_id}].kind")
        expected_sha = _valid_sha(entry["sha256"], f"artifact[{artifact_id}].sha256")
        if expected_sha != expected["sha256"]:
            raise CandidateError(f"artifact-sha-contract-mismatch:{artifact_id}")
        actual_path = _workspace_path(workspace_root, path, f"artifact[{artifact_id}]")
        if not actual_path.is_file():
            raise CandidateError(f"artifact-missing:{artifact_id}:{path}")
        if _sha256(actual_path) != expected_sha:
            raise CandidateError(f"artifact-sha256-mismatch:{artifact_id}")
        result[artifact_id] = actual_path
    if set(result) != set(WORKSPACE_ARTIFACTS):
        raise CandidateError("artifact-set-contract-mismatch")
    return result


def _validate_upstream(workspace_root: Path, value: Any) -> None:
    if not isinstance(value, list) or len(value) != len(UPSTREAM_ARTIFACTS):
        raise CandidateError("upstream-artifact-list-contract-mismatch")
    seen: set[str] = set()
    for index, raw in enumerate(value):
        entry = _object(raw, f"upstream[{index}]", {"id", "path", "sha256"})
        artifact_id = _string(entry, "id", f"upstream[{index}]")
        if artifact_id in seen or artifact_id not in UPSTREAM_ARTIFACTS:
            raise CandidateError(f"unknown-or-duplicate-upstream:{artifact_id}")
        if artifact_id != list(UPSTREAM_ARTIFACTS)[index]:
            raise CandidateError(f"upstream-order-contract-mismatch:{artifact_id}")
        seen.add(artifact_id)
        expected = UPSTREAM_ARTIFACTS[artifact_id]
        if entry["path"] != expected["path"]:
            raise CandidateError(f"upstream-path-contract-mismatch:{artifact_id}")
        if _valid_sha(entry["sha256"], f"upstream[{artifact_id}].sha256") != expected["sha256"]:
            raise CandidateError(f"upstream-sha-contract-mismatch:{artifact_id}")
        path = _workspace_path(workspace_root, entry["path"], f"upstream[{artifact_id}]")
        if not path.is_file() or _sha256(path) != expected["sha256"]:
            raise CandidateError(f"upstream-sha256-mismatch:{artifact_id}")
    if seen != set(UPSTREAM_ARTIFACTS):
        raise CandidateError("upstream-artifact-set-contract-mismatch")


def _validate_runtime_artifacts(value: Any) -> None:
    if not isinstance(value, list) or len(value) != len(RUNTIME_ARTIFACTS):
        raise CandidateError("runtime-artifact-list-contract-mismatch")
    seen: set[str] = set()
    for index, raw in enumerate(value):
        entry = _object(raw, f"conversion.runtime_artifacts[{index}]", {"id", "path", "sha256"})
        artifact_id = _string(entry, "id", f"runtime_artifacts[{index}]")
        if artifact_id in seen or artifact_id not in RUNTIME_ARTIFACTS:
            raise CandidateError(f"unknown-or-duplicate-runtime-artifact:{artifact_id}")
        if artifact_id != list(RUNTIME_ARTIFACTS)[index]:
            raise CandidateError(f"runtime-artifact-order-contract-mismatch:{artifact_id}")
        seen.add(artifact_id)
        expected = RUNTIME_ARTIFACTS[artifact_id]
        if entry["path"] != expected["path"]:
            raise CandidateError(f"runtime-artifact-path-contract-mismatch:{artifact_id}")
        if _valid_sha(entry["sha256"], f"runtime_artifacts[{artifact_id}].sha256") != expected["sha256"]:
            raise CandidateError(f"runtime-artifact-sha-contract-mismatch:{artifact_id}")
        path = _absolute_path(entry["path"], f"runtime_artifacts[{artifact_id}]")
        if not path.is_file() or _sha256(path) != expected["sha256"]:
            raise CandidateError(f"runtime-artifact-sha256-mismatch:{artifact_id}")
    if seen != set(RUNTIME_ARTIFACTS):
        raise CandidateError("runtime-artifact-set-contract-mismatch")


def _validate_mechanism(document: dict[str, Any], artifacts: dict[str, Path]) -> None:
    candidate = _object(
        document["mechanism_candidate"],
        "mechanism_candidate",
        {"id", "status", "role", "source_authority", "artifacts", "license", "declared_domain"},
    )
    _exact(candidate["id"], "ucsd-san-diego-2016-12-14", "mechanism_candidate.id")
    _exact(candidate["status"], "candidate-only", "mechanism_candidate.status")
    _exact(candidate["role"], "gas-phase-detailed-mechanism-reference-candidate", "mechanism_candidate.role")
    authority = _object(candidate["source_authority"], "source_authority", {"landing_page_url", "citation", "update_log_url", "update_log_sha256"})
    _exact(authority["landing_page_url"], EXPECTED_SOURCE_PAGE, "source_authority.landing_page_url")
    _exact(authority["citation"], "Chemical-Kinetic Mechanisms for Combustion Applications; San Diego Mechanism web page; UC San Diego Combustion Research", "source_authority.citation")
    _exact(authority["update_log_url"], EXPECTED_UPDATE_LOG_URL, "source_authority.update_log_url")
    if _valid_sha(authority["update_log_sha256"], "source_authority.update_log_sha256") != EXPECTED_UPDATE_LOG_SHA:
        raise CandidateError("source-authority-update-log-sha-contract-mismatch")
    if candidate["artifacts"] != ["ucsd_mechanism", "ucsd_thermo", "ucsd_transport"]:
        raise CandidateError("mechanism-artifact-order-contract-mismatch")
    license_data = _object(candidate["license"], "mechanism_candidate.license", {"status", "evidence", "fit_use_allowed", "required_before_selection"})
    expected_license = {
        "status": "unverified",
        "evidence": "source page provides citation and copyright notice but no explicit redistribution or validation license",
        "fit_use_allowed": False,
        "required_before_selection": ["obtain explicit use/redistribution terms", "record source version and citation", "independent mechanism review"],
    }
    _exact(license_data, expected_license, "mechanism_candidate.license")
    domain = _object(candidate["declared_domain"], "mechanism_candidate.declared_domain", {"phase_model", "temperature_k", "pressure", "composition", "cryogenic_vle", "solid_carbon"})
    expected_domain = {
        "phase_model": "ideal-gas",
        "temperature_k": [300.0, 3000.0],
        "pressure": "not-declared",
        "composition": "not-declared",
        "cryogenic_vle": "unsupported",
        "solid_carbon": "not-modeled",
    }
    _exact(domain, expected_domain, "mechanism_candidate.declared_domain")
    for artifact_id in candidate["artifacts"]:
        if artifact_id not in artifacts:
            raise CandidateError(f"mechanism-artifact-missing:{artifact_id}")


def _validate_conversion(document: dict[str, Any], artifacts: dict[str, Path]) -> None:
    conversion = _object(
        document["conversion"],
        "conversion",
        {"status", "input_artifact_ids", "generator", "command", "output", "regeneration"},
    )
    _exact(conversion["status"], "reproduced-date-normalized", "conversion.status")
    _exact(conversion["input_artifact_ids"], ["ucsd_mechanism", "ucsd_thermo", "ucsd_transport"], "conversion.input_artifact_ids")
    generator = _object(conversion["generator"], "conversion.generator", {"name", "version", "python_version", "python_root", "wheel_artifact_id", "runtime_artifacts"})
    _exact(generator["name"], "Cantera ck2yaml", "conversion.generator.name")
    _exact(generator["version"], "3.2.0", "conversion.generator.version")
    _exact(generator["python_version"], "3.12.10", "conversion.generator.python_version")
    _exact(generator["python_root"], "F:/tmp/cantera-3.2.0-python312", "conversion.generator.python_root")
    _exact(generator["wheel_artifact_id"], "cantera_wheel", "conversion.generator.wheel_artifact_id")
    _validate_runtime_artifacts(generator["runtime_artifacts"])
    command = _object(conversion["command"], "conversion.command", {"working_directory", "argv", "environment"})
    _exact(command["working_directory"], "Vortex/tools/openfoam", "conversion.command.working_directory")
    _exact(command["argv"], ["python", "-m", "cantera.ck2yaml", "--input=sandiego20161214_mechCK.txt", "--thermo=sandiego20160815_therm.txt", "--transport=sandiego20160815_trans.txt", "--output=sandiego2016_c3h8.yaml"], "conversion.command.argv")
    _exact(command["environment"], {"PYTHONPATH": "F:/tmp/cantera-3.2.0-python312"}, "conversion.command.environment")
    output = _object(conversion["output"], "conversion.output", {"artifact_id", "normalized_sha256", "normalization"})
    _exact(output["artifact_id"], "converted_yaml", "conversion.output.artifact_id")
    _exact(output["normalization"], "UTF-8; CRLF-to-LF; remove top-level date line; final LF", "conversion.output.normalization")
    if _valid_sha(output["normalized_sha256"], "conversion.output.normalized_sha256") != EXPECTED_CANONICAL_YAML_SHA:
        raise CandidateError("conversion-normalized-sha-contract-mismatch")
    if _normalized_yaml_sha(artifacts["converted_yaml"]) != EXPECTED_CANONICAL_YAML_SHA:
        raise CandidateError("conversion-normalized-sha256-mismatch")
    regeneration = _object(conversion["regeneration"], "conversion.regeneration", {"artifact_id", "normalized_sha256", "date_line_may_change"})
    _exact(regeneration["artifact_id"], "regenerated_yaml", "conversion.regeneration.artifact_id")
    if _valid_sha(regeneration["normalized_sha256"], "conversion.regeneration.normalized_sha256") != EXPECTED_CANONICAL_YAML_SHA:
        raise CandidateError("conversion-regeneration-normalized-sha-contract-mismatch")
    if regeneration["date_line_may_change"] is not True:
        raise CandidateError("conversion-regeneration-date-policy-mismatch")
    if _normalized_yaml_sha(artifacts["regenerated_yaml"]) != EXPECTED_CANONICAL_YAML_SHA:
        raise CandidateError("conversion-regeneration-normalized-sha256-mismatch")


def _validate_runtime_probe(document: dict[str, Any]) -> None:
    runtime = _object(document["runtime_probe"], "runtime_probe", {"status", "launcher", "python", "cantera", "mechanism_artifact_id", "observed"})
    _exact(runtime["status"], "observed-not-physical-validation", "runtime_probe.status")
    _exact(runtime["launcher"], {"kind": "wsl", "argv": ["wsl", "python3", "-"], "distribution": "default"}, "runtime_probe.launcher")
    _exact(runtime["python"], {"version": "3.8.10", "executable": "/usr/bin/python3", "sha256": "298a9e830ed52f36c299427565485d717d1ce0179c0597cc16560513eb780b06"}, "runtime_probe.python")
    expected_cantera = {
        "version": "3.1.0",
        "module": "/home/nvoyager/.local/lib/python3.8/site-packages/cantera/__init__.py",
        "module_sha256": "9a7ceae3bae1a7e41dc581e353d44e04af8447114748f5b7c75995e26f8819b9",
        "extension": "/home/nvoyager/.local/lib/python3.8/site-packages/cantera/_cantera.cpython-38-x86_64-linux-gnu.so",
        "extension_sha256": "c57a47c201b0cb8f965b38acc7ebfafbbe95db7464d80359499c64d50555fae8",
        "metadata_sha256": "6593b4da73ead059f489af664cca95f401a041c178eb1a1b3b8fa02d2ae26b0c",
        "record_sha256": "994d183e5592f0dcaef047c5f6c369e821490f58e1e87123fa12d9da9db964b8",
    }
    _exact(runtime["cantera"], expected_cantera, "runtime_probe.cantera")
    _exact(runtime["mechanism_artifact_id"], "converted_yaml", "runtime_probe.mechanism_artifact_id")
    observed = _object(runtime["observed"], "runtime_probe.observed", {"species_count", "reaction_count", "reaction_types", "duplicate_reactions", "temperature_k", "species_order_sha256", "equation_order_sha256", "warnings"})
    expected_observed = {
        "species_count": 59,
        "reaction_count": 271,
        "reaction_types": {"Arrhenius": 233, "falloff-Troe": 23, "three-body-Arrhenius": 15},
        "duplicate_reactions": 6,
        "temperature_k": [300.0, 3000.0],
        "species_order_sha256": EXPECTED_RUNTIME_SPECIES_SHA,
        "equation_order_sha256": EXPECTED_RUNTIME_EQUATION_SHA,
    }
    for key, expected in expected_observed.items():
        _exact(observed[key], expected, f"runtime_probe.observed.{key}")
    warnings_data = _object(observed["warnings"], "runtime_probe.observed.warnings", {"count", "normalization", "records", "sha256"})
    _exact(warnings_data["count"], 14, "runtime_probe.observed.warnings.count")
    _exact(warnings_data["normalization"], "parsed species/property/Tmid records; sorted-key compact JSON plus LF", "runtime_probe.observed.warnings.normalization")
    _exact(warnings_data["records"], EXPECTED_WARNING_RECORDS, "runtime_probe.observed.warnings.records")
    if _valid_sha(warnings_data["sha256"], "runtime_probe.observed.warnings.sha256") != EXPECTED_RUNTIME_WARNING_SHA:
        raise CandidateError("runtime-warning-sha-contract-mismatch")


def _validate_experiments(document: dict[str, Any]) -> None:
    data = _object(document["experimental_evidence"], "experimental_evidence", {"status", "candidate_literature", "training_dataset", "holdout_dataset", "selection_requirements"})
    _exact(data["status"], "unselected", "experimental_evidence.status")
    expected_literature = [
        {"doi": "10.1016/j.proci.2010.06.131", "status": "metadata-only-no-raw-dataset", "url": "https://doi.org/10.1016/j.proci.2010.06.131"},
        {"doi": "10.1016/j.combustflame.2023.112826", "status": "metadata-only-no-raw-dataset", "url": "https://doi.org/10.1016/j.combustflame.2023.112826"},
        {"doi": "10.1016/j.combustflame.2019.10.015", "status": "metadata-only-no-raw-dataset", "url": "https://doi.org/10.1016/j.combustflame.2019.10.015"},
    ]
    _exact(data["candidate_literature"], expected_literature, "experimental_evidence.candidate_literature")
    expected_dataset = {"status": "unselected", "artifact_ids": [], "license": "unselected", "units": "unselected", "conditions": "unselected", "observables": [], "uncertainty": "unselected"}
    _exact(data["training_dataset"], expected_dataset, "experimental_evidence.training_dataset")
    _exact(data["holdout_dataset"], expected_dataset, "experimental_evidence.holdout_dataset")
    expected_requirements = [
        "raw machine-readable data and immutable source hash",
        "explicit license or permission for fitting and redistribution",
        "temperature, pressure, composition and equivalence-ratio units",
        "measurement uncertainty and protocol metadata",
        "disjoint training and holdout partitions frozen before fitting",
        "predeclared observables, metrics and tolerances",
    ]
    _exact(data["selection_requirements"], expected_requirements, "experimental_evidence.selection_requirements")


def _validate_protocol(document: dict[str, Any]) -> None:
    protocol = _object(document["validation_protocol"], "validation_protocol", {"status", "implementation_layer", "physical_layer", "metrics", "applicability_envelope", "fit_parameters", "stop_conditions"})
    _exact(protocol["status"], "not-run", "validation_protocol.status")
    _exact(protocol["implementation_layer"], {"status": "same-mechanism-only-upstream", "required_observables": ["forward_rate", "reverse_rate", "reaction_orders", "rate_units", "formation_energy_qdot", "mass_residual", "CHON_residual", "energy_residual"], "upstream_evidence": "WP0/WP2 same-mechanism instantaneous equivalence; not physical validation"}, "validation_protocol.implementation_layer")
    _exact(protocol["physical_layer"], {"status": "not-run", "required_comparisons": ["detailed-mechanism trajectories", "independent experimental trajectories"], "holdout_required": True}, "validation_protocol.physical_layer")
    _exact(protocol["metrics"], {"status": "not-frozen", "required": ["ignition_delay", "temperature_history", "major_species_history", "heat_release_history"], "tolerances": "must be selected before fitting; no post-hoc tolerance", "basis": "common control-volume and thermodynamic convention required"}, "validation_protocol.metrics")
    _exact(protocol["applicability_envelope"], {"status": "incomplete", "required_dimensions": ["temperature", "pressure", "composition", "equivalence_ratio", "dilution", "residence_time", "uncertainty"], "temperature_k": [300.0, 3000.0], "pressure_pa": "not-declared", "equivalence_ratio": "not-declared", "dilution": "not-declared", "residence_time_s": "not-declared", "phase_partition": "gas-only candidate; VLE/near-critical/solid-risk excluded"}, "validation_protocol.applicability_envelope")
    _exact(protocol["fit_parameters"], {"status": "not-run", "parameters": ["pre-exponential factors", "activation energies", "reaction orders", "formation-energy basis"], "method": "not-selected"}, "validation_protocol.fit_parameters")
    _exact(protocol["stop_conditions"], ["missing or ambiguous source/license", "changed mechanism/runtime fingerprint", "missing raw experiment or partition", "noncomparable energy or rate basis", "holdout metric failure", "out-of-envelope or phase-state mismatch"], "validation_protocol.stop_conditions")


def _validate_cross_scope(document: dict[str, Any]) -> None:
    cross = _object(document["cross_scope"], "cross_scope", {"kij", "hem", "condensed_carbon", "comparison"})
    _exact(cross["kij"], {"owner": "WP3", "status": "incomplete", "missing_pair_policy": "fail-closed", "refprop_is_pr_reference": False}, "cross_scope.kij")
    _exact(cross["hem"], {"owner": "RFAM_RHO_E_Y_FLASH", "status": "not-owned-by-detailed-mechanism-candidate", "quality_uds_authority": False}, "cross_scope.hem")
    _exact(cross["condensed_carbon"], {"owner": "A408_DILUTE_ZERO_SLIP_LEDGER_ONLY", "status": "wiring-or-diagnostic-only", "physics_validation": "not-run"}, "cross_scope.condensed_carbon")
    _exact(cross["comparison"], {"quantitative_domain": "validated-single-phase-fluid-only", "excluded": ["quality", "latent_energy", "rhoSoot", "number_density", "wall_carbon", "two_phase_mixed_density", "two_phase_phase_composition"]}, "cross_scope.comparison")


def validate(workspace_root: Path, manifest_path: Path) -> tuple[dict[str, Any], dict[str, Path], int]:
    try:
        text = manifest_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise CandidateError(f"manifest-open-failed:{exc}") from exc
    document = _decode_json(text, "manifest")
    document = _object(document, "root", EXPECTED_ROOT_KEYS)
    _exact(document["schema"], SCHEMA, "root.schema")
    if _int(document, "contract_version", "root") != 1:
        raise CandidateError("contract-version-mismatch")
    _exact(document["status"], "blocked", "root.status")
    _exact(document["scope"], "provenance-screening-only-no-fitting-or-physical-validation", "root.scope")
    approvals = _object(document["approvals"], "approvals", {"mechanism_selected_for_fit", "r2_fitting_approved", "physical_validation_approved", "production_ready", "reacting_hot_start_approved"})
    for key in approvals:
        if _bool(approvals, key, "approvals"):
            raise CandidateError(f"approval-must-remain-false:{key}")
    _validate_upstream(workspace_root, document["upstream"])
    artifacts = _validate_artifacts(workspace_root, document["artifacts"])
    _validate_mechanism(document, artifacts)
    _validate_conversion(document, artifacts)
    _validate_runtime_probe(document)
    _validate_experiments(document)
    _validate_protocol(document)
    _validate_cross_scope(document)
    _exact(document["limitations"], ["candidate is ideal-gas detailed chemistry only", "14 NASA7 continuity warnings are recorded and unresolved", "300-3000 K mechanism intersection does not cover cryogenic 130 K injection", "no experimental raw data, license, partition or tolerance is selected", "Cantera is an implementation/reference tool, not proof of physical R2 validity", "no HEM quality, PR kij, soot physics, Fluent or hot-start capability is enabled"], "limitations")
    return document, artifacts, 1


def _network_sha(url: str) -> str:
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            data = response.read()
    except Exception as exc:
        raise CandidateError(f"official-download-failed:{url}:{exc}") from exc
    return hashlib.sha256(data).hexdigest()


def verify_official_downloads() -> int:
    for artifact_id in ("ucsd_mechanism", "ucsd_thermo", "ucsd_transport"):
        expected = WORKSPACE_ARTIFACTS[artifact_id]
        if _network_sha(expected["url"]) != expected["sha256"]:
            raise CandidateError(f"official-sha256-mismatch:{artifact_id}")
    if _network_sha(EXPECTED_UPDATE_LOG_URL) != EXPECTED_UPDATE_LOG_SHA:
        raise CandidateError("official-sha256-mismatch:update_log")
    return 4


def _wsl_path(path: Path) -> str:
    normalized = path.resolve().as_posix()
    match = re.match(r"^([A-Za-z]):/(.*)$", normalized)
    if not match:
        raise CandidateError(f"cannot-map-wsl-path:{normalized}")
    return f"/mnt/{match.group(1).lower()}/{match.group(2)}"


def verify_runtime_probe(document: dict[str, Any], mechanism_path: Path) -> int:
    script = r'''
import cantera as ct
import hashlib
import importlib.metadata
import json
import platform
import re
import sys
import warnings
from collections import Counter

mechanism = __MECHANISM__

def sha(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def line_sha(values):
    return hashlib.sha256(('\n'.join(values) + '\n').encode()).hexdigest()

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    gas = ct.Solution(mechanism)

records = []
for warning in caught:
    text = str(warning.message)
    species = re.search(r'For species ([A-Za-z0-9_-]+)', text)
    prop = re.search(r'discontinuity in ([^ ]+) detected', text)
    tmid = re.search(r'Tmid = ([0-9.]+)', text)
    if not (species and prop and tmid):
        raise RuntimeError('unparseable-warning')
    records.append({'species': species.group(1), 'property': prop.group(1), 'Tmid_K': float(tmid.group(1))})

distribution = importlib.metadata.distribution('cantera')
extension = next(str(distribution.locate_file(f)) for f in distribution.files if str(f).endswith('_cantera.cpython-38-x86_64-linux-gnu.so'))
metadata = str(distribution.locate_file(next(f for f in distribution.files if str(f).endswith('.dist-info/METADATA'))))
record = str(distribution.locate_file(next(f for f in distribution.files if str(f).endswith('.dist-info/RECORD'))))
warning_sha = hashlib.sha256((json.dumps(records, sort_keys=True, separators=(',', ':')) + '\n').encode()).hexdigest()

output = {
    'python': {'version': platform.python_version(), 'executable': sys.executable, 'sha256': sha(sys.executable)},
    'cantera': {
        'version': ct.__version__, 'module': ct.__file__, 'module_sha256': sha(ct.__file__),
        'extension': extension, 'extension_sha256': sha(extension),
        'metadata_sha256': sha(metadata), 'record_sha256': sha(record),
    },
    'observed': {
        'species_count': gas.n_species, 'reaction_count': gas.n_reactions,
        'reaction_types': dict(sorted(Counter(r.reaction_type for r in gas.reactions()).items())),
        'duplicate_reactions': sum(bool(r.duplicate) for r in gas.reactions()),
        'temperature_k': [gas.min_temp, gas.max_temp],
        'species_order_sha256': line_sha(list(gas.species_names)),
        'equation_order_sha256': line_sha([r.equation for r in gas.reactions()]),
        'warnings': {
            'count': len(records),
            'normalization': 'parsed species/property/Tmid records; sorted-key compact JSON plus LF',
            'records': records, 'sha256': warning_sha,
        },
    },
}
print(json.dumps(output, sort_keys=True))
'''.replace("__MECHANISM__", repr(_wsl_path(mechanism_path)))
    process = subprocess.run(
        ["wsl", "python3", "-"], input=script, text=True,
        capture_output=True, check=False,
    )
    if process.returncode != 0:
        raise CandidateError(f"runtime-probe-command-failed:{process.returncode}:{process.stderr.strip()}")
    actual = _decode_json(process.stdout.strip(), "runtime-probe")
    runtime = document["runtime_probe"]
    expected = {
        "python": runtime["python"],
        "cantera": runtime["cantera"],
        "observed": runtime["observed"],
    }
    if actual != expected:
        raise CandidateError("runtime-probe-contract-mismatch")
    return 1


def verify_regeneration(repo_root: Path) -> int:
    python_path = Path(RUNTIME_ARTIFACTS["cantera320_python"]["path"])
    input_dir = (repo_root.parent / "Vortex/tools/openfoam").resolve()
    with tempfile.TemporaryDirectory(prefix="rfam-r2-regeneration-") as temp:
        output = Path(temp) / "sandiego2016_c3h8.yaml"
        argv = [
            str(python_path), "-m", "cantera.ck2yaml",
            "--input=sandiego20161214_mechCK.txt",
            "--thermo=sandiego20160815_therm.txt",
            "--transport=sandiego20160815_trans.txt",
            f"--output={output}",
        ]
        env = os.environ.copy()
        env["PYTHONPATH"] = "F:/tmp/cantera-3.2.0-python312"
        process = subprocess.run(
            argv, cwd=input_dir, env=env, text=True,
            capture_output=True, check=False,
        )
        if process.returncode != 0:
            raise CandidateError(f"regeneration-command-failed:{process.returncode}:{process.stderr.strip()}")
        if _normalized_yaml_sha(output) != EXPECTED_CANONICAL_YAML_SHA:
            raise CandidateError("regeneration-normalized-sha256-mismatch")
    return 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--contract", type=Path, default=Path("tests/reference/r2_physics_candidate_v1.json"))
    parser.add_argument("--verify-official-downloads", action="store_true")
    parser.add_argument("--verify-runtime-probe", action="store_true")
    parser.add_argument("--verify-regeneration", action="store_true")
    args = parser.parse_args(argv)
    try:
        repo_root = args.root.resolve()
        contract = args.contract
        if not contract.is_absolute():
            contract = (repo_root / contract).resolve()
        document, artifacts, checks = validate(repo_root.parent, contract)
        if args.verify_official_downloads:
            checks += verify_official_downloads()
        if args.verify_runtime_probe:
            checks += verify_runtime_probe(document, artifacts["converted_yaml"])
        if args.verify_regeneration:
            checks += verify_regeneration(repo_root)
    except Exception as exc:  # Fail closed without a traceback for malformed input.
        reason = str(exc).replace("\r", " ").replace("\n", " ")
        print(f"R2_PHYSICS_CANDIDATE_FAIL reason={reason}")
        return 2
    print(f"R2_PHYSICS_CANDIDATE_BLOCKED checks={checks}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
