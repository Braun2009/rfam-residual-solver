#!/usr/bin/env python3
"""Validate the current-source synthetic RFAM kij objective probe."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tools"))
import kij_fit_harness as harness  # noqa: E402


SCHEMA = "rfam.kij-fit-objective-contract/v2"
REFERENCE_ONLY_MARKER = "RFAM_KIJ_FIT_OBJECTIVE_V2_REFERENCE_ONLY_BLOCKED"
PASS_MARKER = "RFAM_KIJ_FIT_OBJECTIVE_V2_SYNTHETIC_TEST_ONLY_PASS"
FAIL_MARKER = "RFAM_KIJ_FIT_OBJECTIVE_V2_CONTRACT_FAIL"

EXPECTED_HARNESS_TOOL_SHA256 = (
    "e22ae302b485fb7075f5149538a4ad9889de5f46b7ef284f67910b6fad316726"
)
EXPECTED_HARNESS_MANIFEST_SHA256 = (
    "a60d3adf8a02e6fae57af37f49f1c8652737bb58a156c5785c7e117dadfd8480"
)
EXPECTED_REFERENCE_SHA256 = (
    "45eb0ebc984d95d4b8cf5c499d466574fa0c6bb4f8fe29bebb22aff89398738e"
)
EXPECTED_EXECUTABLE_SHA256 = (
    "f55a71b353abfa8d2f994a4b6af74dcd6437225182c068f812ee4b7bde648b7f"
)
EXPECTED_SOURCE_BUNDLE_SHA256 = (
    "6a97953eedf2a1c931ac53c4da2a40d81ec0a9ec6cb12463ea2c9ac827294883"
)

EXPECTED_TP_PH_OWNER_REPORT_SHA256 = (
    "a1392524c04615e6f2b13172d1916eab49c943aebe5c54baf1bef47c820e02ee"
)
EXPECTED_TP_PH_REVIEW_SHA256 = (
    "988643161a70ded8b8519325a33b5a2640077d74b96e3c1afe0db5b86305af56"
)
EXPECTED_CALORIC_80K_REPORT_SHA256 = (
    "19a6c536616bd53b6b465fcd8aa985aea2a522f13c321f43598c646ffb85f934"
)

HARNESS_MANIFEST_PATH = (
    "realFluidAllMach/tests/reference/kij_fit_harness_synthetic_manifest_v1.json"
)
REFERENCE_PATH = (
    "realFluidAllMach/tests/reference/kij_fit_objective_synthetic_v2.csv"
)
CONTRACT_PATH = (
    "realFluidAllMach/tests/reference/kij_fit_objective_contract_v2.json"
)

SOURCE_PATHS = (
    "realFluidAllMach/tools/kij_fit_objective_probe.cpp",
    "realFluidAllMach/src/thermo/ThermoDataProvenance.C",
    "realFluidAllMach/src/thermo/ThermoDataProvenance.H",
    "realFluidAllMach/src/thermo/SpeciesRegistry.H",
    "realFluidAllMach/src/thermo/NASA7Thermo.H",
    "realFluidAllMach/src/thermo/PengRobinsonAlpha.H",
    "realFluidAllMach/src/thermo/EquationOfState/PengRobinsonEOS.C",
    "realFluidAllMach/src/thermo/EquationOfState/PengRobinsonEOS.H",
    "realFluidAllMach/src/thermo/MixingRules/VdWOneFluidMixing.C",
    "realFluidAllMach/src/thermo/MixingRules/VdWOneFluidMixing.H",
    "realFluidAllMach/src/thermo/MixingRules/SpeciesCriticalDB.C",
    "realFluidAllMach/src/thermo/MixingRules/SpeciesCriticalDB.H",
    "realFluidAllMach/src/thermo/MixingRules/BinaryInteractionDB.C",
    "realFluidAllMach/src/thermo/MixingRules/BinaryInteractionDB.H",
    "realFluidAllMach/src/thermo/RootSelection/CubicRoots.C",
    "realFluidAllMach/src/thermo/RootSelection/CubicRoots.H",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/ConservedPhaseMixture.H",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/PRMixtureKernel.C",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/PRMixtureKernel.H",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/PRTPFlash.C",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/PRTPFlash.H",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/TabulatedIdealCaloric.C",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/TabulatedIdealCaloric.H",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/PRPhaseCaloric.C",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/PRPhaseCaloric.H",
)

TRANSLATION_UNITS = tuple(
    path for path in SOURCE_PATHS if path.endswith((".cpp", ".C"))
)

INCLUDE_DIRECTORIES = (
    "realFluidAllMach/src",
    "realFluidAllMach/src/thermo",
    "realFluidAllMach/src/thermo/MixingRules",
    "realFluidAllMach/src/thermo/EquationOfState",
    "realFluidAllMach/src/thermo/RootSelection",
    "realFluidAllMach/src/thermo/PhaseEquilibrium",
)

INPUT_HEADER = (
    "watermark",
    "case_id",
    "split",
    "pressure_pa",
    "temperature_k",
    "pressure_standard_uncertainty_pa",
    "liquid_x_o2",
    "liquid_x_o2_standard_uncertainty",
    "vapor_y_o2",
    "vapor_y_o2_standard_uncertainty",
)

OUTPUT_HEADER = (
    "case_id",
    "split",
    "pressure_pa",
    "temperature_k",
    "observed_liquid_x_o2",
    "observed_vapor_y_o2",
    "reference_liquid_x_o2",
    "reference_vapor_y_o2",
    "trial_liquid_x_o2",
    "trial_vapor_y_o2",
    "reference_tie_line_max_abs_residual",
    "trial_log_fugacity_residual_o2",
    "trial_log_fugacity_residual_c3h8",
    "reference_rho_kg_m3",
    "trial_rho_kg_m3",
    "reference_e_j_kg",
    "trial_e_j_kg",
    "phase_equilibrium_chi2",
    "liquid_x_chi2",
    "vapor_y_chi2",
    "rho_chi2",
    "e_chi2",
    "total_chi2",
    "reference_material_residual",
    "trial_material_residual",
    "reference_energy_residual",
    "trial_energy_residual",
)

EXPECTED_APPROVALS = {
    "physical_data_selected": False,
    "fit_execution_approved": False,
    "kij_value_validated": False,
    "runtime_kij_emission_approved": False,
    "binary_interaction_db_write_approved": False,
    "hem_enabled": False,
    "fluent_execution_approved": False,
    "hot_start_approved": False,
    "production_ready": False,
}

EXPECTED_OBJECTIVE = {
    "data_class": "synthetic-nonphysical-test-only",
    "reference_test_kij": -0.05,
    "trial_test_kij": 0.15,
    "test_values_role": "unvalidated-endpoints-from-prior-offline-sensitivity-screen",
    "vapor_mole_fraction_for_reconstructed_feed": 0.5,
    "phase_equilibrium_residual": "log(x_i*phi_i_liquid)-log(y_i*phi_i_vapor)-at-observed-p-t-x-y",
    "phase_equilibrium_standard_scale": "pressure_standard_uncertainty_pa/pressure_pa",
    "liquid_x_standard_scale": "liquid_x_o2_standard_uncertainty",
    "vapor_y_standard_scale": "vapor_y_o2_standard_uncertainty",
    "density_relative_standard_uncertainty": 0.02,
    "energy_standard_uncertainty_j_kg": 5000.0,
    "aggregation": "sum-of-squared-normalized-residuals",
    "split_reporting": "training-and-holdout-separate-no-selection-from-holdout",
    "fit_executed": False,
    "parameter_values_publishable": False,
}

EXPECTED_ACCEPTANCE = {
    "required_training_records": 2,
    "required_holdout_records": 1,
    "maximum_reference_tie_line_absolute_residual": 5.0e-9,
    "maximum_material_residual": 1.0e-8,
    "maximum_energy_identity_residual": 1.0e-8,
    "require_finite_nonnegative_objective_components": True,
    "require_byte_identical_regeneration": True,
}

EXPECTED_LIMITATIONS = [
    "all input rows and both test kij values are synthetic or unvalidated test-only quantities",
    "the pressure-uncertainty-to-log-fugacity scale is a scaffold convention, not a validated statistical model",
    "no optimizer, parameter covariance, physical training data, independent physical holdout, or fit acceptance threshold is provided",
    "the probe mutates only local in-memory BinaryInteractionDB copies and cannot publish a runtime entry",
    "Cantera is not used and cannot supply PR kij evidence",
    "current TP/PH repair requires a separately bound external topology-envelope admission",
    "the 80 K C3H8 ideal-caloric extension is provenance and does not validate kij",
    "no HEM, Fluent, zero-advance, hot-start, or production capability is enabled",
]

EXPECTED_CORE_PROVENANCE = {
    "tp_ph_flash_source": {
        "path": "realFluidAllMach/src/thermo/PhaseEquilibrium/PRTPFlash.C",
        "sha256": "e16a747f84c6d6891c3154bb385cad8ce39c6233952d8b8e1ce4fee0e07130dc",
        "role": "accepted-tp-ph-stability-repair-current-source",
        "admission": "external-topology-envelope-required-non-three-phase-subset-only",
    },
    "ideal_caloric_80k_source": {
        "path": "realFluidAllMach/src/thermo/PhaseEquilibrium/TabulatedIdealCaloric.C",
        "sha256": "6a0e016e905375fa3f6f387d5b80e5bcf373c051707b452b5ee759872513d5a4",
        "role": "current-c3h8-ideal-caloric-lower-bound-80-k",
        "admission": "caloric-range-extension-not-kij-validation",
    },
}

ROOT_KEYS = (
    "schema",
    "contract_version",
    "status",
    "scope",
    "approvals",
    "upstream",
    "core_provenance",
    "source_bundle",
    "objective",
    "reference",
    "build",
    "acceptance",
    "limitations",
)

SUMMARY_RE = re.compile(
    r"^# RFAM_KIJ_FIT_OBJECTIVE_SYNTHETIC_TEST_ONLY_PASS "
    r"training_records=(\d+) holdout_records=(\d+) "
    r"training_objective=([^ ]+) holdout_objective=([^ ]+) "
    r"fit_executed=0 runtime_kij_emitted=0 "
    r"binary_interaction_db_modified=0 production_ready=0$"
)
SHA_RE = re.compile(r"^[0-9a-f]{64}$")


class ObjectiveContractError(ValueError):
    """A structured fail-close objective contract error."""


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ObjectiveContractError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _exact(actual: Any, expected: Any, label: str) -> None:
    if type(actual) is not type(expected):
        raise ObjectiveContractError(f"type-mismatch:{label}")
    if isinstance(expected, dict):
        if tuple(actual) != tuple(expected):
            missing = [key for key in expected if key not in actual]
            unknown = [key for key in actual if key not in expected]
            if missing:
                raise ObjectiveContractError(f"missing-keys:{label}:" + ",".join(missing))
            if unknown:
                raise ObjectiveContractError(f"unknown-keys:{label}:" + ",".join(unknown))
            raise ObjectiveContractError(f"key-order-mismatch:{label}")
        for key in expected:
            _exact(actual[key], expected[key], f"{label}.{key}")
        return
    if isinstance(expected, list):
        if len(actual) != len(expected):
            raise ObjectiveContractError(f"list-length-mismatch:{label}")
        for index, expected_item in enumerate(expected):
            _exact(actual[index], expected_item, f"{label}[{index}]")
        return
    if actual != expected:
        raise ObjectiveContractError(f"value-mismatch:{label}")


def _workspace_path(root: Path, relative: str, label: str) -> Path:
    if not isinstance(relative, str) or not relative or "\\" in relative:
        raise ObjectiveContractError(f"invalid-path:{label}")
    if relative.startswith("/") or re.match(r"^[A-Za-z]:", relative):
        raise ObjectiveContractError(f"path-not-relative:{label}")
    resolved_root = root.resolve()
    candidate = (resolved_root / relative).resolve()
    if candidate != resolved_root and resolved_root not in candidate.parents:
        raise ObjectiveContractError(f"path-escapes-root:{label}")
    return candidate


def _read(path: Path, label: str) -> bytes:
    try:
        return path.read_bytes()
    except OSError as exc:
        raise ObjectiveContractError(f"read-failed:{label}:{exc}") from exc


def _decode_contract(path: Path) -> dict[str, Any]:
    payload = _read(path, "contract")
    if b"\r" in payload or not payload.endswith(b"\n"):
        raise ObjectiveContractError("contract-not-canonical-lf")
    try:
        text = payload.decode("ascii")
        document = json.loads(text, object_pairs_hook=_reject_duplicate_pairs)
    except ObjectiveContractError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ObjectiveContractError(f"contract-json-invalid:{exc}") from exc
    if not isinstance(document, dict):
        raise ObjectiveContractError("contract-root-not-object")
    if text != json.dumps(document, indent=2, ensure_ascii=True) + "\n":
        raise ObjectiveContractError("contract-not-canonical-json")
    return document


def source_bundle_sha(root: Path) -> str:
    digest = hashlib.sha256()
    for relative in SOURCE_PATHS:
        payload = _read(_workspace_path(root, relative, f"source:{relative}"), relative)
        encoded = relative.encode("ascii")
        digest.update(encoded)
        digest.update(b"\0")
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return digest.hexdigest()


def _expected_contract() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "contract_version": 2,
        "status": "current-source-synthetic-objective-scaffold-no-fit",
        "scope": "offline-current-source-rfam-pr-pxy-rho-e-objective-training-holdout-dry-run",
        "approvals": EXPECTED_APPROVALS,
        "upstream": [
            {
                "id": "kij_fit_harness_validator",
                "path": "realFluidAllMach/tools/kij_fit_harness.py",
                "sha256": EXPECTED_HARNESS_TOOL_SHA256,
            },
            {
                "id": "kij_fit_harness_synthetic_manifest",
                "path": HARNESS_MANIFEST_PATH,
                "sha256": EXPECTED_HARNESS_MANIFEST_SHA256,
            },
            {
                "id": "tp_ph_repair_owner_report",
                "path": "tmp/rfam-tp-ph-stability-repair-replacement-v2-20260723.md",
                "sha256": EXPECTED_TP_PH_OWNER_REPORT_SHA256,
            },
            {
                "id": "tp_ph_repair_independent_review",
                "path": "tmp/rfam-tp-ph-stability-repair-replacement-v2-independent-review-20260723.md",
                "sha256": EXPECTED_TP_PH_REVIEW_SHA256,
            },
            {
                "id": "caloric_80k_diagnosis",
                "path": "tmp/rfam-prph-three-phase-readonly-diagnosis-v1-20260723.md",
                "sha256": EXPECTED_CALORIC_80K_REPORT_SHA256,
            },
        ],
        "core_provenance": EXPECTED_CORE_PROVENANCE,
        "source_bundle": {
            "algorithm": "sha256(path-nul-size-u64be-bytes)",
            "files": list(SOURCE_PATHS),
            "sha256": EXPECTED_SOURCE_BUNDLE_SHA256,
        },
        "objective": EXPECTED_OBJECTIVE,
        "reference": {
            "path": REFERENCE_PATH,
            "sha256": EXPECTED_REFERENCE_SHA256,
            "byte_count": 2439,
            "rows": 3,
        },
        "build": {
            "platform": "WSL-Ubuntu-20.04-x86_64",
            "compiler_path": "/usr/bin/x86_64-linux-gnu-g++-9",
            "compiler_version": "g++ (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0",
            "compiler_sha256": "0c0b987719385f8e242819dfef6a3e97b4b2540ce4035dbc11565df2966bd69e",
            "language": "c++17",
            "optimization": "-O2",
            "openmp": "-fopenmp",
            "include_directories": list(INCLUDE_DIRECTORIES),
            "translation_units": list(TRANSLATION_UNITS),
            "expected_executable_sha256": EXPECTED_EXECUTABLE_SHA256,
        },
        "acceptance": EXPECTED_ACCEPTANCE,
        "limitations": EXPECTED_LIMITATIONS,
    }


def _validate_bound_file(root: Path, item: dict[str, Any]) -> None:
    path = _workspace_path(root, item["path"], f"upstream:{item['id']}")
    if _sha256(_read(path, item["id"])) != item["sha256"]:
        raise ObjectiveContractError(f"upstream-sha256-mismatch:{item['id']}")


def validate_contract(root: Path, contract_path: Path) -> dict[str, Any]:
    document = _decode_contract(contract_path)
    expected = _expected_contract()
    _exact(document, expected, "root")
    for item in document["upstream"]:
        _validate_bound_file(root, item)
    for provenance_id, item in document["core_provenance"].items():
        path = _workspace_path(root, item["path"], f"core:{provenance_id}")
        if _sha256(_read(path, f"core:{provenance_id}")) != item["sha256"]:
            raise ObjectiveContractError(
                f"core-provenance-sha256-mismatch:{provenance_id}"
            )
    if source_bundle_sha(root) != EXPECTED_SOURCE_BUNDLE_SHA256:
        raise ObjectiveContractError("source-bundle-sha256-mismatch")
    reference_path = _workspace_path(root, REFERENCE_PATH, "reference")
    payload = _read(reference_path, "reference")
    if len(payload) != document["reference"]["byte_count"]:
        raise ObjectiveContractError("reference-byte-count-mismatch")
    if _sha256(payload) != EXPECTED_REFERENCE_SHA256:
        raise ObjectiveContractError("reference-sha256-mismatch")
    return document


def _manifest_and_input(root: Path) -> tuple[dict[str, Any], bytes, list[tuple[str, str]]]:
    manifest_path = _workspace_path(root, HARNESS_MANIFEST_PATH, "harness-manifest")
    audit = harness.validate_manifest(
        root,
        manifest_path,
        EXPECTED_HARNESS_MANIFEST_SHA256,
        allow_synthetic_test_fixture=True,
    )
    if audit["data_class"] != harness.SYNTHETIC_CLASS:
        raise ObjectiveContractError("objective-requires-synthetic-harness-input")

    manifest = json.loads(manifest_path.read_text(encoding="ascii"))
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(INPUT_HEADER)
    identities: list[tuple[str, str]] = []
    for authority in manifest["artifacts"]:
        path = _workspace_path(root, authority["path"], authority["artifact_id"])
        with path.open("r", encoding="ascii", newline="") as stream:
            reader = csv.DictReader(stream)
            for row in reader:
                writer.writerow(
                    (
                        row["watermark"],
                        row["record_id"],
                        row["split"],
                        row["pressure_pa"],
                        row["temperature_k"],
                        row["pressure_standard_uncertainty_pa"],
                        row["liquid_x_o2"],
                        row["liquid_x_o2_standard_uncertainty"],
                        row["vapor_y_o2"],
                        row["vapor_y_o2_standard_uncertainty"],
                    )
                )
                identities.append((row["record_id"], row["split"]))
    return manifest, output.getvalue().encode("ascii"), identities


def _finite(text: str, label: str) -> float:
    try:
        value = float(text)
    except ValueError as exc:
        raise ObjectiveContractError(f"invalid-number:{label}") from exc
    if not math.isfinite(value):
        raise ObjectiveContractError(f"nonfinite-number:{label}")
    return value


def _close(lhs: float, rhs: float, tolerance: float = 1.0e-11) -> bool:
    return abs(lhs - rhs) <= tolerance*max(abs(lhs), abs(rhs), 1.0)


def validate_output(
    payload: bytes, identities: list[tuple[str, str]]
) -> dict[str, Any]:
    if b"\r" in payload or not payload.endswith(b"\n"):
        raise ObjectiveContractError("output-not-canonical-lf")
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as exc:
        raise ObjectiveContractError("output-not-ascii") from exc
    lines = text.splitlines()
    if len(lines) != len(identities) + 4:
        raise ObjectiveContractError("output-line-count-mismatch")
    _exact(lines[0], "# rfam.kij-fit-objective-probe/v1", "output.schema")
    _exact(
        lines[1],
        "# SYNTHETIC_NONPHYSICAL_TEST_ONLY reference_test_kij=-0.050000000000000003 trial_test_kij=0.14999999999999999 not_fitted=1 runtime_publishable=0",
        "output.watermark",
    )
    header = next(csv.reader([lines[2]], strict=True))
    _exact(header, list(OUTPUT_HEADER), "output.header")

    training_total = 0.0
    holdout_total = 0.0
    training_records = 0
    holdout_records = 0
    component_positive = {key: False for key in range(17, 22)}
    for index, expected_identity in enumerate(identities):
        fields = next(csv.reader([lines[index + 3]], strict=True))
        if len(fields) != len(OUTPUT_HEADER):
            raise ObjectiveContractError(f"output-column-count:row-{index + 1}")
        _exact((fields[0], fields[1]), expected_identity, f"output.identity[{index}]")
        values = [_finite(value, f"row-{index + 1}:{OUTPUT_HEADER[column]}") for column, value in enumerate(fields[2:], start=2)]
        del values  # Parsed for finiteness; named accesses below use original fields.
        tie_line = float(fields[10])
        if tie_line > EXPECTED_ACCEPTANCE["maximum_reference_tie_line_absolute_residual"]:
            raise ObjectiveContractError(f"reference-tie-line-residual:row-{index + 1}")
        components = [float(fields[column]) for column in range(17, 22)]
        if any(value < 0.0 for value in components):
            raise ObjectiveContractError(f"negative-objective-component:row-{index + 1}")
        for component_index, value in zip(range(17, 22), components):
            component_positive[component_index] |= value > 0.0
        total = float(fields[22])
        if not _close(total, sum(components)):
            raise ObjectiveContractError(f"objective-sum-mismatch:row-{index + 1}")
        for column in (23, 24):
            if float(fields[column]) > EXPECTED_ACCEPTANCE["maximum_material_residual"]:
                raise ObjectiveContractError(f"material-residual:row-{index + 1}")
        for column in (25, 26):
            if float(fields[column]) > EXPECTED_ACCEPTANCE["maximum_energy_identity_residual"]:
                raise ObjectiveContractError(f"energy-residual:row-{index + 1}")
        if fields[1] == "training":
            training_total += total
            training_records += 1
        else:
            holdout_total += total
            holdout_records += 1
    if not all(component_positive.values()):
        raise ObjectiveContractError("objective-component-not-exercised")

    match = SUMMARY_RE.fullmatch(lines[-1])
    if match is None:
        raise ObjectiveContractError("summary-marker-mismatch")
    if int(match.group(1)) != training_records or int(match.group(2)) != holdout_records:
        raise ObjectiveContractError("summary-split-count-mismatch")
    if not _close(float(match.group(3)), training_total) or not _close(float(match.group(4)), holdout_total):
        raise ObjectiveContractError("summary-objective-mismatch")
    if training_records != EXPECTED_ACCEPTANCE["required_training_records"]:
        raise ObjectiveContractError("training-count-mismatch")
    if holdout_records != EXPECTED_ACCEPTANCE["required_holdout_records"]:
        raise ObjectiveContractError("holdout-count-mismatch")
    return {
        "rows": len(identities),
        "training_records": training_records,
        "holdout_records": holdout_records,
        "training_objective": training_total,
        "holdout_objective": holdout_total,
    }


def run_probe(executable: Path, input_payload: bytes) -> bytes:
    payload = _read(executable, "probe-executable")
    if _sha256(payload) != EXPECTED_EXECUTABLE_SHA256:
        raise ObjectiveContractError("probe-executable-sha256-mismatch")
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    try:
        result = subprocess.run(
            [
                str(executable),
                "--reference-test-kij",
                str(EXPECTED_OBJECTIVE["reference_test_kij"]),
                "--trial-test-kij",
                str(EXPECTED_OBJECTIVE["trial_test_kij"]),
                "--density-relative-standard-uncertainty",
                str(EXPECTED_OBJECTIVE["density_relative_standard_uncertainty"]),
                "--energy-standard-uncertainty-j-kg",
                str(EXPECTED_OBJECTIVE["energy_standard_uncertainty_j_kg"]),
            ],
            input=input_payload,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=120,
            env=environment,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ObjectiveContractError(f"probe-execution-failed:{exc}") from exc
    if result.returncode != 0:
        error = result.stderr.decode("ascii", errors="replace").strip()
        raise ObjectiveContractError(f"probe-returned-{result.returncode}:{error}")
    if result.stderr:
        raise ObjectiveContractError("probe-stderr-not-empty")
    return result.stdout


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--probe-executable", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        root = args.root.resolve()
        contract = args.contract.resolve()
        validate_contract(root, contract)
        _, input_payload, identities = _manifest_and_input(root)
        reference = _read(_workspace_path(root, REFERENCE_PATH, "reference"), "reference")
        summary = validate_output(reference, identities)
        executed = args.probe_executable is not None
        if executed:
            actual = run_probe(args.probe_executable.resolve(), input_payload)
            validate_output(actual, identities)
            if actual != reference:
                raise ObjectiveContractError("probe-output-not-byte-identical")
    except Exception as exc:  # All malformed or unexpected input fails closed.
        reason = str(exc) if isinstance(exc, ObjectiveContractError) else f"unexpected:{type(exc).__name__}:{exc}"
        print(f"{FAIL_MARKER} reason={reason.replace(chr(13), ' ').replace(chr(10), ' ')}")
        return 2

    marker = PASS_MARKER if executed else REFERENCE_ONLY_MARKER
    print(
        f"{marker} rows={summary['rows']} training={summary['training_records']} "
        f"holdout={summary['holdout_records']} fit_executed=0 validated_kij=0 "
        "runtime_kij_emitted=0 binary_interaction_db_modified=0 production_ready=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
