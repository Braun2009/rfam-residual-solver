#!/usr/bin/env python3
"""Fail-closed offline topology-envelope screen for the O2/C3H8 PR path."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import itertools
import json
import math
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable


CONTRACT_SCHEMA = "rfam.o2-c3h8-topology-envelope-contract/v1"
RISK_SCHEMA = "rfam.o2-c3h8-topology-model-risk/v1"
TRAJECTORY_SCHEMA = "rfam.hem-inlet-mixing-path-diagnostic/v1"
EVIDENCE_SCHEMA = "rfam.o2-c3h8-topology-envelope-evidence/v1"
FAIL_MARKER = "RFAM_O2_C3H8_TOPOLOGY_ENVELOPE_FAIL"

STATUS_THREE_PHASE = "THREE_PHASE_REQUIRED"
STATUS_OUTSIDE = "OUT_OF_ENVELOPE"
STATUS_CLEAR = "SCREEN_ONLY_NO_MODEL_RISK_DETECTED"

CONTRACT_PATH = (
    "realFluidAllMach/tests/reference/"
    "o2_c3h8_topology_envelope_contract_v1.json"
)
TRAJECTORY_PATH = (
    "realFluidAllMach/tests/reference/"
    "o2_c3h8_predicted_trajectory_512_v1.csv"
)
RISK_PATH = (
    "realFluidAllMach/tests/reference/"
    "o2_c3h8_topology_model_risk_v1.json"
)
TOPOLOGY_REGRESSION_PATH = (
    "realFluidAllMach/tests/test_pr_binary_phase_topology_probe.py"
)

EXPECTED_TRAJECTORY_SHA256 = (
    "a6d498887225a685b2f7b5d40248913ad2e7ba93158e12bbe64963a01783cdc0"
)
EXPECTED_TRAJECTORY_BYTES = 37328
EXPECTED_TRAJECTORY_ROWS = 225
EXPECTED_RISK_SHA256 = (
    "3ef2f001f56af5cbdb8fa851992736ec517cb17c5ed1656f32b85d9946d072b2"
)
EXPECTED_RISK_BYTES = 3361
EXPECTED_TOPOLOGY_REGRESSION_SHA256 = (
    "857dca5af09200aa75dda022e547e8f00dbcb9d7a7de8f24c40a55814f590043"
)
EXPECTED_TOPOLOGY_REGRESSION_BYTES = 16404

EXPECTED_TRAJECTORY_EXECUTABLE_SHA256 = (
    "8e9a5b79c6edaa39f5aa6d01f8453665cdabae9d302b3f258856946eb0f7c4c6"
)
EXPECTED_TOPOLOGY_EXECUTABLE_SHA256 = (
    "e48518f7bc36a24f28ec8b67bd8d83e0858c4bb50143c42a1371b281664a9ddc"
)
EXPECTED_COMPILER_SHA256 = (
    "0c0b987719385f8e242819dfef6a3e97b4b2540ce4035dbc11565df2966bd69e"
)
EXPECTED_TRAJECTORY_BUNDLE_SHA256 = (
    "69edadb5d4885a1bdc262af3e1de7ad81d360810ba43882f44d23dca2851766a"
)
EXPECTED_TOPOLOGY_BUNDLE_SHA256 = (
    "fc8b1e05d9e140b055c209b7b92af21a44db6caaab599800cb830635f535d176"
)

COMMON_SOURCE_PATHS = (
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
    "realFluidAllMach/src/thermo/PhaseEquilibrium/PRMixtureKernel.C",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/PRMixtureKernel.H",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/PRTPFlash.C",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/PRTPFlash.H",
)

TOPOLOGY_SOURCE_PATHS = (
    "realFluidAllMach/tools/pr_binary_phase_topology_probe.cpp",
    *COMMON_SOURCE_PATHS,
)

TRAJECTORY_SOURCE_PATHS = (
    "realFluidAllMach/tools/hem_inlet_mixing_path_probe.cpp",
    *COMMON_SOURCE_PATHS,
    "realFluidAllMach/src/thermo/PhaseEquilibrium/ConservedPhaseMixture.H",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/TabulatedIdealCaloric.C",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/TabulatedIdealCaloric.H",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/PRPhaseCaloric.C",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/PRPhaseCaloric.H",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/PRPHFlash.C",
    "realFluidAllMach/src/thermo/PhaseEquilibrium/PRPHFlash.H",
)

INCLUDE_DIRECTORIES = (
    "realFluidAllMach/src",
    "realFluidAllMach/src/thermo",
    "realFluidAllMach/src/thermo/MixingRules",
    "realFluidAllMach/src/thermo/EquationOfState",
    "realFluidAllMach/src/thermo/RootSelection",
    "realFluidAllMach/src/thermo/PhaseEquilibrium",
)

TRAJECTORY_HEADER = (
    "pressure_pa",
    "kij",
    "fuel_mass_fraction",
    "o2_mole_fraction",
    "target_h_j_kg",
    "status",
    "phase_status",
    "temperature_k",
    "region",
    "vapor_mass_fraction",
    "enthalpy_relative_residual",
    "valid_temperature_samples",
    "invalid_temperature_samples",
    "phase_topology_validated",
    "production_ready",
)

EXPECTED_APPROVALS = {
    "physical_vlle_validated": False,
    "kij_fitted": False,
    "three_phase_owner_implemented": False,
    "hem_enabled": False,
    "fluent_execution_approved": False,
    "hot_start_approved": False,
    "production_ready": False,
}

EXPECTED_UNCERTAINTY = {
    "basis": "declared-screening-intervals-not-physical-uncertainty",
    "kij_offset": [-0.005, 0.005],
    "temperature_k_offset": [-0.2, 0.2],
    "pressure_pa_offset": [-50000.0, 50000.0],
    "overall_z_o2_offset": [-0.002, 0.002],
}

EXPECTED_NUMERICAL_ERROR = {
    "basis": "absolute-bounds-added-by-minkowski-sum",
    "kij_offset": [-1.0e-5, 1.0e-5],
    "temperature_k_offset": [-0.01, 0.01],
    "pressure_pa_offset": [-1000.0, 1000.0],
    "overall_z_o2_offset": [-1.0e-5, 1.0e-5],
}

EXPECTED_SUPPORTED = {
    "pressure_pa": [4.0e6, 6.0e6],
    "temperature_k": [110.0, 300.0],
    "kij": [-0.05, 0.15],
    "overall_z_o2": [0.0, 1.0],
}

EXPECTED_GUARD = {
    "pressure_pa_offset": [-50000.0, 50000.0],
    "temperature_k_offset": [-0.05, 0.05],
    "kij_offset": [-0.001, 0.001],
    "overall_z_o2_offset": [0.0, 0.0],
}

EXPECTED_RISK_RECORDS = (
    ("risk-5mpa-kij-0.055", 5.0e6, 156.37717703, 0.055),
    ("risk-5mpa-kij-0.10", 5.0e6, 155.6846784, 0.10),
    ("risk-5mpa-kij-0.15", 5.0e6, 155.1625280618, 0.15),
)

EXPECTED_LIMITATIONS = [
    "the registry records possible-three-phase model risk, not physical VLLE",
    "overall composition is conservatively unlocalized and therefore spans zero through one",
    "the declared intervals are screening inputs and are not fitted parameter or measurement uncertainty",
    "a clear screen does not prove global phase-topology uniqueness or physical validity",
    "no kij fit, three-phase owner, HEM capability, Fluent execution, hot-start, or production flag is enabled",
]

TRAJECTORY_SUMMARY = re.compile(
    r"# HEM_INLET_MIXING_PATH_DIAGNOSTIC_BLOCKED rows=225 successes=220 "
    r"failures=5 unbracketed=3 two_phase=126 scan_intervals=512 "
    r"phase_topology_validated=0 production_ready=0"
)


class EnvelopeError(ValueError):
    """A malformed authority or unsupported envelope fails closed."""


@dataclass(frozen=True)
class Interval:
    low: float
    high: float

    def shifted(self, offset: "Interval") -> "Interval":
        return Interval(self.low + offset.low, self.high + offset.high)

    def hull(self, other: "Interval") -> "Interval":
        return Interval(min(self.low, other.low), max(self.high, other.high))

    def overlaps(self, other: "Interval") -> bool:
        return self.low <= other.high and other.low <= self.high

    def contained_by(self, other: "Interval") -> bool:
        return other.low <= self.low and self.high <= other.high

    def samples(self) -> tuple[float, ...]:
        if self.low == self.high:
            return (self.low,)
        return (self.low, 0.5 * (self.low + self.high), self.high)


@dataclass(frozen=True)
class Anchor:
    identifier: str
    group: str
    order: float
    pressure: float
    temperature: float
    kij: float
    z_o2: float
    committable: bool = True


@dataclass(frozen=True)
class RiskRecord:
    identifier: str
    pressure: float
    temperature: float
    kij: float
    z_interval: Interval


Box = Dict[str, Interval]


def _reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise EnvelopeError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def _read_bytes(path: Path, label: str) -> bytes:
    try:
        return path.read_bytes()
    except OSError as exc:
        raise EnvelopeError(f"{label}-read-failed:{exc}") from exc


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _load_json(path: Path, label: str) -> tuple[dict[str, Any], bytes]:
    payload = _read_bytes(path, label)
    if b"\r" in payload:
        raise EnvelopeError(f"{label}-not-canonical-lf")
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as exc:
        raise EnvelopeError(f"{label}-not-ascii") from exc
    try:
        value = json.loads(text, object_pairs_hook=_reject_duplicates)
    except EnvelopeError:
        raise
    except (json.JSONDecodeError, ValueError) as exc:
        raise EnvelopeError(f"{label}-invalid-json:{exc}") from exc
    if not isinstance(value, dict):
        raise EnvelopeError(f"{label}-root-not-object")
    return value, payload


def _object(value: Any, label: str, keys: Iterable[str]) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise EnvelopeError(f"type-mismatch:{label}:object")
    expected = set(keys)
    actual = set(value)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise EnvelopeError(f"key-mismatch:{label}:missing={missing}:extra={extra}")
    return value


def _exact(actual: Any, expected: Any, label: str) -> None:
    if type(actual) is not type(expected) or actual != expected:
        raise EnvelopeError(f"value-mismatch:{label}")


def _number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise EnvelopeError(f"type-mismatch:{label}:number")
    result = float(value)
    if not math.isfinite(result):
        raise EnvelopeError(f"nonfinite:{label}")
    return result


def _interval(value: Any, label: str, contains_zero: bool = False) -> Interval:
    if not isinstance(value, list) or len(value) != 2:
        raise EnvelopeError(f"type-mismatch:{label}:two-number-array")
    low = _number(value[0], f"{label}[0]")
    high = _number(value[1], f"{label}[1]")
    if low > high:
        raise EnvelopeError(f"reversed-interval:{label}")
    if contains_zero and not (low <= 0.0 <= high):
        raise EnvelopeError(f"offset-interval-must-contain-zero:{label}")
    return Interval(low, high)


def _workspace_path(root: Path, relative: str, label: str) -> Path:
    if not isinstance(relative, str) or not relative or "\\" in relative:
        raise EnvelopeError(f"invalid-relative-path:{label}")
    candidate = Path(relative)
    if candidate.is_absolute() or ".." in candidate.parts:
        raise EnvelopeError(f"path-escapes-root:{label}")
    resolved = (root / candidate).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise EnvelopeError(f"path-escapes-root:{label}") from exc
    return resolved


def source_bundle_sha(root: Path, paths: tuple[str, ...]) -> str:
    digest = hashlib.sha256()
    for relative in paths:
        payload = _read_bytes(_workspace_path(root, relative, relative), relative)
        digest.update(relative.encode("ascii"))
        digest.update(b"\0")
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return digest.hexdigest()


def _validate_build(
    value: Any,
    label: str,
    expected_sources: tuple[str, ...],
    expected_bundle_sha: str,
    expected_executable_sha: str,
) -> None:
    build = _object(
        value,
        label,
        {
            "platform",
            "compiler_path",
            "compiler_sha256",
            "compiler_version",
            "language",
            "optimization",
            "openmp",
            "include_directories",
            "source_bundle_algorithm",
            "source_files",
            "source_bundle_sha256",
            "expected_executable_sha256",
        },
    )
    _exact(build["platform"], "WSL-Ubuntu-20.04-x86_64", f"{label}.platform")
    _exact(build["compiler_path"], "/usr/bin/x86_64-linux-gnu-g++-9", f"{label}.compiler_path")
    _exact(build["compiler_sha256"], EXPECTED_COMPILER_SHA256, f"{label}.compiler_sha256")
    _exact(
        build["compiler_version"],
        "x86_64-linux-gnu-g++-9 (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0",
        f"{label}.compiler_version",
    )
    _exact(build["language"], "c++17", f"{label}.language")
    _exact(build["optimization"], "-O2", f"{label}.optimization")
    _exact(build["openmp"], "-fopenmp", f"{label}.openmp")
    _exact(build["include_directories"], list(INCLUDE_DIRECTORIES), f"{label}.include_directories")
    _exact(
        build["source_bundle_algorithm"],
        "sha256(path-nul-size-u64be-bytes)",
        f"{label}.source_bundle_algorithm",
    )
    _exact(build["source_files"], list(expected_sources), f"{label}.source_files")
    _exact(build["source_bundle_sha256"], expected_bundle_sha, f"{label}.source_bundle_sha256")
    _exact(build["expected_executable_sha256"], expected_executable_sha, f"{label}.expected_executable_sha256")


def validate_risk_registry(root: Path, path: Path) -> tuple[list[RiskRecord], bytes]:
    document, payload = _load_json(path, "risk-registry")
    root_object = _object(
        document,
        "risk-registry",
        {
            "schema",
            "status",
            "pair",
            "classification",
            "backend",
            "approach_guard_offsets",
            "records",
            "approvals",
            "limitations",
        },
    )
    _exact(root_object["schema"], RISK_SCHEMA, "risk-registry.schema")
    _exact(root_object["status"], "possible-three-phase-model-risk-only", "risk-registry.status")
    _exact(root_object["pair"], ["O2", "C3H8"], "risk-registry.pair")
    _exact(
        root_object["classification"],
        {
            "name": "possible-three-phase-model-risk",
            "physical_vlle_validated": False,
            "overall_z_scope": "conservative-zero-through-one-until-validated-localization",
        },
        "risk-registry.classification",
    )
    _exact(
        root_object["backend"],
        {
            "probe_source_path": "realFluidAllMach/tools/pr_binary_phase_topology_probe.cpp",
            "probe_source_sha256": "a0e9f8499775dc32f597ca455d8580ec397b5484240c22a8651dd2321cb8b746",
            "regression_path": TOPOLOGY_REGRESSION_PATH,
            "regression_sha256": EXPECTED_TOPOLOGY_REGRESSION_SHA256,
            "regression_byte_count": EXPECTED_TOPOLOGY_REGRESSION_BYTES,
            "pr_tp_source_path": "realFluidAllMach/src/thermo/PhaseEquilibrium/PRTPFlash.C",
            "pr_tp_source_sha256": "e16a747f84c6d6891c3154bb385cad8ce39c6233952d8b8e1ce4fee0e07130dc",
            "composition_intervals": 20000,
            "expected_executable_sha256": EXPECTED_TOPOLOGY_EXECUTABLE_SHA256,
        },
        "risk-registry.backend",
    )
    regression_path = _workspace_path(
        root, TOPOLOGY_REGRESSION_PATH, "topology-regression"
    )
    regression_payload = _read_bytes(regression_path, "topology-regression")
    if len(regression_payload) != EXPECTED_TOPOLOGY_REGRESSION_BYTES:
        raise EnvelopeError("topology-regression-byte-count-mismatch")
    if _sha256(regression_payload) != EXPECTED_TOPOLOGY_REGRESSION_SHA256:
        raise EnvelopeError("topology-regression-sha256-mismatch")
    _exact(root_object["approach_guard_offsets"], EXPECTED_GUARD, "risk-registry.approach_guard_offsets")
    _exact(root_object["approvals"], EXPECTED_APPROVALS, "risk-registry.approvals")
    _exact(root_object["limitations"], EXPECTED_LIMITATIONS, "risk-registry.limitations")

    records_value = root_object["records"]
    if not isinstance(records_value, list) or len(records_value) != len(EXPECTED_RISK_RECORDS):
        raise EnvelopeError("risk-record-count-mismatch")
    records: list[RiskRecord] = []
    for index, expected in enumerate(EXPECTED_RISK_RECORDS):
        item = _object(
            records_value[index],
            f"risk-registry.records[{index}]",
            {
                "id",
                "pressure_pa",
                "temperature_k",
                "kij",
                "overall_z_o2_interval",
                "tie_segments",
                "three_phase_tangent_pairs",
                "three_phase_tie_segments",
                "possible_three_phase",
                "physical_vlle_validated",
                "production_ready",
            },
        )
        _exact(item["id"], expected[0], f"risk-registry.records[{index}].id")
        _exact(_number(item["pressure_pa"], "pressure"), expected[1], f"risk-registry.records[{index}].pressure")
        _exact(_number(item["temperature_k"], "temperature"), expected[2], f"risk-registry.records[{index}].temperature")
        _exact(_number(item["kij"], "kij"), expected[3], f"risk-registry.records[{index}].kij")
        z_interval = _interval(item["overall_z_o2_interval"], "risk-z")
        _exact([z_interval.low, z_interval.high], [0.0, 1.0], f"risk-registry.records[{index}].overall_z_o2_interval")
        _exact(item["tie_segments"], 2, f"risk-registry.records[{index}].tie_segments")
        _exact(item["three_phase_tangent_pairs"], 1, f"risk-registry.records[{index}].three_phase_tangent_pairs")
        _exact(item["three_phase_tie_segments"], 2, f"risk-registry.records[{index}].three_phase_tie_segments")
        _exact(item["possible_three_phase"], True, f"risk-registry.records[{index}].possible_three_phase")
        _exact(item["physical_vlle_validated"], False, f"risk-registry.records[{index}].physical_vlle_validated")
        _exact(item["production_ready"], False, f"risk-registry.records[{index}].production_ready")
        records.append(RiskRecord(expected[0], expected[1], expected[2], expected[3], z_interval))
    return records, payload


def validate_contract(root: Path, path: Path) -> tuple[dict[str, Any], list[RiskRecord], bytes, bytes]:
    expected_path = _workspace_path(root, CONTRACT_PATH, "contract")
    if path.resolve() != expected_path:
        raise EnvelopeError("contract-path-mismatch")
    document, _ = _load_json(path, "contract")
    root_object = _object(
        document,
        "contract",
        {
            "schema",
            "contract_version",
            "status",
            "scope",
            "authorities",
            "trajectory",
            "uncertainty_intervals",
            "numerical_interpolation_error_intervals",
            "sampling",
            "supported_screen_envelope",
            "risk_policy",
            "approvals",
            "limitations",
        },
    )
    _exact(root_object["schema"], CONTRACT_SCHEMA, "contract.schema")
    _exact(root_object["contract_version"], 1, "contract.contract_version")
    _exact(root_object["status"], "offline-model-risk-gate-blocked", "contract.status")
    _exact(
        root_object["scope"],
        "current-rfam-pr-o2-c3h8-predicted-trajectory-uncertainty-tube-no-physical-vlle",
        "contract.scope",
    )
    _exact(root_object["approvals"], EXPECTED_APPROVALS, "contract.approvals")
    _exact(root_object["limitations"], EXPECTED_LIMITATIONS, "contract.limitations")
    _exact(root_object["uncertainty_intervals"], EXPECTED_UNCERTAINTY, "contract.uncertainty_intervals")
    _exact(
        root_object["numerical_interpolation_error_intervals"],
        EXPECTED_NUMERICAL_ERROR,
        "contract.numerical_interpolation_error_intervals",
    )
    _exact(root_object["supported_screen_envelope"], EXPECTED_SUPPORTED, "contract.supported_screen_envelope")
    _exact(
        root_object["sampling"],
        {
            "group_key": ["pressure_pa", "kij"],
            "order_key": "fuel_mass_fraction",
            "segment_bound": "axis-aligned-hull-of-adjacent-effective-anchor-intervals",
            "sample_points_per_axis": ["lower", "midpoint", "upper"],
            "closed_interval_intersection": True,
            "uncertainty_widening_must_be_monotone": True,
        },
        "contract.sampling",
    )
    _exact(
        root_object["trajectory"],
        {
            "artifact_id": "predicted-trajectory",
            "schema": TRAJECTORY_SCHEMA,
            "scan_intervals": 512,
            "rows": EXPECTED_TRAJECTORY_ROWS,
            "uncommittable_status_policy": STATUS_OUTSIDE,
        },
        "contract.trajectory",
    )
    _exact(
        root_object["risk_policy"],
        {
            "registry_artifact_id": "model-risk-registry",
            "intersection_or_approach_status": STATUS_THREE_PHASE,
            "outside_or_uncommittable_status": STATUS_OUTSIDE,
            "clear_status": STATUS_CLEAR,
            "status_precedence": [STATUS_THREE_PHASE, STATUS_OUTSIDE, STATUS_CLEAR],
            "physical_vlle_claim": False,
            "z_localization_policy": "conservative-all-compositions-until-validated",
        },
        "contract.risk_policy",
    )

    authorities = _object(
        root_object["authorities"],
        "contract.authorities",
        {
            "trajectory_reference",
            "model_risk_registry",
            "trajectory_probe_build",
            "topology_probe_build",
        },
    )
    trajectory_authority = _object(
        authorities["trajectory_reference"],
        "contract.authorities.trajectory_reference",
        {"id", "path", "sha256", "byte_count"},
    )
    _exact(trajectory_authority["id"], "predicted-trajectory", "trajectory-authority.id")
    _exact(trajectory_authority["path"], TRAJECTORY_PATH, "trajectory-authority.path")
    _exact(trajectory_authority["sha256"], EXPECTED_TRAJECTORY_SHA256, "trajectory-authority.sha256")
    _exact(trajectory_authority["byte_count"], EXPECTED_TRAJECTORY_BYTES, "trajectory-authority.byte_count")
    trajectory_path = _workspace_path(root, TRAJECTORY_PATH, "trajectory-reference")
    trajectory_payload = _read_bytes(trajectory_path, "trajectory-reference")
    if len(trajectory_payload) != EXPECTED_TRAJECTORY_BYTES:
        raise EnvelopeError("trajectory-reference-byte-count-mismatch")
    if _sha256(trajectory_payload) != EXPECTED_TRAJECTORY_SHA256:
        raise EnvelopeError("trajectory-reference-sha256-mismatch")

    risk_authority = _object(
        authorities["model_risk_registry"],
        "contract.authorities.model_risk_registry",
        {"id", "path", "sha256", "byte_count"},
    )
    _exact(risk_authority["id"], "model-risk-registry", "risk-authority.id")
    _exact(risk_authority["path"], RISK_PATH, "risk-authority.path")
    _exact(risk_authority["sha256"], EXPECTED_RISK_SHA256, "risk-authority.sha256")
    _exact(risk_authority["byte_count"], EXPECTED_RISK_BYTES, "risk-authority.byte_count")
    risk_path = _workspace_path(root, RISK_PATH, "risk-registry")
    if _sha256(_read_bytes(risk_path, "risk-registry")) != EXPECTED_RISK_SHA256:
        raise EnvelopeError("risk-registry-sha256-mismatch")
    risks, risk_payload = validate_risk_registry(root, risk_path)

    _validate_build(
        authorities["trajectory_probe_build"],
        "trajectory_probe_build",
        TRAJECTORY_SOURCE_PATHS,
        EXPECTED_TRAJECTORY_BUNDLE_SHA256,
        EXPECTED_TRAJECTORY_EXECUTABLE_SHA256,
    )
    _validate_build(
        authorities["topology_probe_build"],
        "topology_probe_build",
        TOPOLOGY_SOURCE_PATHS,
        EXPECTED_TOPOLOGY_BUNDLE_SHA256,
        EXPECTED_TOPOLOGY_EXECUTABLE_SHA256,
    )
    if source_bundle_sha(root, TRAJECTORY_SOURCE_PATHS) != EXPECTED_TRAJECTORY_BUNDLE_SHA256:
        raise EnvelopeError("trajectory-source-bundle-mismatch")
    if source_bundle_sha(root, TOPOLOGY_SOURCE_PATHS) != EXPECTED_TOPOLOGY_BUNDLE_SHA256:
        raise EnvelopeError("topology-source-bundle-mismatch")
    return document, risks, trajectory_payload, risk_payload


def parse_trajectory(payload: bytes) -> list[Anchor]:
    if b"\r" in payload:
        raise EnvelopeError("trajectory-reference-not-canonical-lf")
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as exc:
        raise EnvelopeError("trajectory-reference-not-ascii") from exc
    lines = text.splitlines()
    if len(lines) != EXPECTED_TRAJECTORY_ROWS + 3:
        raise EnvelopeError(f"trajectory-line-count:{len(lines)}")
    if lines[0] != f"# {TRAJECTORY_SCHEMA}":
        raise EnvelopeError("trajectory-schema-marker-mismatch")
    if TRAJECTORY_SUMMARY.fullmatch(lines[-1]) is None:
        raise EnvelopeError("trajectory-summary-marker-mismatch")
    reader = csv.DictReader(io.StringIO("\n".join(lines[1:-1])))
    if tuple(reader.fieldnames or ()) != TRAJECTORY_HEADER:
        raise EnvelopeError("trajectory-header-mismatch")
    rows = list(reader)
    if len(rows) != EXPECTED_TRAJECTORY_ROWS:
        raise EnvelopeError("trajectory-row-count-mismatch")

    anchors: list[Anchor] = []
    status_counts: dict[str, int] = {}
    two_phase_successes = 0
    for index, row in enumerate(rows):
        try:
            pressure = float(row["pressure_pa"])
            kij = float(row["kij"])
            order = float(row["fuel_mass_fraction"])
            z_o2 = float(row["o2_mole_fraction"])
            temperature = float(row["temperature_k"])
        except (KeyError, ValueError) as exc:
            raise EnvelopeError(f"trajectory-number:row-{index + 1}") from exc
        if not all(math.isfinite(value) for value in (pressure, kij, order, z_o2, temperature)):
            raise EnvelopeError(f"trajectory-nonfinite:row-{index + 1}")
        if not (0.0 <= order <= 1.0 and 0.0 <= z_o2 <= 1.0):
            raise EnvelopeError(f"trajectory-composition-bound:row-{index + 1}")
        status = row["status"]
        if status not in {"SUCCESS", "NO_ENTHALPY_BRACKET", "ENTHALPY_RESIDUAL_FAILURE"}:
            raise EnvelopeError(f"trajectory-status:row-{index + 1}")
        if row["phase_topology_validated"] != "0" or row["production_ready"] != "0":
            raise EnvelopeError(f"trajectory-approval-marker:row-{index + 1}")
        status_counts[status] = status_counts.get(status, 0) + 1
        if status == "SUCCESS" and row["region"] == "TWO_PHASE":
            two_phase_successes += 1
        group = f"p={row['pressure_pa']}|kij={row['kij']}"
        anchors.append(
            Anchor(
                identifier=f"row-{index + 1}",
                group=group,
                order=order,
                pressure=pressure,
                temperature=temperature,
                kij=kij,
                z_o2=z_o2,
                committable=status == "SUCCESS",
            )
        )
    if status_counts != {
        "SUCCESS": 220,
        "NO_ENTHALPY_BRACKET": 3,
        "ENTHALPY_RESIDUAL_FAILURE": 2,
    }:
        raise EnvelopeError(f"trajectory-status-counts:{status_counts}")
    if two_phase_successes != 126:
        raise EnvelopeError(f"trajectory-two-phase-count:{two_phase_successes}")
    groups: dict[str, list[Anchor]] = {}
    for anchor in anchors:
        groups.setdefault(anchor.group, []).append(anchor)
    if len(groups) != 15 or any(len(group) != 15 for group in groups.values()):
        raise EnvelopeError("trajectory-group-shape-mismatch")
    for group in groups.values():
        orders = [anchor.order for anchor in group]
        if orders != sorted(orders) or len(set(orders)) != len(orders):
            raise EnvelopeError("trajectory-group-order-mismatch")
    return anchors


def offset_set(value: dict[str, Any]) -> dict[str, Interval]:
    return {
        "pressure_pa": _interval(value["pressure_pa_offset"], "pressure-offset", True),
        "temperature_k": _interval(value["temperature_k_offset"], "temperature-offset", True),
        "kij": _interval(value["kij_offset"], "kij-offset", True),
        "overall_z_o2": _interval(value["overall_z_o2_offset"], "z-offset", True),
    }


def supported_set(value: dict[str, Any]) -> dict[str, Interval]:
    return {
        "pressure_pa": _interval(value["pressure_pa"], "supported-pressure"),
        "temperature_k": _interval(value["temperature_k"], "supported-temperature"),
        "kij": _interval(value["kij"], "supported-kij"),
        "overall_z_o2": _interval(value["overall_z_o2"], "supported-z"),
    }


def anchor_box(
    anchor: Anchor,
    uncertainty: dict[str, Interval],
    numerical_error: dict[str, Interval],
) -> Box:
    nominal = {
        "pressure_pa": Interval(anchor.pressure, anchor.pressure),
        "temperature_k": Interval(anchor.temperature, anchor.temperature),
        "kij": Interval(anchor.kij, anchor.kij),
        "overall_z_o2": Interval(anchor.z_o2, anchor.z_o2),
    }
    return {
        name: nominal[name].shifted(uncertainty[name]).shifted(numerical_error[name])
        for name in nominal
    }


def hull_box(left: Box, right: Box) -> Box:
    return {name: left[name].hull(right[name]) for name in left}


def risk_box(record: RiskRecord, guard: dict[str, Interval]) -> Box:
    nominal = {
        "pressure_pa": Interval(record.pressure, record.pressure),
        "temperature_k": Interval(record.temperature, record.temperature),
        "kij": Interval(record.kij, record.kij),
        "overall_z_o2": record.z_interval,
    }
    return {name: nominal[name].shifted(guard[name]) for name in nominal}


def boxes_overlap(left: Box, right: Box) -> bool:
    return all(left[name].overlaps(right[name]) for name in left)


def box_contained(box: Box, supported: dict[str, Interval]) -> bool:
    return all(box[name].contained_by(supported[name]) for name in box)


def point_in_box(point: dict[str, float], box: Box) -> bool:
    return all(box[name].low <= point[name] <= box[name].high for name in box)


def sample_box(box: Box) -> Iterable[dict[str, float]]:
    names = tuple(box)
    axes = [box[name].samples() for name in names]
    for values in itertools.product(*axes):
        yield dict(zip(names, values))


def evaluate_tube(
    anchors: list[Anchor],
    risks: list[RiskRecord],
    uncertainty: dict[str, Interval],
    numerical_error: dict[str, Interval],
    guard: dict[str, Interval],
    supported: dict[str, Interval],
) -> dict[str, Any]:
    if not anchors:
        raise EnvelopeError("empty-trajectory")
    anchor_boxes = [(anchor.identifier, anchor_box(anchor, uncertainty, numerical_error)) for anchor in anchors]
    groups: dict[str, list[tuple[Anchor, Box]]] = {}
    for anchor, (_, box) in zip(anchors, anchor_boxes):
        groups.setdefault(anchor.group, []).append((anchor, box))
    segment_boxes: list[tuple[str, Box]] = []
    for group_name, values in groups.items():
        ordered = sorted(values, key=lambda item: item[0].order)
        for index in range(len(ordered) - 1):
            left, right = ordered[index], ordered[index + 1]
            segment_boxes.append(
                (
                    f"{group_name}|segment={index}",
                    hull_box(left[1], right[1]),
                )
            )
    all_boxes = [("anchor", *entry) for entry in anchor_boxes]
    all_boxes.extend(("segment", *entry) for entry in segment_boxes)
    risk_boxes = [(record.identifier, risk_box(record, guard)) for record in risks]

    risk_hits: list[str] = []
    sampled_risk_hits = 0
    sample_count = 0
    outside: list[str] = []
    for kind, identifier, box in all_boxes:
        if not box_contained(box, supported):
            outside.append(f"{kind}:{identifier}")
        for risk_id, candidate in risk_boxes:
            if boxes_overlap(box, candidate):
                risk_hits.append(f"{kind}:{identifier}->{risk_id}")
        for point in sample_box(box):
            sample_count += 1
            if any(point_in_box(point, candidate) for _, candidate in risk_boxes):
                sampled_risk_hits += 1

    uncommittable = [anchor.identifier for anchor in anchors if not anchor.committable]
    if risk_hits:
        status = STATUS_THREE_PHASE
    elif outside or uncommittable:
        status = STATUS_OUTSIDE
    else:
        status = STATUS_CLEAR
    return {
        "schema": EVIDENCE_SCHEMA,
        "status": status,
        "anchors": len(anchor_boxes),
        "segments": len(segment_boxes),
        "sample_points": sample_count,
        "model_risk_records": len(risk_boxes),
        "risk_box_hits": len(risk_hits),
        "sampled_risk_hits": sampled_risk_hits,
        "out_of_envelope_boxes": len(outside),
        "uncommittable_anchors": len(uncommittable),
        "first_risk_hits": risk_hits[:8],
        "first_out_of_envelope": outside[:8],
        "physical_vlle_validated": False,
        "kij_fitted": False,
        "three_phase_owner_implemented": False,
        "production_ready": False,
    }


def _verify_executable(path: Path, expected_sha: str, label: str) -> None:
    if _sha256(_read_bytes(path, label)) != expected_sha:
        raise EnvelopeError(f"{label}-sha256-mismatch")


def verify_trajectory_probe(path: Path, reference: bytes) -> None:
    _verify_executable(path, EXPECTED_TRAJECTORY_EXECUTABLE_SHA256, "trajectory-probe")
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    try:
        result = subprocess.run(
            [str(path), "--scan-intervals", "512"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=180,
            env=environment,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise EnvelopeError(f"trajectory-probe-execution:{exc}") from exc
    if result.returncode != 0 or result.stderr:
        raise EnvelopeError(f"trajectory-probe-returned:{result.returncode}")
    if result.stdout != reference:
        raise EnvelopeError("trajectory-probe-reference-byte-mismatch")


def verify_topology_probe(path: Path, risks: list[RiskRecord]) -> None:
    _verify_executable(path, EXPECTED_TOPOLOGY_EXECUTABLE_SHA256, "topology-probe")
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    for record in risks:
        try:
            result = subprocess.run(
                [
                    str(path),
                    "--state",
                    repr(record.pressure),
                    repr(record.temperature),
                    repr(record.kij),
                    "--composition-intervals",
                    "20000",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
                timeout=180,
                env=environment,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise EnvelopeError(f"topology-probe-execution:{record.identifier}:{exc}") from exc
        if result.returncode != 0 or result.stderr:
            raise EnvelopeError(f"topology-probe-returned:{record.identifier}:{result.returncode}")
        try:
            lines = result.stdout.decode("ascii").splitlines()
        except UnicodeDecodeError as exc:
            raise EnvelopeError(f"topology-probe-not-ascii:{record.identifier}") from exc
        if len(lines) != 4 or lines[0] != "# rfam.pr-binary-phase-topology-screen/v1":
            raise EnvelopeError(f"topology-probe-schema:{record.identifier}")
        rows = list(csv.DictReader(io.StringIO("\n".join(lines[1:3]))))
        if len(rows) != 1:
            raise EnvelopeError(f"topology-probe-row-count:{record.identifier}")
        row = rows[0]
        if not (
            row.get("possible_three_phase") == "1"
            and row.get("tie_segments") == "2"
            and row.get("three_phase_tangent_pairs") == "1"
            and row.get("three_phase_tie_segments") == "2"
            and row.get("production_ready") == "0"
        ):
            raise EnvelopeError(f"topology-probe-classification:{record.identifier}")
        if "possible_three_phase_rows=1" not in lines[-1]:
            raise EnvelopeError(f"topology-probe-summary:{record.identifier}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--trajectory-probe", type=Path)
    parser.add_argument("--topology-probe", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        root = args.root.resolve()
        document, risks, trajectory_payload, _ = validate_contract(
            root, args.contract.resolve()
        )
        anchors = parse_trajectory(trajectory_payload)
        uncertainty = offset_set(document["uncertainty_intervals"])
        numerical_error = offset_set(
            document["numerical_interpolation_error_intervals"]
        )
        guard = offset_set(EXPECTED_GUARD)
        supported = supported_set(document["supported_screen_envelope"])
        evidence = evaluate_tube(
            anchors, risks, uncertainty, numerical_error, guard, supported
        )
        evidence["trajectory_reference_sha256"] = EXPECTED_TRAJECTORY_SHA256
        evidence["risk_registry_sha256"] = EXPECTED_RISK_SHA256
        evidence["trajectory_probe_verified"] = args.trajectory_probe is not None
        evidence["topology_probe_verified"] = args.topology_probe is not None
        if args.trajectory_probe is not None:
            verify_trajectory_probe(args.trajectory_probe.resolve(), trajectory_payload)
        if args.topology_probe is not None:
            verify_topology_probe(args.topology_probe.resolve(), risks)
    except Exception as exc:  # Any malformed or unexpected state fails closed.
        reason = str(exc) if isinstance(exc, EnvelopeError) else f"unexpected:{type(exc).__name__}:{exc}"
        reason = reason.replace("\r", " ").replace("\n", " ")
        print(f"{FAIL_MARKER} reason={reason}")
        return 2

    print(json.dumps(evidence, sort_keys=True, separators=(",", ":")))
    print(
        f"RFAM_O2_C3H8_TOPOLOGY_ENVELOPE_{evidence['status']} "
        f"anchors={evidence['anchors']} segments={evidence['segments']} "
        f"risk_box_hits={evidence['risk_box_hits']} "
        f"out_of_envelope_boxes={evidence['out_of_envelope_boxes']} "
        "physical_vlle_validated=0 kij_fitted=0 "
        "three_phase_owner_implemented=0 production_ready=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
