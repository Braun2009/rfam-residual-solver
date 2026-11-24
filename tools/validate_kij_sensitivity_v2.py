#!/usr/bin/env python3
"""Validate the current-source O2/C3H8 PR kij sensitivity contract."""

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


SCHEMA = "rfam.kij-sensitivity-contract/v2"
CSV_SCHEMA = "rfam.kij-sensitivity-probe/v1"
SUCCESS_MARKER = "RFAM_KIJ_SENSITIVITY_V2_SCREEN_PASS"
BLOCKED_MARKER = "RFAM_KIJ_SENSITIVITY_V2_CONTRACT_BLOCKED"
FAIL_MARKER = "RFAM_KIJ_SENSITIVITY_V2_FAIL"

EXPECTED_SOURCE_BUNDLE_SHA256 = (
    "0da258a828d59604dfb98191753788a88bd0bdffb37fcb372e721e2b6a98a53b"
)
EXPECTED_REFERENCE_SHA256 = (
    "535c970495126eff832ddf73b36102d7725415563f8c9f3a94ec557fdaa6b671"
)
EXPECTED_EXECUTABLE_SHA256 = (
    "2b14dc31c2cbbab6b471ad4c972cbf6ba50d66692f052e9d5e6badb2c5d84348"
)

EXPECTED_ROOT_KEYS = {
    "schema", "contract_version", "status", "scope", "approvals",
    "upstream", "pair", "model", "state_grid", "outputs", "tolerances",
    "core_provenance", "source_bundle", "reference", "build", "conclusion",
    "limitations",
}

EXPECTED_APPROVALS = {
    "kij_value_validated": False,
    "production_ready": False,
    "hem_enabled": False,
    "fluent_execution_approved": False,
    "hot_start_approved": False,
}

EXPECTED_UPSTREAM = {
    "wp0_contract": (
        "realFluidAllMach/tests/reference/vortex_physics_validation_contract_v1.json",
        "34b15eae9c89d3c44638748e3283fed59429912f16faa33a224f70f78710439a",
    ),
    "wp0_owner_report": (
        "tmp/a408-wp0-wp2-replacement-v4-20260722.md",
        "05ddf8a971cc79ca76196505ead2be37258cf7797a5f404eb14ec72fd5067420",
    ),
    "wp0_independent_review": (
        "tmp/a408-wp0-wp2-replacement-v4-independent-review-20260722.md",
        "b05980810f2db2eaddfa48839479e9cf5cf2abd80940f700e9f1db986541d626",
    ),
    "tp_ph_repair_owner_report": (
        "tmp/rfam-tp-ph-stability-repair-replacement-v2-20260723.md",
        "a1392524c04615e6f2b13172d1916eab49c943aebe5c54baf1bef47c820e02ee",
    ),
    "tp_ph_repair_independent_review": (
        "tmp/rfam-tp-ph-stability-repair-replacement-v2-independent-review-20260723.md",
        "988643161a70ded8b8519325a33b5a2640077d74b96e3c1afe0db5b86305af56",
    ),
    "caloric_80k_diagnosis": (
        "tmp/rfam-prph-three-phase-readonly-diagnosis-v1-20260723.md",
        "19a6c536616bd53b6b465fcd8aa985aea2a522f13c321f43598c646ffb85f934",
    ),
}

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

KIJ_VALUES = (-0.05, 0.0, 0.055, 0.10, 0.15)
STATE_GRID = (
    ("o2_rich_low_2mpa", 2.0e6, 180.0, 0.90, "pair-sensitivity", "TWO_PHASE"),
    ("global_low_2mpa", 2.0e6, 220.0, 2.8 / 3.85, "pair-sensitivity", "TWO_PHASE"),
    ("fuel_rich_low_2mpa", 2.0e6, 260.0, 0.20, "pair-sensitivity", "TWO_PHASE"),
    ("global_transition_4mpa", 4.0e6, 300.0, 2.8 / 3.85, "pair-sensitivity", "SINGLE_VAPOR"),
    ("global_inlet_mix_4mpa", 4.0e6, 450.0, 2.8 / 3.85, "pair-sensitivity", "SINGLE_VAPOR"),
    ("global_background_4mpa", 4.0e6, 700.0, 2.8 / 3.85, "pair-sensitivity", "SINGLE_VAPOR"),
    ("global_ignition_4mpa", 4.0e6, 1200.0, 2.8 / 3.85, "pair-sensitivity", "SINGLE_VAPOR"),
    ("global_pilot_4mpa", 4.0e6, 2400.0, 2.8 / 3.85, "pair-sensitivity", "SINGLE_VAPOR"),
    ("o2_supercritical_control_6mpa", 6.0e6, 200.0, 1.0, "phase-semantics-control", "SINGLE_SUPERCRITICAL"),
    ("global_dense_6mpa", 6.0e6, 220.0, 2.8 / 3.85, "pair-sensitivity", "TWO_PHASE"),
)

FIELDS = (
    "case_id", "kij", "pressure_pa", "temperature_k", "y_o2", "y_c3h8",
    "status", "region", "quality_defined", "vapor_mass_fraction",
    "vapor_volume_fraction", "rho_kg_m3", "e_j_kg", "liquid_x_o2",
    "vapor_x_o2", "vapor_y_o2", "vapor_y_c3h8",
    "vapor_rho_o2_kg_m3", "vapor_rho_c3h8_kg_m3",
    "single_rho_o2_kg_m3", "single_rho_c3h8_kg_m3",
    "material_residual", "energy_identity_residual",
)

SOURCE_PATHS = (
    "realFluidAllMach/tools/kij_sensitivity_probe.cpp",
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

EXPECTED_LIMITATIONS = [
    "screening states are conservative anchors, not a completed Vortex trajectory envelope",
    "no VLE/PVT data source, fit, holdout, uncertainty, or validated kij is supplied",
    "the historical kij=0.055 value remains legacy unprovenanced diagnostic input",
    "CoolProp and REFPROP parameters are not PR kij values and are not copied here",
    "current TP/PH repair requires a separately bound external topology-envelope admission",
    "the 80 K C3H8 ideal-caloric extension is provenance and does not validate kij",
    "no HEM, Fluent, zero-advance, hot-start, or production capability is enabled",
]


class ContractError(ValueError):
    """A fail-closed contract or reference violation."""


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ContractError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def _decode_json(path: Path) -> Any:
    try:
        text = path.read_text(encoding="ascii")
        return json.loads(text, object_pairs_hook=_reject_duplicate_pairs)
    except ContractError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ContractError(f"contract-json-invalid:{exc}") from exc


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _read_bytes(path: Path, label: str) -> bytes:
    try:
        return path.read_bytes()
    except OSError as exc:
        raise ContractError(f"file-read-failed:{label}:{exc}") from exc


def _valid_sha(value: Any, label: str) -> str:
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ContractError(f"invalid-sha256:{label}")
    if value == "0" * 64:
        raise ContractError(f"invalid-sha256:{label}")
    return value


def _object(value: Any, label: str, keys: set[str]) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ContractError(f"not-object:{label}")
    missing = keys - set(value)
    unknown = set(value) - keys
    if missing:
        raise ContractError(f"missing-keys:{label}:" + ",".join(sorted(missing)))
    if unknown:
        raise ContractError(f"unknown-keys:{label}:" + ",".join(sorted(unknown)))
    return value


def _exact(actual: Any, expected: Any, label: str) -> None:
    if type(actual) is not type(expected) or actual != expected:
        raise ContractError(f"contract-mismatch:{label}")


def _number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"invalid-number:{label}")
    result = float(value)
    if not math.isfinite(result):
        raise ContractError(f"nonfinite-number:{label}")
    return result


def _workspace_path(root: Path, relative: Any, label: str) -> Path:
    if not isinstance(relative, str) or not relative or "\\" in relative:
        raise ContractError(f"invalid-path:{label}")
    if relative.startswith("/") or re.match(r"^[A-Za-z]:", relative):
        raise ContractError(f"path-not-relative:{label}")
    root = root.resolve()
    candidate = (root / relative).resolve()
    if candidate != root and root not in candidate.parents:
        raise ContractError(f"path-escapes-root:{label}")
    return candidate


def source_bundle_sha(root: Path, paths: tuple[str, ...] = SOURCE_PATHS) -> str:
    digest = hashlib.sha256()
    for relative in paths:
        path = _workspace_path(root, relative, f"source:{relative}")
        payload = _read_bytes(path, f"source:{relative}")
        digest.update(relative.encode("ascii"))
        digest.update(b"\0")
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return digest.hexdigest()


def _parse_float(value: str, label: str, *, allow_nan: bool = False) -> float:
    if allow_nan and value == "nan":
        return math.nan
    try:
        result = float(value)
    except ValueError as exc:
        raise ContractError(f"invalid-float:{label}") from exc
    if not math.isfinite(result):
        raise ContractError(f"nonfinite-float:{label}")
    return result


def _is_close(actual: float, expected: float, tolerance: float = 1.0e-14) -> bool:
    return abs(actual - expected) <= tolerance * max(1.0, abs(expected))


def _require_nan(value: float, label: str) -> None:
    if not math.isnan(value):
        raise ContractError(f"expected-nan:{label}")


def _require_finite(value: float, label: str) -> None:
    if not math.isfinite(value):
        raise ContractError(f"expected-finite:{label}")


def _validate_csv(
    payload: bytes,
    material_tolerance: float,
    energy_tolerance: float,
    sensitivity_thresholds: dict[str, Any],
) -> dict[str, float]:
    if b"\r" in payload or not payload.endswith(b"\n"):
        raise ContractError("reference-not-canonical-lf")
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as exc:
        raise ContractError("reference-not-ascii") from exc
    lines = text.splitlines()
    expected_rows = len(STATE_GRID) * len(KIJ_VALUES)
    if len(lines) != expected_rows + 3:
        raise ContractError(f"reference-line-count:{len(lines)}")
    if lines[0] != f"# {CSV_SCHEMA}":
        raise ContractError("reference-schema-marker-mismatch")
    if lines[-1] != f"# RFAM_KIJ_SENSITIVITY_PASS rows={expected_rows}":
        raise ContractError("reference-success-marker-mismatch")

    reader = csv.DictReader(io.StringIO("\n".join(lines[1:-1]) + "\n"))
    if tuple(reader.fieldnames or ()) != FIELDS:
        raise ContractError("reference-field-contract-mismatch")
    raw_rows = list(reader)
    if len(raw_rows) != expected_rows:
        raise ContractError("reference-row-count-mismatch")

    rows_by_case: dict[str, list[dict[str, float | str]]] = {}
    row_index = 0
    for state in STATE_GRID:
        case_id, pressure, temperature, y_o2, _role, expected_region = state
        for expected_kij in KIJ_VALUES:
            raw = raw_rows[row_index]
            label = f"row[{row_index}]"
            row_index += 1
            if None in raw or set(raw) != set(FIELDS):
                raise ContractError(f"reference-row-shape:{label}")
            if raw["case_id"] != case_id:
                raise ContractError(f"case-order-mismatch:{label}")
            kij = _parse_float(raw["kij"], f"{label}.kij")
            if not _is_close(kij, expected_kij):
                raise ContractError(f"kij-grid-mismatch:{label}")
            for field, expected in (
                ("pressure_pa", pressure),
                ("temperature_k", temperature),
                ("y_o2", y_o2),
                ("y_c3h8", 1.0 - y_o2),
            ):
                actual = _parse_float(raw[field], f"{label}.{field}")
                if not _is_close(actual, expected):
                    raise ContractError(f"state-grid-mismatch:{label}.{field}")
            if raw["status"] != "SUCCESS":
                raise ContractError(f"probe-status-not-success:{label}")
            if raw["region"] != expected_region:
                raise ContractError(f"region-mismatch:{label}")
            if raw["quality_defined"] not in {"0", "1"}:
                raise ContractError(f"invalid-quality-flag:{label}")

            values = {
                field: _parse_float(raw[field], f"{label}.{field}", allow_nan=True)
                for field in FIELDS[9:]
            }
            rho = values["rho_kg_m3"]
            energy = values["e_j_kg"]
            _require_finite(rho, f"{label}.rho_kg_m3")
            _require_finite(energy, f"{label}.e_j_kg")
            if rho <= 0.0:
                raise ContractError(f"nonpositive-density:{label}")
            material = values["material_residual"]
            energy_residual = values["energy_identity_residual"]
            if not (0.0 <= material <= material_tolerance):
                raise ContractError(f"material-closure-failed:{label}")
            if not (0.0 <= energy_residual <= energy_tolerance):
                raise ContractError(f"energy-closure-failed:{label}")

            if expected_region == "TWO_PHASE":
                if raw["quality_defined"] != "1":
                    raise ContractError(f"two-phase-quality-undefined:{label}")
                for field in (
                    "vapor_mass_fraction", "vapor_volume_fraction",
                    "liquid_x_o2", "vapor_x_o2", "vapor_y_o2",
                    "vapor_y_c3h8", "vapor_rho_o2_kg_m3",
                    "vapor_rho_c3h8_kg_m3",
                ):
                    _require_finite(values[field], f"{label}.{field}")
                if not 0.0 < values["vapor_mass_fraction"] < 1.0:
                    raise ContractError(f"two-phase-quality-out-of-range:{label}")
                if not 0.0 < values["vapor_volume_fraction"] < 1.0:
                    raise ContractError(f"two-phase-volume-fraction-out-of-range:{label}")
                _require_nan(values["single_rho_o2_kg_m3"], f"{label}.single_rho_o2")
                _require_nan(values["single_rho_c3h8_kg_m3"], f"{label}.single_rho_c3h8")
            elif expected_region == "SINGLE_VAPOR":
                if raw["quality_defined"] != "1":
                    raise ContractError(f"single-vapor-quality-undefined:{label}")
                if values["vapor_mass_fraction"] != 1.0 or values["vapor_volume_fraction"] != 1.0:
                    raise ContractError(f"single-vapor-fraction-mismatch:{label}")
                _require_nan(values["liquid_x_o2"], f"{label}.liquid_x_o2")
                for field in (
                    "vapor_x_o2", "vapor_y_o2", "vapor_y_c3h8",
                    "vapor_rho_o2_kg_m3", "vapor_rho_c3h8_kg_m3",
                ):
                    _require_finite(values[field], f"{label}.{field}")
                _require_nan(values["single_rho_o2_kg_m3"], f"{label}.single_rho_o2")
                _require_nan(values["single_rho_c3h8_kg_m3"], f"{label}.single_rho_c3h8")
            else:
                if raw["quality_defined"] != "0":
                    raise ContractError(f"supercritical-quality-invented:{label}")
                for field in (
                    "vapor_mass_fraction", "vapor_volume_fraction",
                    "liquid_x_o2", "vapor_x_o2", "vapor_y_o2",
                    "vapor_y_c3h8", "vapor_rho_o2_kg_m3",
                    "vapor_rho_c3h8_kg_m3",
                ):
                    _require_nan(values[field], f"{label}.{field}")
                _require_finite(values["single_rho_o2_kg_m3"], f"{label}.single_rho_o2")
                _require_finite(values["single_rho_c3h8_kg_m3"], f"{label}.single_rho_c3h8")

            if expected_region in {"TWO_PHASE", "SINGLE_VAPOR"}:
                if not _is_close(values["vapor_y_o2"] + values["vapor_y_c3h8"], 1.0):
                    raise ContractError(f"vapor-composition-closure-failed:{label}")
            if expected_region == "SINGLE_VAPOR":
                partial_sum = values["vapor_rho_o2_kg_m3"] + values["vapor_rho_c3h8_kg_m3"]
                if not _is_close(partial_sum, rho):
                    raise ContractError(f"single-vapor-density-closure-failed:{label}")
            if expected_region == "SINGLE_SUPERCRITICAL":
                partial_sum = values["single_rho_o2_kg_m3"] + values["single_rho_c3h8_kg_m3"]
                if not _is_close(partial_sum, rho):
                    raise ContractError(f"supercritical-density-closure-failed:{label}")

            parsed: dict[str, float | str] = {"kij": kij, "region": expected_region}
            parsed.update(values)
            rows_by_case.setdefault(case_id, []).append(parsed)

    control = rows_by_case["o2_supercritical_control_6mpa"]
    control_signature = {
        key: value for key, value in control[0].items() if key != "kij"
    }
    for index, row in enumerate(control[1:], start=1):
        signature = {key: value for key, value in row.items() if key != "kij"}
        for key in control_signature:
            left = control_signature[key]
            right = signature[key]
            if isinstance(left, float) and math.isnan(left):
                if not isinstance(right, float) or not math.isnan(right):
                    raise ContractError(f"control-kij-dependence:{index}:{key}")
            elif right != left:
                raise ContractError(f"control-kij-dependence:{index}:{key}")

    metrics: dict[str, float] = {}
    for case_id, threshold_name in (
        ("global_low_2mpa", "global_low_2mpa"),
        ("global_dense_6mpa", "global_dense_6mpa"),
    ):
        rows = rows_by_case[case_id]
        qualities = [float(row["vapor_mass_fraction"]) for row in rows]
        densities = [float(row["rho_kg_m3"]) for row in rows]
        baseline = next(
            float(row["rho_kg_m3"])
            for row in rows
            if _is_close(float(row["kij"]), 0.055)
        )
        quality_span = max(qualities) - min(qualities)
        density_relative_span = (max(densities) - min(densities)) / baseline
        threshold = sensitivity_thresholds[threshold_name]
        if quality_span < _number(threshold["minimum_quality_span"], f"threshold.{case_id}.quality"):
            raise ContractError(f"sensitivity-quality-threshold-failed:{case_id}")
        if density_relative_span < _number(threshold["minimum_density_relative_span"], f"threshold.{case_id}.density"):
            raise ContractError(f"sensitivity-density-threshold-failed:{case_id}")
        metrics[f"{case_id}.quality_span"] = quality_span
        metrics[f"{case_id}.density_relative_span"] = density_relative_span
    return metrics


def validate(root: Path, contract_path: Path, probe: Path | None) -> dict[str, float]:
    document = _object(_decode_json(contract_path), "root", EXPECTED_ROOT_KEYS)
    _exact(document["schema"], SCHEMA, "schema")
    _exact(document["contract_version"], 2, "contract_version")
    _exact(document["status"], "current-source-screening-only", "status")
    _exact(
        document["scope"],
        "offline-o2-c3h8-kij-current-source-priority-screen-no-fit",
        "scope",
    )
    _exact(document["approvals"], EXPECTED_APPROVALS, "approvals")

    upstream = document["upstream"]
    if not isinstance(upstream, list) or len(upstream) != len(EXPECTED_UPSTREAM):
        raise ContractError("upstream-list-contract-mismatch")
    for index, (expected_id, expected) in enumerate(EXPECTED_UPSTREAM.items()):
        item = _object(upstream[index], f"upstream[{index}]", {"id", "path", "sha256"})
        _exact(item["id"], expected_id, f"upstream[{index}].id")
        _exact(item["path"], expected[0], f"upstream[{index}].path")
        _exact(item["sha256"], expected[1], f"upstream[{index}].sha256")
        path = _workspace_path(root, item["path"], f"upstream:{expected_id}")
        if _sha256_bytes(_read_bytes(path, f"upstream:{expected_id}")) != expected[1]:
            raise ContractError(f"upstream-sha256-mismatch:{expected_id}")

    provenance = _object(
        document["core_provenance"],
        "core_provenance",
        set(EXPECTED_CORE_PROVENANCE),
    )
    for provenance_id, expected in EXPECTED_CORE_PROVENANCE.items():
        item = _object(
            provenance[provenance_id],
            f"core_provenance.{provenance_id}",
            {"path", "sha256", "role", "admission"},
        )
        for key, value in expected.items():
            _exact(item[key], value, f"core_provenance.{provenance_id}.{key}")
        path = _workspace_path(root, item["path"], f"core:{provenance_id}")
        if _sha256_bytes(_read_bytes(path, f"core:{provenance_id}")) != expected["sha256"]:
            raise ContractError(f"core-provenance-sha256-mismatch:{provenance_id}")

    _exact(document["pair"], {
        "species": ["O2", "C3H8"],
        "species_order": ["O2", "C3H8"],
        "composition_basis": "mass_fraction",
        "symmetry_required": True,
    }, "pair")

    model = _object(document["model"], "model", {
        "eos", "mixing_rule", "kij_form", "historical_value",
        "historical_value_status", "sweep",
    })
    _exact(model["eos"], "Peng-Robinson", "model.eos")
    _exact(model["mixing_rule"], "van-der-Waals-one-fluid", "model.mixing_rule")
    _exact(model["kij_form"], "constant-symmetric", "model.kij_form")
    _exact(model["historical_value"], 0.055, "model.historical_value")
    _exact(model["historical_value_status"], "legacy-unprovenanced-not-validated", "model.historical_value_status")
    if not isinstance(model["sweep"], list) or len(model["sweep"]) != len(KIJ_VALUES):
        raise ContractError("kij-sweep-contract-mismatch")
    for index, expected in enumerate(KIJ_VALUES):
        if not _is_close(_number(model["sweep"][index], f"model.sweep[{index}]"), expected):
            raise ContractError(f"kij-sweep-contract-mismatch:{index}")

    grid = document["state_grid"]
    if not isinstance(grid, list) or len(grid) != len(STATE_GRID):
        raise ContractError("state-grid-contract-mismatch")
    for index, expected in enumerate(STATE_GRID):
        item = _object(grid[index], f"state_grid[{index}]", {
            "case_id", "pressure_pa", "temperature_k", "y_o2", "y_c3h8",
            "role", "expected_region",
        })
        case_id, pressure, temperature, y_o2, role, region = expected
        _exact(item["case_id"], case_id, f"state_grid[{index}].case_id")
        for key, value in (
            ("pressure_pa", pressure), ("temperature_k", temperature),
            ("y_o2", y_o2), ("y_c3h8", 1.0 - y_o2),
        ):
            if not _is_close(_number(item[key], f"state_grid[{index}].{key}"), value):
                raise ContractError(f"state-grid-contract-mismatch:{index}:{key}")
        _exact(item["role"], role, f"state_grid[{index}].role")
        _exact(item["expected_region"], region, f"state_grid[{index}].expected_region")

    _exact(document["outputs"], list(FIELDS), "outputs")
    tolerances = _object(document["tolerances"], "tolerances", {
        "material_residual_abs", "energy_identity_residual_abs",
        "sensitivity_thresholds",
    })
    material_tolerance = _number(tolerances["material_residual_abs"], "tolerances.material")
    energy_tolerance = _number(tolerances["energy_identity_residual_abs"], "tolerances.energy")
    _exact(material_tolerance, 2.0e-12, "tolerances.material")
    _exact(energy_tolerance, 2.0e-12, "tolerances.energy")
    thresholds = _object(tolerances["sensitivity_thresholds"], "sensitivity_thresholds", {
        "global_low_2mpa", "global_dense_6mpa",
    })
    _exact(thresholds, {
        "global_low_2mpa": {
            "minimum_quality_span": 0.01,
            "minimum_density_relative_span": 0.02,
        },
        "global_dense_6mpa": {
            "minimum_quality_span": 0.05,
            "minimum_density_relative_span": 0.08,
        },
    }, "sensitivity_thresholds")

    bundle = _object(document["source_bundle"], "source_bundle", {
        "algorithm", "files", "sha256",
    })
    _exact(bundle["algorithm"], "sha256(path-nul-size-u64be-bytes)", "source_bundle.algorithm")
    _exact(bundle["files"], list(SOURCE_PATHS), "source_bundle.files")
    bundle_sha = _valid_sha(bundle["sha256"], "source_bundle.sha256")

    reference = _object(document["reference"], "reference", {
        "path", "sha256", "schema", "row_count",
    })
    expected_reference_path = "realFluidAllMach/tests/reference/kij_sensitivity_o2_c3h8_v2.csv"
    _exact(reference["path"], expected_reference_path, "reference.path")
    _exact(reference["schema"], CSV_SCHEMA, "reference.schema")
    _exact(reference["row_count"], 50, "reference.row_count")
    reference_sha = _valid_sha(reference["sha256"], "reference.sha256")
    reference_path = _workspace_path(root, reference["path"], "reference.path")
    reference_payload = _read_bytes(reference_path, "reference")
    if _sha256_bytes(reference_payload) != reference_sha:
        raise ContractError("reference-sha256-mismatch")
    metrics = _validate_csv(
        reference_payload, material_tolerance, energy_tolerance, thresholds
    )
    _exact(reference_sha, EXPECTED_REFERENCE_SHA256, "reference.sha256")

    build = _object(document["build"], "build", {
        "platform", "compiler_path", "compiler_version", "compiler_sha256",
        "flags", "include_directories", "translation_units",
        "expected_executable_sha256",
    })
    _exact(build["platform"], "WSL-x86_64-linux-gnu", "build.platform")
    _exact(build["compiler_path"], "/usr/bin/x86_64-linux-gnu-g++-9", "build.compiler_path")
    _exact(build["compiler_version"], "g++ (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0", "build.compiler_version")
    _exact(build["compiler_sha256"], "0c0b987719385f8e242819dfef6a3e97b4b2540ce4035dbc11565df2966bd69e", "build.compiler_sha256")
    _exact(build["flags"], ["-std=c++17", "-O2", "-fopenmp"], "build.flags")
    _exact(build["include_directories"], list(INCLUDE_DIRECTORIES), "build.include_directories")
    _exact(build["translation_units"], list(TRANSLATION_UNITS), "build.translation_units")
    executable_sha = _valid_sha(build["expected_executable_sha256"], "build.expected_executable_sha256")
    _exact(
        executable_sha,
        EXPECTED_EXECUTABLE_SHA256,
        "build.expected_executable_sha256",
    )

    _exact(document["conclusion"], {
        "pair_priority": "high-data-priority",
        "sensitivity_detected": True,
        "fit_performed": False,
        "external_data_selected": False,
        "validated_kij": False,
        "missing_pair_policy": "fail-closed-until-sourced-and-validated",
    }, "conclusion")
    _exact(document["limitations"], EXPECTED_LIMITATIONS, "limitations")

    if probe is not None:
        probe = probe.resolve()
        if not probe.is_file():
            raise ContractError("probe-executable-missing")
        if _sha256_bytes(_read_bytes(probe, "probe-executable")) != executable_sha:
            raise ContractError("probe-executable-sha256-mismatch")

    # Keep stale live sources fail-closed, but only after the complete static
    # contract, reference, and optional executable identity have been checked.
    # This preserves targeted diagnostics without allowing a stale bundle to
    # reach probe execution or a success marker.
    if source_bundle_sha(root) != bundle_sha:
        raise ContractError("source-bundle-sha256-mismatch")
    _exact(bundle_sha, EXPECTED_SOURCE_BUNDLE_SHA256, "source_bundle.sha256")

    if probe is not None:
        try:
            environment = dict(os.environ)
            environment["OMP_NUM_THREADS"] = "1"
            result = subprocess.run(
                [str(probe)],
                capture_output=True,
                check=False,
                timeout=60,
                env=environment,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise ContractError(f"probe-execution-failed:{exc}") from exc
        if result.returncode != 0:
            raise ContractError(f"probe-return-code:{result.returncode}")
        if result.stderr:
            raise ContractError("probe-stderr-not-empty")
        if result.stdout != reference_payload:
            raise ContractError("probe-reference-byte-mismatch")
    return metrics


def parse_args() -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[2]
    default_contract = (
        Path(__file__).resolve().parents[1]
        / "tests/reference/kij_sensitivity_contract_v2.json"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument("--contract", type=Path, default=default_contract)
    parser.add_argument("--probe-executable", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        metrics = validate(args.root, args.contract, args.probe_executable)
    except Exception as exc:
        reason = str(exc) if isinstance(exc, ContractError) else f"unexpected:{type(exc).__name__}:{exc}"
        print(f"{FAIL_MARKER} reason={reason}")
        return 2

    detail = " ".join(
        f"{key}={value:.9g}" for key, value in sorted(metrics.items())
    )
    if args.probe_executable is None:
        print(f"{BLOCKED_MARKER} rows=50 {detail} validated_kij=0 production_ready=0")
    else:
        print(f"{SUCCESS_MARKER} rows=50 {detail} validated_kij=0 production_ready=0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
