#!/usr/bin/env python3
"""Generate and screen a frozen C3Mech 44I0 0-D reference under WSL.

This is an internal route screen. It does not ingest experimental data, fit
the RFAM two-step mechanism, or enable any runtime chemistry capability.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Iterable, Sequence


class ScreenError(RuntimeError):
    """Fail-closed pipeline error."""


TOOLS_DIR = Path(__file__).resolve().parent
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

import artifact_composite  # noqa: E402
from run_r2_0d_rfam_protocol import parse_runner_output  # noqa: E402


REFERENCE_SCHEMA = "rfam.r2-c3mech-wsl-reference/v1"
SCREEN_SCHEMA = "rfam.r2-c3mech-wsl-screen/v1"
REFERENCE_STATUS = "OUT_OF_MECHANISM_THERMO_ENVELOPE"
SCREEN_STATUS = "INTERNAL_DETAILED_MECHANISM_MISMATCH_SCREEN_ONLY"
FAIL_MARKER = "R2_C3MECH_WSL_SCREEN_FAIL"
REFERENCE_MARKER = "R2_C3MECH_WSL_REFERENCE_BLOCKED"
SCREEN_MARKER = "R2_C3MECH_WSL_SCREEN_BLOCKED"

INPUT_MANIFEST_REL = Path(
    "realFluidAllMach/tests/reference/"
    "r2_c3mech_wsl_screen_inputs_v1_composite.json"
)
INPUT_COMPOSITE_SHA256 = "c04c7699636bf2b2b5dd2b6e900e7eb2b53a6c07d25cdb24438b9b48b5ecb59f"

MECHANISM_REL = Path(
    "tmp/c3mech-v4.0.1-rights-screen/"
    "C3MechV4.0.1_44I0_C0-C3-C4_LT-HT_Cantera.yaml"
)
LICENSE_REL = Path("tmp/c3mech-v4.0.1-rights-screen/LICENSE")
RIGHTS_REPORT_REL = Path(
    "tmp/a408-c3mech-v4.0.1-propane-mechanism-rights-screen-v1-20260728.md"
)
RUNNER_REL = Path("realFluidAllMach/build-phase-foundation/rfam_r2_0d_runner")
MATRIX_REL = Path(
    "realFluidAllMach/tests/reference/r2_c3mech_wsl_cases_v1.json"
)
CASE_IDS = {
    "protocol-training": (
        "train_phi0p8_diluted",
        "train_phi1p0",
        "train_phi1p2",
    ),
    "protocol-holdout": (
        "holdout_phi0p6_diluted",
        "holdout_phi1p1",
        "holdout_phi1p4",
    ),
}
SPECIES = ("C3H8", "O2", "CO", "CO2", "H2O", "N2")
SAMPLE_TIMES = (
    0.0,
    1.0e-14,
    3.0e-14,
    1.0e-13,
    3.0e-13,
    1.0e-12,
    3.0e-12,
    1.0e-11,
    3.0e-11,
    1.0e-10,
    3.0e-10,
    1.0e-9,
    3.0e-9,
    1.0e-8,
    1.5e-8,
    2.0e-8,
    3.0e-8,
    4.0e-8,
    6.0e-8,
    8.0e-8,
    1.0e-7,
    1.5e-7,
    2.0e-7,
    3.0e-7,
    4.0e-7,
    6.0e-7,
    8.0e-7,
    1.0e-6,
    1.5e-6,
    2.0e-6,
    3.0e-6,
    4.0e-6,
    6.0e-6,
    8.0e-6,
    1.0e-5,
    1.5e-5,
    2.0e-5,
    3.0e-5,
    4.0e-5,
    6.0e-5,
    8.0e-5,
    1.0e-4,
    1.5e-4,
    2.0e-4,
    3.0e-4,
    4.0e-4,
    6.0e-4,
    8.0e-4,
    1.0e-3,
    1.5e-3,
    2.0e-3,
    3.0e-3,
    4.0e-3,
    5.0e-3,
)
RUNTIME_IDENTITY = {
    "python_executable": "/usr/bin/python3.8",
    "python_sha256": "298a9e830ed52f36c299427565485d717d1ce0179c0597cc16560513eb780b06",
    "python_size": 5490456,
    "cantera_version": "3.1.0",
    "cantera_module": "/home/nvoyager/.local/lib/python3.8/site-packages/cantera/__init__.py",
    "cantera_module_sha256": "9a7ceae3bae1a7e41dc581e353d44e04af8447114748f5b7c75995e26f8819b9",
    "cantera_module_size": 819,
    "cantera_extension": (
        "/home/nvoyager/.local/lib/python3.8/site-packages/cantera/"
        "_cantera.cpython-38-x86_64-linux-gnu.so"
    ),
    "cantera_extension_sha256": "c57a47c201b0cb8f965b38acc7ebfafbbe95db7464d80359499c64d50555fae8",
    "cantera_extension_size": 13069169,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def strict_pairs(pairs: Iterable[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ScreenError("duplicate-json-key:" + key)
        result[key] = value
    return result


def reject_constant(value: str) -> None:
    raise ScreenError("nonfinite-json-constant:" + value)


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=strict_pairs,
            parse_constant=reject_constant,
        )
    except ScreenError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ScreenError("invalid-json:" + label) from error
    if not isinstance(value, dict):
        raise ScreenError("json-root-not-object:" + label)
    return value


def finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ScreenError("not-a-number:" + label)
    try:
        result = float(value)
    except ValueError as error:
        raise ScreenError("not-a-number:" + label) from error
    if not math.isfinite(result):
        raise ScreenError("nonfinite:" + label)
    return result


def stable(value: float, label: str) -> float:
    if not math.isfinite(value):
        raise ScreenError("nonfinite:" + label)
    return float(format(float(value), ".12g"))


def stable_event(value: float, label: str) -> float:
    if not math.isfinite(value):
        raise ScreenError("nonfinite:" + label)
    return float(format(float(value), ".10g"))


def canonical_bytes(document: dict[str, Any]) -> bytes:
    return (
        json.dumps(document, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n"
    ).encode("ascii")


def publish_no_replace(path: Path, payload: bytes) -> None:
    if os.path.lexists(str(path)):
        raise ScreenError("output-path-already-exists")
    if not path.parent.is_dir():
        raise ScreenError("output-parent-missing")
    descriptor, temporary_name = tempfile.mkstemp(
        prefix="." + path.name + ".", suffix=".tmp", dir=str(path.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as error:
            raise ScreenError("output-publication-race") from error
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def exact_file(root: Path, relative: Path, label: str) -> Path:
    declared = root / relative
    if declared.is_symlink() or not declared.is_file():
        raise ScreenError("missing-or-aliased-authority:" + label)
    resolved = declared.resolve()
    if resolved != declared.absolute():
        raise ScreenError("authority-path-mismatch:" + label)
    return resolved


def verify_authorities(root: Path) -> dict[str, str]:
    manifest = exact_file(root, INPUT_MANIFEST_REL, "input-composite-manifest")
    try:
        actual = artifact_composite.verify_manifest(
            root, manifest, INPUT_COMPOSITE_SHA256
        )
    except artifact_composite.CompositeError as error:
        raise ScreenError("input-composite:" + str(error)) from error
    if actual != INPUT_COMPOSITE_SHA256:
        raise ScreenError("input-composite-return-mismatch")
    return {
        "manifest_path": INPUT_MANIFEST_REL.as_posix(),
        "composite_sha256": INPUT_COMPOSITE_SHA256,
    }


def validate_matrix(root: Path) -> tuple[Path, dict[str, Any]]:
    path = exact_file(root, MATRIX_REL, "case-matrix")
    document = load_json(path, "case-matrix")
    if set(document) != {
        "schema",
        "partition_locked",
        "sample_times_s",
        "species_order",
        "cases",
        "approvals",
    }:
        raise ScreenError("matrix-schema-open")
    if document["schema"] != "rfam.r2-c3mech-wsl-case-matrix/v1":
        raise ScreenError("matrix-schema-mismatch")
    if document["partition_locked"] is not True:
        raise ScreenError("matrix-partition-not-locked")
    if document["species_order"] != list(SPECIES):
        raise ScreenError("matrix-species-order-mismatch")
    times = document["sample_times_s"]
    if not isinstance(times, list) or tuple(finite(item, "sample-time") for item in times) != SAMPLE_TIMES:
        raise ScreenError("matrix-sample-grid-mismatch")
    approvals = document["approvals"]
    if not isinstance(approvals, dict) or not approvals or any(
        value is not False for value in approvals.values()
    ):
        raise ScreenError("matrix-approval-enabled")
    cases = document["cases"]
    expected_ids = CASE_IDS["protocol-training"] + CASE_IDS["protocol-holdout"]
    if not isinstance(cases, list) or tuple(
        case.get("id") for case in cases if isinstance(case, dict)
    ) != expected_ids:
        raise ScreenError("matrix-case-order-mismatch")
    allowed = {
        "id",
        "partition",
        "temperature_k",
        "pressure_pa",
        "equivalence_ratio",
        "fuel",
        "oxidizer",
    }
    for index, case in enumerate(cases):
        if not isinstance(case, dict) or set(case) != allowed:
            raise ScreenError("matrix-case-schema-open")
        expected_partition = (
            "protocol-training" if index < 3 else "protocol-holdout"
        )
        if case["partition"] != expected_partition:
            raise ScreenError("matrix-case-partition-mismatch")
        if finite(case["temperature_k"], "case-temperature") <= 0.0:
            raise ScreenError("matrix-case-temperature-invalid")
        if finite(case["pressure_pa"], "case-pressure") <= 0.0:
            raise ScreenError("matrix-case-pressure-invalid")
        if finite(case["equivalence_ratio"], "case-phi") <= 0.0:
            raise ScreenError("matrix-case-phi-invalid")
        if not isinstance(case["fuel"], str) or not isinstance(case["oxidizer"], str):
            raise ScreenError("matrix-case-composition-invalid")
    return path, document


def load_cantera() -> tuple[Any, dict[str, Any]]:
    try:
        ct = importlib.import_module("cantera")
    except Exception as error:
        raise ScreenError("cantera-import-failed") from error
    executable = Path(sys.executable).resolve()
    module = Path(ct.__file__).resolve()
    extension = Path(ct._cantera.__file__).resolve()
    actual = {
        "python_executable": str(executable),
        "python_sha256": sha256(executable),
        "python_size": executable.stat().st_size,
        "cantera_version": ct.__version__,
        "cantera_module": str(module),
        "cantera_module_sha256": sha256(module),
        "cantera_module_size": module.stat().st_size,
        "cantera_extension": str(extension),
        "cantera_extension_sha256": sha256(extension),
        "cantera_extension_size": extension.stat().st_size,
    }
    if actual != RUNTIME_IDENTITY:
        raise ScreenError("cantera-runtime-identity-mismatch")
    if hasattr(ct, "suppress_thermo_warnings"):
        ct.suppress_thermo_warnings()
    return ct, actual


def load_mechanism(ct: Any, path: Path) -> tuple[Any, dict[str, Any]]:
    try:
        gas = ct.Solution(str(path))
    except Exception as error:
        raise ScreenError("mechanism-load-failed") from error
    missing = [name for name in SPECIES if name not in gas.species_names]
    inventory = {
        "model_id": "C3Mech-V4.0.1-44I0",
        "species_count": gas.n_species,
        "reaction_count": gas.n_reactions,
        "thermo_model": gas.thermo_model,
        "thermo_min_k": stable(float(gas.min_temp), "thermo-min"),
        "thermo_max_k": stable(float(gas.max_temp), "thermo-max"),
    }
    expected = {
        "model_id": "C3Mech-V4.0.1-44I0",
        "species_count": 729,
        "reaction_count": 4080,
        "thermo_model": "ideal-gas",
        "thermo_min_k": 400.0,
        "thermo_max_k": 2000.0,
    }
    if inventory != expected or missing:
        raise ScreenError("mechanism-inventory-mismatch")
    return gas, inventory


def mass_fractions(phase: Any) -> dict[str, float]:
    values = {name: stable(float(phase[name].Y[0]), "Y." + name) for name in SPECIES}
    if any(value < 0.0 or value > 1.0 for value in values.values()):
        raise ScreenError("species-mass-fraction-out-of-range")
    return values


def state_record(time_s: float, reactor: Any) -> dict[str, Any]:
    phase = reactor.thermo
    return {
        "time_s": stable(time_s, "sample-time"),
        "status": "eligible-within-mechanism-thermo-intersection",
        "temperature_k": stable(float(phase.T), "sample-temperature"),
        "pressure_pa": stable(float(phase.P), "sample-pressure"),
        "density_kg_m3": stable(float(phase.density), "sample-density"),
        "heat_release_rate_w_m3": stable(
            float(phase.heat_release_rate), "sample-qdot"
        ),
        "mass_fractions": mass_fractions(phase),
    }


def crossing_time(
    previous_time: float,
    previous_temperature: float,
    current_time: float,
    current_temperature: float,
    threshold: float,
    label: str,
) -> float:
    span = current_temperature - previous_temperature
    if span <= 0.0:
        raise ScreenError("temperature-crossing-not-bracketed:" + label)
    fraction = (threshold - previous_temperature) / span
    if fraction < 0.0 or fraction > 1.0:
        raise ScreenError("temperature-crossing-fraction-invalid:" + label)
    return stable_event(previous_time + fraction * (current_time - previous_time), label)


def run_case(ct: Any, mechanism: Path, definition: dict[str, Any]) -> dict[str, Any]:
    gas, _ = load_mechanism(ct, mechanism)
    gas.TP = float(definition["temperature_k"]), float(definition["pressure_pa"])
    gas.set_equivalence_ratio(
        float(definition["equivalence_ratio"]),
        definition["fuel"],
        definition["oxidizer"],
    )
    try:
        reactor = ct.IdealGasConstPressureReactor(gas, energy="on")
        network = ct.ReactorNet([reactor])
        network.rtol = 1.0e-10
        network.atol = 1.0e-18
        reactor.set_advance_limit("temperature", 1.0)
    except Exception as error:
        raise ScreenError("reactor-initialization-failed") from error

    initial_temperature = float(reactor.T)
    ignition_threshold = initial_temperature + 400.0
    thermo_max = 2000.0
    ignition_delay: float | None = None
    ignition_status = (
        "threshold-coincides-with-thermo-bound-unscored"
        if ignition_threshold >= thermo_max
        else "not-reached-within-thermo-envelope"
    )
    exit_record: dict[str, Any] | None = None
    samples: list[dict[str, Any]] = [state_record(0.0, reactor)]
    if abs(sum(samples[0]["mass_fractions"].values()) - 1.0) > 2.0e-9:
        raise ScreenError("initial-six-species-mass-fraction-sum-mismatch")

    for target in SAMPLE_TIMES[1:]:
        if exit_record is not None:
            samples.append(
                {
                    "time_s": stable(target, "sample-time"),
                    "status": "not-evaluated-after-thermo-envelope-exit",
                }
            )
            continue
        while network.time < target:
            previous_time = float(network.time)
            previous_temperature = float(reactor.T)
            try:
                reached = float(network.advance(target, apply_limit=True))
            except Exception as error:
                raise ScreenError("cantera-advance-failed:" + definition["id"]) from error
            current_temperature = float(reactor.T)
            if (
                not math.isfinite(reached)
                or reached <= previous_time
                or reached > target
                or not math.isfinite(current_temperature)
            ):
                raise ScreenError("cantera-time-or-temperature-progress-invalid")

            if (
                ignition_delay is None
                and ignition_threshold < thermo_max
                and previous_temperature < ignition_threshold <= current_temperature
            ):
                ignition_delay = crossing_time(
                    previous_time,
                    previous_temperature,
                    reached,
                    current_temperature,
                    ignition_threshold,
                    "ignition-delay",
                )
                ignition_status = "scored-within-mechanism-thermo-intersection"

            if previous_temperature < thermo_max <= current_temperature:
                exit_record = {
                    "status": "thermo-envelope-exit-detected-no-endpoint-state-published",
                    "time_s": crossing_time(
                        previous_time,
                        previous_temperature,
                        reached,
                        current_temperature,
                        thermo_max,
                        "thermo-exit",
                    ),
                    "temperature_k": thermo_max,
                    "publishable_chemistry_state": False,
                    "interpolation": "1-K-advance-limited-linear-temperature",
                }
                break

        if exit_record is None:
            samples.append(state_record(target, reactor))
        else:
            samples.append(
                {
                    "time_s": stable(target, "sample-time"),
                    "status": "outside-thermo-envelope-unscored",
                }
            )

    if exit_record is None:
        raise ScreenError("thermo-envelope-exit-not-observed:" + definition["id"])
    if ignition_status == "not-reached-within-thermo-envelope":
        raise ScreenError("ignition-not-reached-before-thermo-exit:" + definition["id"])

    return {
        **definition,
        "initial_mass_fractions": samples[0]["mass_fractions"],
        "ignition_delay": {
            "definition": "first-1K-limited-linear-crossing-of-T-minus-T0-equals-400K",
            "status": ignition_status,
            "time_s": ignition_delay,
        },
        "thermo_envelope_exit": exit_record,
        "samples": samples,
    }


def generate_reference(root: Path) -> dict[str, Any]:
    authority = verify_authorities(root)
    ct, runtime = load_cantera()
    mechanism = exact_file(root, MECHANISM_REL, "c3mech-44i0")
    exact_file(root, LICENSE_REL, "c3mech-license")
    exact_file(root, RIGHTS_REPORT_REL, "c3mech-rights-report")
    _, inventory = load_mechanism(ct, mechanism)
    matrix_path, matrix = validate_matrix(root)
    cases = [run_case(ct, mechanism, definition) for definition in matrix["cases"]]
    return {
        "schema": REFERENCE_SCHEMA,
        "status": REFERENCE_STATUS,
        "generator_sha256": sha256(Path(__file__).resolve()),
        "authority": authority,
        "runtime_identity": runtime,
        "mechanism": inventory,
        "rights_boundary": {
            "mechanism_code_license": "CC-BY-4.0-candidate",
            "experimental_training_rights_cleared": False,
            "experimental_holdout_rights_cleared": False,
            "experimental_values_ingested": False,
            "physical_fit_or_validation_claim": False,
        },
        "reactor": {
            "model": "IdealGasConstPressureReactor",
            "energy": "on",
            "rtol": 1.0e-10,
            "atol": 1.0e-18,
            "temperature_advance_limit_k": 1.0,
            "thermo_exit_policy": "stop-before-publishing-2000K-endpoint",
        },
        "matrices": [
            {
                "path": matrix_path.relative_to(root).as_posix(),
                "training_case_ids": list(CASE_IDS["protocol-training"]),
                "holdout_case_ids": list(CASE_IDS["protocol-holdout"]),
                "partition_locked_before_screen": True,
            }
        ],
        "sample_times_s": list(SAMPLE_TIMES),
        "species_order": list(SPECIES),
        "observables": [
            "ignition_delay_s",
            "temperature_k",
            "heat_release_rate_w_m3",
            *["Y_" + name for name in SPECIES],
        ],
        "cases": cases,
        "approvals": {
            "r1_r2_fit_approved": False,
            "physical_validation_approved": False,
            "phase_aware_chemistry_enabled": False,
            "fluent_execution_approved": False,
            "production_ready": False,
        },
    }


def require_exact_keys(value: Any, keys: set[str], label: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != keys:
        raise ScreenError("schema-open:" + label)
    return value


def validate_reference(document: dict[str, Any]) -> None:
    root = require_exact_keys(
        document,
        {
            "schema",
            "status",
            "generator_sha256",
            "authority",
            "runtime_identity",
            "mechanism",
            "rights_boundary",
            "reactor",
            "matrices",
            "sample_times_s",
            "species_order",
            "observables",
            "cases",
            "approvals",
        },
        "reference",
    )
    if root["schema"] != REFERENCE_SCHEMA or root["status"] != REFERENCE_STATUS:
        raise ScreenError("reference-schema-or-status-mismatch")
    if root["generator_sha256"] != sha256(Path(__file__).resolve()):
        raise ScreenError("reference-generator-identity-mismatch")
    if root["authority"] != {
        "manifest_path": INPUT_MANIFEST_REL.as_posix(),
        "composite_sha256": INPUT_COMPOSITE_SHA256,
    }:
        raise ScreenError("reference-authority-mismatch")
    if root["runtime_identity"] != RUNTIME_IDENTITY:
        raise ScreenError("reference-runtime-identity-mismatch")
    if root["species_order"] != list(SPECIES):
        raise ScreenError("reference-species-order-mismatch")
    if root["sample_times_s"] != list(SAMPLE_TIMES):
        raise ScreenError("reference-sample-grid-mismatch")
    approvals = root["approvals"]
    if not isinstance(approvals, dict) or not approvals or any(
        value is not False for value in approvals.values()
    ):
        raise ScreenError("reference-approval-enabled")
    rights = root["rights_boundary"]
    if not isinstance(rights, dict) or any(
        rights.get(key) is not False
        for key in (
            "experimental_training_rights_cleared",
            "experimental_holdout_rights_cleared",
            "experimental_values_ingested",
            "physical_fit_or_validation_claim",
        )
    ):
        raise ScreenError("reference-rights-boundary-weakened")
    cases = root["cases"]
    expected_ids = CASE_IDS["protocol-training"] + CASE_IDS["protocol-holdout"]
    if not isinstance(cases, list) or tuple(
        case.get("id") for case in cases if isinstance(case, dict)
    ) != expected_ids:
        raise ScreenError("reference-case-order-mismatch")
    for case in cases:
        samples = case.get("samples")
        if not isinstance(samples, list) or len(samples) != len(SAMPLE_TIMES):
            raise ScreenError("reference-case-sample-count-mismatch")
        initial = case.get("initial_mass_fractions")
        if not isinstance(initial, dict) or set(initial) != set(SPECIES):
            raise ScreenError("reference-initial-composition-shape")
        if abs(sum(finite(initial[name], "initial-Y") for name in SPECIES) - 1.0) > 2.0e-9:
            raise ScreenError("reference-initial-composition-sum")
        saw_exit = False
        for index, (sample, expected_time) in enumerate(zip(samples, SAMPLE_TIMES)):
            if not isinstance(sample, dict):
                raise ScreenError("reference-sample-not-object")
            if finite(sample.get("time_s"), "sample-time") != expected_time:
                raise ScreenError("reference-sample-time-mismatch")
            status = sample.get("status")
            if status == "eligible-within-mechanism-thermo-intersection":
                if saw_exit or set(sample) != {
                    "time_s",
                    "status",
                    "temperature_k",
                    "pressure_pa",
                    "density_kg_m3",
                    "heat_release_rate_w_m3",
                    "mass_fractions",
                }:
                    raise ScreenError("reference-eligible-sample-shape")
                composition = sample["mass_fractions"]
                if not isinstance(composition, dict) or set(composition) != set(SPECIES):
                    raise ScreenError("reference-sample-composition-shape")
                finite(sample["temperature_k"], "sample-temperature")
                finite(sample["heat_release_rate_w_m3"], "sample-qdot")
            elif status == "outside-thermo-envelope-unscored":
                if set(sample) != {"time_s", "status"} or saw_exit:
                    raise ScreenError("reference-exit-sample-shape")
                saw_exit = True
            elif status == "not-evaluated-after-thermo-envelope-exit":
                if set(sample) != {"time_s", "status"} or not saw_exit:
                    raise ScreenError("reference-post-exit-sample-shape")
            else:
                raise ScreenError("reference-sample-status-invalid:" + str(index))
        if not saw_exit:
            raise ScreenError("reference-case-missing-thermo-exit")


def run_rfam_case(runner: Path, case: dict[str, Any]) -> tuple[dict[str, str], list[dict[str, Any]]]:
    initial = case["initial_mass_fractions"]
    command = [
        str(runner),
        case["id"],
        format(float(case["temperature_k"]), ".17g"),
        format(float(case["pressure_pa"]), ".17g"),
        str(len(SAMPLE_TIMES)),
        "0",
        *[format(float(initial[name]), ".17g") for name in SPECIES],
        *[format(value, ".17g") for value in SAMPLE_TIMES],
    ]
    result = subprocess.run(
        command, check=False, text=True, capture_output=True, timeout=60
    )
    if result.returncode != 0:
        raise ScreenError(
            "rfam-runner-rejected:"
            + case["id"]
            + ":rc="
            + str(result.returncode)
        )
    try:
        return parse_runner_output(
            result.stdout, result.stderr, case["id"], len(SAMPLE_TIMES)
        )
    except Exception as error:
        raise ScreenError("rfam-runner-output-invalid:" + case["id"]) from error


def rms(values: Sequence[float]) -> float:
    if not values:
        raise ScreenError("empty-rms-input")
    return math.sqrt(math.fsum(value * value for value in values) / len(values))


def compare_case(
    case: dict[str, Any],
    headers: dict[str, str],
    rfam_samples: list[dict[str, Any]],
) -> dict[str, Any]:
    paired: list[tuple[dict[str, Any], dict[str, Any]]] = []
    for detailed, rfam in zip(case["samples"], rfam_samples):
        if detailed["status"] != "eligible-within-mechanism-thermo-intersection":
            continue
        temperature = finite(rfam["temperature_k"], "rfam-temperature")
        if 300.0 <= temperature <= 3400.0:
            paired.append((detailed, rfam))
    if len(paired) < 2:
        raise ScreenError("insufficient-joint-prefix:" + case["id"])

    temperature_errors = [
        finite(rfam["temperature_k"], "rfam-temperature")
        - finite(detailed["temperature_k"], "detailed-temperature")
        for detailed, rfam in paired
    ]
    qdot_reference = [
        finite(detailed["heat_release_rate_w_m3"], "detailed-qdot")
        for detailed, _ in paired
    ]
    qdot_errors = [
        finite(rfam["qdot_w_m3"], "rfam-qdot") - reference
        for reference, (_, rfam) in zip(qdot_reference, paired)
    ]
    qdot_scale = max(max(abs(value) for value in qdot_reference), 1.0e8)
    species_metrics: dict[str, Any] = {}
    for species_index, name in enumerate(SPECIES):
        errors = [
            finite(rfam["mass_fractions"][species_index], "rfam-Y")
            - finite(detailed["mass_fractions"][name], "detailed-Y")
            for detailed, rfam in paired
        ]
        species_metrics[name] = {
            "rmse_mass_fraction": stable(rms(errors), "species-rmse"),
            "max_abs_mass_fraction": stable(
                max(abs(value) for value in errors), "species-max"
            ),
        }

    detailed_idt = case["ignition_delay"]
    rfam_detected = headers.get("ignition_detected") == "1"
    rfam_idt = finite(headers["ignition_delay_s"], "rfam-idt") if rfam_detected else None
    idt_ratio = None
    if detailed_idt["status"] == "scored-within-mechanism-thermo-intersection":
        if not rfam_detected or rfam_idt is None or rfam_idt <= 0.0:
            raise ScreenError("rfam-ignition-not-detected:" + case["id"])
        idt_ratio = rfam_idt / finite(detailed_idt["time_s"], "detailed-idt")

    detailed_qdot0 = finite(paired[0][0]["heat_release_rate_w_m3"], "qdot0")
    rfam_qdot0 = finite(paired[0][1]["qdot_w_m3"], "rfam-qdot0")
    return {
        "case_id": case["id"],
        "partition": case["partition"],
        "joint_prefix_sample_count": len(paired),
        "joint_prefix_last_time_s": stable(
            finite(paired[-1][0]["time_s"], "last-time"), "last-time"
        ),
        "detailed_ignition_delay": detailed_idt,
        "rfam_ignition_detected": rfam_detected,
        "rfam_ignition_delay_s": stable(rfam_idt, "rfam-idt") if rfam_idt else None,
        "rfam_to_detailed_idt_ratio": stable(idt_ratio, "idt-ratio") if idt_ratio else None,
        "initial_qdot_sign_matches": (
            (detailed_qdot0 == 0.0 and rfam_qdot0 == 0.0)
            or (detailed_qdot0 * rfam_qdot0 > 0.0)
        ),
        "initial_qdot_absolute_ratio": stable(
            abs(rfam_qdot0) / max(abs(detailed_qdot0), 1.0), "qdot0-ratio"
        ),
        "temperature_rmse_k": stable(rms(temperature_errors), "temperature-rmse"),
        "temperature_max_abs_error_k": stable(
            max(abs(value) for value in temperature_errors), "temperature-max"
        ),
        "qdot_reference_scale_w_m3": stable(qdot_scale, "qdot-scale"),
        "qdot_normalized_rmse": stable(
            rms(qdot_errors) / qdot_scale, "qdot-normalized-rmse"
        ),
        "species": species_metrics,
    }


def partition_summary(cases: list[dict[str, Any]], partition: str) -> dict[str, Any]:
    selected = [case for case in cases if case["partition"] == partition]
    ratios = [
        case["rfam_to_detailed_idt_ratio"]
        for case in selected
        if case["rfam_to_detailed_idt_ratio"] is not None
    ]
    return {
        "case_count": len(selected),
        "joint_prefix_sample_count": sum(
            case["joint_prefix_sample_count"] for case in selected
        ),
        "scored_idt_count": len(ratios),
        "rfam_to_detailed_idt_ratio_min": stable(min(ratios), "idt-min"),
        "rfam_to_detailed_idt_ratio_max": stable(max(ratios), "idt-max"),
        "maximum_temperature_error_k": stable(
            max(case["temperature_max_abs_error_k"] for case in selected),
            "partition-temperature-max",
        ),
        "maximum_qdot_normalized_rmse": stable(
            max(case["qdot_normalized_rmse"] for case in selected),
            "partition-qdot-max",
        ),
        "initial_qdot_sign_mismatch_count": sum(
            not case["initial_qdot_sign_matches"] for case in selected
        ),
    }


def generate_screen(root: Path, reference: dict[str, Any]) -> dict[str, Any]:
    verify_authorities(root)
    validate_reference(reference)
    runner = exact_file(root, RUNNER_REL, "rfam-runner")
    comparisons = []
    for case in reference["cases"]:
        headers, samples = run_rfam_case(runner, case)
        comparisons.append(compare_case(case, headers, samples))
    return {
        "schema": SCREEN_SCHEMA,
        "status": SCREEN_STATUS,
        "authority": reference["authority"],
        "reference_sha256": hashlib.sha256(canonical_bytes(reference)).hexdigest(),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "comparison_basis": {
            "reactor": "ideal-gas-constant-pressure",
            "sample_axis": "same-predeclared-fixed-grid",
            "joint_prefix": "detailed-eligible-and-rfam-300-to-3400K",
            "heat_release_basis_fit_approved": False,
            "mechanism_or_experiment_truth_claim": False,
        },
        "partitions": {
            partition: partition_summary(comparisons, partition)
            for partition in ("protocol-training", "protocol-holdout")
        },
        "cases": comparisons,
        "approvals": {
            "r1_r2_fit_approved": False,
            "physical_validation_approved": False,
            "phase_aware_chemistry_enabled": False,
            "fluent_execution_approved": False,
            "production_ready": False,
        },
        "limitations": [
            "C3Mech mechanism code is a candidate under CC BY 4.0; experimental data rights remain unresolved",
            "this comparison is detailed-mechanism screening, not experimental validation",
            "states at or above the 2000 K C3Mech thermo intersection are unscored",
            "no parameter fit, uncertainty quantification, HEM chemistry, Fluent run, or production approval is implied",
        ],
    }


def run_pipeline(root: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    reference = generate_reference(root)
    screen = generate_screen(root, reference)
    return reference, screen


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workspace-root", type=Path, default=Path(__file__).resolve().parents[2]
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    reference_parser = subparsers.add_parser("reference")
    reference_parser.add_argument("--output", type=Path, required=True)
    screen_parser = subparsers.add_parser("screen")
    screen_parser.add_argument("--reference", type=Path, required=True)
    screen_parser.add_argument("--output", type=Path, required=True)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--reference-output", type=Path, required=True)
    run_parser.add_argument("--screen-output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        root = args.workspace_root.resolve()
        if args.command == "reference":
            document = generate_reference(root)
            publish_no_replace(args.output.resolve(), canonical_bytes(document))
            print(
                REFERENCE_MARKER
                + " cases="
                + str(len(document["cases"]))
                + " samples="
                + str(len(SAMPLE_TIMES))
            )
        elif args.command == "screen":
            reference = load_json(args.reference.resolve(), "reference")
            document = generate_screen(root, reference)
            publish_no_replace(args.output.resolve(), canonical_bytes(document))
            print(
                SCREEN_MARKER
                + " cases="
                + str(len(document["cases"]))
                + " fit_approved=0 production_ready=0"
            )
        else:
            reference, screen = run_pipeline(root)
            reference_output = args.reference_output.resolve()
            screen_output = args.screen_output.resolve()
            if os.path.lexists(str(reference_output)) or os.path.lexists(str(screen_output)):
                raise ScreenError("output-path-already-exists")
            publish_no_replace(reference_output, canonical_bytes(reference))
            publish_no_replace(screen_output, canonical_bytes(screen))
            print(
                SCREEN_MARKER
                + " cases="
                + str(len(screen["cases"]))
                + " reference_cases="
                + str(len(reference["cases"]))
                + " fit_approved=0 production_ready=0"
            )
        return 0
    except Exception as error:
        print(
            FAIL_MARKER + " " + type(error).__name__ + ":" + str(error),
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
