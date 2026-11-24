#!/usr/bin/env python3
"""Fail-closed validator for the blocked local Cantera R2 0D protocol.

This gate establishes only reproducible local-reference provenance and a
predeclared protocol split.  It deliberately cannot approve a fit, select a
detailed mechanism, select experimental data, or claim physical validation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path, PurePosixPath
from typing import Any, Callable


class ContractError(RuntimeError):
    pass


EXPECTED_UPSTREAM = [
    {
        "id": "wp2_same_mechanism_instantaneous_reference",
        "path": "tests/reference/r2_instantaneous_reference_v1.json",
        "sha256": "8bd3f2c05d70bb59af7668aa37c35bab3726c79d5ea8f3b650ba065980c1840a",
        "scope": "implementation-equivalence-only",
    },
    {
        "id": "wp2_detailed_mechanism_candidate",
        "path": "tests/reference/r2_physics_candidate_v1.json",
        "sha256": "cdd6046bab363ca5c9b0c637d341b9ff4ca3501b6e0f9dd142fdca1d17ba59d5",
        "scope": "blocked-provenance-screening-only",
    },
]
EXPECTED_ARTIFACTS = [
    {
        "id": "local_cantera_mechanism",
        "path": "tests/reference/cantera_gri30_v3_0.yaml",
        "sha256": "06650b1e0ee0012f6903d5328b1bb218cb6007d07f8ebe375d18f24811039345",
        "kind": "checked-in-local-cantera-ideal-gas-mechanism",
    },
    {
        "id": "reference_generator",
        "path": "scripts/generate_r2_0d_calibration_reference.py",
        "sha256": "43ad404b9eb63a1fb6dd343023a6f726cb5efd4f3e84229c9a50336fb460ed49",
        "kind": "local-cantera-0d-reference-generator",
    },
    {
        "id": "reference_output",
        "path": "tests/reference/r2_0d_cantera_reference_v1.json",
        "sha256": "11d1135a0d0beae0138a3b1dbdde9cae6ccdcee092c103d537311e72f497981d",
        "kind": "generated-local-cantera-0d-reference",
    },
]
EXPECTED_APPROVALS = {
    "local_reference_selected_for_fit": False,
    "detailed_mechanism_selected_for_fit": False,
    "experimental_data_selected": False,
    "r2_fitting_approved": False,
    "physical_validation_approved": False,
    "production_ready": False,
    "reacting_hot_start_approved": False,
}
EXPECTED_CASES = [
    ("train_phi0p8_diluted", "protocol-training"),
    ("train_phi1p0", "protocol-training"),
    ("train_phi1p2", "protocol-training"),
    ("holdout_phi0p6_diluted", "protocol-holdout"),
    ("holdout_phi1p1", "protocol-holdout"),
    ("holdout_phi1p4", "protocol-holdout"),
]
EXPECTED_SAMPLE_TIMES = [
    0.0,
    1.0e-8,
    3.0e-8,
    1.0e-7,
    3.0e-7,
    1.0e-6,
    3.0e-6,
    1.0e-5,
    3.0e-5,
    1.0e-4,
    3.0e-4,
    1.0e-3,
    3.0e-3,
    5.0e-3,
]
EXPECTED_MAJOR_SPECIES = ["C3H8", "O2", "CO", "CO2", "H2O", "N2"]
EXPECTED_OBSERVABLES = [
    "ignition_delay",
    "temperature_history",
    "C3H8_mass_fraction_history",
    "O2_mass_fraction_history",
    "CO_mass_fraction_history",
    "CO2_mass_fraction_history",
    "H2O_mass_fraction_history",
    "Cantera_heat_release_rate_history",
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _strict_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ContractError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def _load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_strict_object)
    except ContractError:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContractError(f"json-load-failed:{path.name}:{error}") from error


def _object(value: Any, location: str, keys: set[str]) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ContractError(f"invalid-object:{location}")
    actual = set(value)
    missing = sorted(keys - actual)
    unknown = sorted(actual - keys)
    if missing:
        raise ContractError(f"missing-keys:{location}:{','.join(missing)}")
    if unknown:
        raise ContractError(f"unknown-keys:{location}:{','.join(unknown)}")
    return value


def _array(value: Any, location: str) -> list[Any]:
    if not isinstance(value, list):
        raise ContractError(f"invalid-array:{location}")
    return value


def _exact(value: Any, expected: Any, location: str) -> None:
    if value != expected:
        raise ContractError(f"contract-value-mismatch:{location}")


def _bool(value: Any, location: str) -> bool:
    if not isinstance(value, bool):
        raise ContractError(f"invalid-bool:{location}")
    return value


def _finite(value: Any, location: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"invalid-number:{location}")
    number = float(value)
    if not math.isfinite(number):
        raise ContractError(f"nonfinite-number:{location}")
    return number


def _relative_path(root: Path, value: Any, location: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ContractError(f"invalid-path:{location}")
    normalized = PurePosixPath(value)
    if normalized.is_absolute() or ".." in normalized.parts or "." in normalized.parts:
        raise ContractError(f"path-escapes-root:{location}")
    candidate = (root / Path(*normalized.parts)).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as error:
        raise ContractError(f"path-escapes-root:{location}") from error
    return candidate


def _verify_file_binding(root: Path, item: dict[str, Any], expected: dict[str, str], location: str) -> Path:
    _object(item, location, {"id", "path", "sha256", "kind"})
    for key, expected_value in expected.items():
        if item[key] != expected_value:
            if key == "path":
                raise ContractError(f"artifact-path-contract-mismatch:{expected['id']}")
            raise ContractError(f"artifact-{key}-contract-mismatch:{expected['id']}")
    path = _relative_path(root, item["path"], f"{location}.path")
    if not path.is_file():
        raise ContractError(f"artifact-missing:{expected['id']}")
    actual_hash = _sha256(path)
    if actual_hash != expected["sha256"]:
        raise ContractError(f"artifact-sha256-mismatch:{expected['id']}")
    return path


def _verify_upstream(root: Path, items: Any) -> None:
    records = _array(items, "upstream")
    if len(records) != len(EXPECTED_UPSTREAM):
        raise ContractError("upstream-count-mismatch")
    for index, expected in enumerate(EXPECTED_UPSTREAM):
        item = _object(records[index], f"upstream[{index}]", {"id", "path", "sha256", "scope"})
        if item != expected:
            if item.get("path") != expected["path"]:
                raise ContractError("upstream-path-contract-mismatch:" + expected["id"])
            raise ContractError("upstream-contract-mismatch:" + expected["id"])
        path = _relative_path(root, item["path"], f"upstream[{index}].path")
        if not path.is_file():
            raise ContractError("upstream-missing:" + expected["id"])
        if _sha256(path) != expected["sha256"]:
            raise ContractError("upstream-sha256-mismatch:" + expected["id"])


def _verify_reference(root: Path, artifact_paths: dict[str, Path]) -> None:
    reference = _load_json(artifact_paths["reference_output"])
    _object(
        reference,
        "reference",
        {
            "schema",
            "status",
            "cantera_version",
            "runtime_identity",
            "generator_sha256",
            "mechanism_sha256",
            "reactor",
            "sample_times_s",
            "major_species",
            "units",
            "cases",
        },
    )
    _exact(reference["schema"], "rfam.r2-0d-cantera-reference/v1", "reference.schema")
    _exact(reference["status"], "local-cantera-reference-not-fit-eligible", "reference.status")
    _exact(reference["cantera_version"], "3.1.0", "reference.cantera_version")
    _exact(
        reference["runtime_identity"],
        {
            "python_executable": "/usr/bin/python3.8",
            "python_sha256": "298a9e830ed52f36c299427565485d717d1ce0179c0597cc16560513eb780b06",
            "cantera_module": "/home/nvoyager/.local/lib/python3.8/site-packages/cantera/__init__.py",
            "cantera_module_sha256": "9a7ceae3bae1a7e41dc581e353d44e04af8447114748f5b7c75995e26f8819b9",
            "cantera_extension": "/home/nvoyager/.local/lib/python3.8/site-packages/cantera/_cantera.cpython-38-x86_64-linux-gnu.so",
            "cantera_extension_sha256": "c57a47c201b0cb8f965b38acc7ebfafbbe95db7464d80359499c64d50555fae8",
        },
        "reference.runtime_identity",
    )
    _exact(reference["generator_sha256"], EXPECTED_ARTIFACTS[1]["sha256"], "reference.generator_sha256")
    _exact(reference["mechanism_sha256"], EXPECTED_ARTIFACTS[0]["sha256"], "reference.mechanism_sha256")
    reactor = _object(reference["reactor"], "reference.reactor", {"model", "phase", "energy", "rtol", "atol", "phase_scope"})
    _exact(reactor["model"], "IdealGasConstPressureReactor", "reference.reactor.model")
    _exact(reactor["phase"], "ideal-gas", "reference.reactor.phase")
    _exact(reactor["energy"], "on", "reference.reactor.energy")
    _exact(_finite(reactor["rtol"], "reference.reactor.rtol"), 1.0e-10, "reference.reactor.rtol")
    _exact(_finite(reactor["atol"], "reference.reactor.atol"), 1.0e-18, "reference.reactor.atol")
    _exact(
        reactor["phase_scope"],
        "high-temperature gas-only; VLE, near-critical and solid-risk states excluded",
        "reference.reactor.phase_scope",
    )
    times = [_finite(value, "reference.sample_times_s") for value in _array(reference["sample_times_s"], "reference.sample_times_s")]
    _exact(times, EXPECTED_SAMPLE_TIMES, "reference.sample_times_s")
    _exact(reference["major_species"], EXPECTED_MAJOR_SPECIES, "reference.major_species")
    _exact(reference["units"], {"temperature": "K", "pressure": "Pa", "density": "kg/m3", "heat_release_rate": "W/m3", "mass_fraction": "1", "time": "s"}, "reference.units")

    cases = _array(reference["cases"], "reference.cases")
    if len(cases) != len(EXPECTED_CASES):
        raise ContractError("reference-case-count-mismatch")
    for index, (case_id, partition) in enumerate(EXPECTED_CASES):
        case = _object(
            cases[index],
            f"reference.cases[{index}]",
            {
                "id", "partition", "temperature_k", "pressure_pa", "equivalence_ratio", "fuel", "oxidizer",
                "ignition_delay_definition", "ignition_delay_s", "initial_major_species_mass_fractions", "samples",
            },
        )
        _exact(case["id"], case_id, f"reference.cases[{index}].id")
        _exact(case["partition"], partition, f"reference.cases[{index}].partition")
        if _finite(case["temperature_k"], f"reference.cases[{index}].temperature_k") < 1000.0:
            raise ContractError("reference-case-temperature-out-of-scope:" + case_id)
        if _finite(case["pressure_pa"], f"reference.cases[{index}].pressure_pa") <= 0.0:
            raise ContractError("reference-case-pressure-invalid:" + case_id)
        if _finite(case["equivalence_ratio"], f"reference.cases[{index}].equivalence_ratio") <= 0.0:
            raise ContractError("reference-case-phi-invalid:" + case_id)
        _exact(case["ignition_delay_definition"], "first adaptive integrator time where T - T0 >= 400 K", f"reference.cases[{index}].ignition_delay_definition")
        if _finite(case["ignition_delay_s"], f"reference.cases[{index}].ignition_delay_s") <= 0.0:
            raise ContractError("reference-case-no-ignition:" + case_id)
        initial = _object(case["initial_major_species_mass_fractions"], f"reference.cases[{index}].initial_major_species_mass_fractions", set(EXPECTED_MAJOR_SPECIES))
        if abs(sum(_finite(value, f"reference.cases[{index}].initial.{name}") for name, value in initial.items()) - 1.0) > 1.0e-10:
            raise ContractError("reference-initial-mass-fraction-sum:" + case_id)
        samples = _array(case["samples"], f"reference.cases[{index}].samples")
        if len(samples) != len(EXPECTED_SAMPLE_TIMES):
            raise ContractError("reference-sample-count-mismatch:" + case_id)
        for sample_index, sample in enumerate(samples):
            item = _object(sample, f"reference.cases[{index}].samples[{sample_index}]", {"time_s", "temperature_k", "pressure_pa", "density_kg_m3", "heat_release_rate_w_m3", "major_species_mass_fractions"})
            _exact(_finite(item["time_s"], "reference.sample.time_s"), EXPECTED_SAMPLE_TIMES[sample_index], "reference.sample.time_s")
            if _finite(item["temperature_k"], "reference.sample.temperature_k") <= 0.0:
                raise ContractError("reference-sample-temperature-invalid:" + case_id)
            if _finite(item["pressure_pa"], "reference.sample.pressure_pa") <= 0.0:
                raise ContractError("reference-sample-pressure-invalid:" + case_id)
            if _finite(item["density_kg_m3"], "reference.sample.density_kg_m3") <= 0.0:
                raise ContractError("reference-sample-density-invalid:" + case_id)
            _finite(item["heat_release_rate_w_m3"], "reference.sample.heat_release_rate_w_m3")
            fractions = _object(item["major_species_mass_fractions"], "reference.sample.major_species_mass_fractions", set(EXPECTED_MAJOR_SPECIES))
            for name, value in fractions.items():
                if _finite(value, "reference.sample.mass_fraction." + name) < -1.0e-14:
                    raise ContractError("reference-sample-negative-mass-fraction:" + case_id)


def _verify_protocol(document: dict[str, Any]) -> None:
    local_reference = _object(document["local_reference"], "local_reference", {"mechanism", "runtime", "reference_output_status"})
    mechanism = _object(local_reference["mechanism"], "local_reference.mechanism", {"id", "status", "phase_model", "fit_use_allowed", "license_status", "excluded_domains"})
    _exact(mechanism["id"], "gri-mech-3.0-local-protocol-reference", "local_reference.mechanism.id")
    _exact(mechanism["status"], "local-reference-only-not-selected-for-fit", "local_reference.mechanism.status")
    _exact(mechanism["phase_model"], "ideal-gas", "local_reference.mechanism.phase_model")
    if _bool(mechanism["fit_use_allowed"], "local_reference.mechanism.fit_use_allowed"):
        raise ContractError("local-reference-must-not-be-fit-eligible")
    _exact(mechanism["license_status"], "not-reviewed-by-this-contract", "local_reference.mechanism.license_status")
    _exact(mechanism["excluded_domains"], ["cryogenic-injection", "vapor-liquid-equilibrium", "near-critical-real-fluid", "solid-risk", "two-phase-quality"], "local_reference.mechanism.excluded_domains")
    runtime = _object(local_reference["runtime"], "local_reference.runtime", {"cantera_version", "python_executable", "python_sha256", "cantera_module", "cantera_module_sha256", "cantera_extension", "cantera_extension_sha256", "execution_scope", "network_access", "fluent_access"})
    _exact(runtime["cantera_version"], "3.1.0", "local_reference.runtime.cantera_version")
    _exact(runtime["python_executable"], "/usr/bin/python3.8", "local_reference.runtime.python_executable")
    _exact(runtime["python_sha256"], "298a9e830ed52f36c299427565485d717d1ce0179c0597cc16560513eb780b06", "local_reference.runtime.python_sha256")
    _exact(runtime["cantera_module"], "/home/nvoyager/.local/lib/python3.8/site-packages/cantera/__init__.py", "local_reference.runtime.cantera_module")
    _exact(runtime["cantera_module_sha256"], "9a7ceae3bae1a7e41dc581e353d44e04af8447114748f5b7c75995e26f8819b9", "local_reference.runtime.cantera_module_sha256")
    _exact(runtime["cantera_extension"], "/home/nvoyager/.local/lib/python3.8/site-packages/cantera/_cantera.cpython-38-x86_64-linux-gnu.so", "local_reference.runtime.cantera_extension")
    _exact(runtime["cantera_extension_sha256"], "c57a47c201b0cb8f965b38acc7ebfafbbe95db7464d80359499c64d50555fae8", "local_reference.runtime.cantera_extension_sha256")
    _exact(runtime["execution_scope"], "local-offline-reference-generation-only", "local_reference.runtime.execution_scope")
    if _bool(runtime["network_access"], "local_reference.runtime.network_access"):
        raise ContractError("local-reference-network-access-must-remain-false")
    if _bool(runtime["fluent_access"], "local_reference.runtime.fluent_access"):
        raise ContractError("local-reference-fluent-access-must-remain-false")
    _exact(local_reference["reference_output_status"], "generated-but-not-a-fit-or-physical-validation-target", "local_reference.reference_output_status")

    protocol = _object(document["protocol"], "protocol", {"reactor", "protocol_training_case_ids", "protocol_holdout_case_ids", "declared_observables", "common_basis_status", "rfam_0d_runner_status", "no_posthoc_partition_changes"})
    reactor = _object(protocol["reactor"], "protocol.reactor", {"model", "energy", "reactor_scope", "sampling", "sample_times_s"})
    _exact(reactor["model"], "IdealGasConstPressureReactor", "protocol.reactor.model")
    _exact(reactor["energy"], "on", "protocol.reactor.energy")
    _exact(reactor["reactor_scope"], "high-temperature-gas-only", "protocol.reactor.reactor_scope")
    _exact(reactor["sampling"], "fixed-output-times-with-adaptive-internal-integrator", "protocol.reactor.sampling")
    _exact([_finite(value, "protocol.reactor.sample_times_s") for value in _array(reactor["sample_times_s"], "protocol.reactor.sample_times_s")], EXPECTED_SAMPLE_TIMES, "protocol.reactor.sample_times_s")
    training = _array(protocol["protocol_training_case_ids"], "protocol.protocol_training_case_ids")
    holdout = _array(protocol["protocol_holdout_case_ids"], "protocol.protocol_holdout_case_ids")
    _exact(training, [case_id for case_id, partition in EXPECTED_CASES if partition == "protocol-training"], "protocol.protocol_training_case_ids")
    _exact(holdout, [case_id for case_id, partition in EXPECTED_CASES if partition == "protocol-holdout"], "protocol.protocol_holdout_case_ids")
    if set(training) & set(holdout):
        raise ContractError("protocol-partition-overlap")
    _exact(protocol["declared_observables"], EXPECTED_OBSERVABLES, "protocol.declared_observables")
    _exact(protocol["common_basis_status"], "not-established-for-rfam-vs-cantera-trajectory-comparison", "protocol.common_basis_status")
    _exact(protocol["rfam_0d_runner_status"], "not-implemented", "protocol.rfam_0d_runner_status")
    if not _bool(protocol["no_posthoc_partition_changes"], "protocol.no_posthoc_partition_changes"):
        raise ContractError("protocol-no-posthoc-partition-change-required")


def _verify_gates(document: dict[str, Any]) -> None:
    gates = _object(document["fit_and_experiment_gates"], "fit_and_experiment_gates", {"status", "detailed_mechanism", "experimental_training_dataset", "experimental_holdout_dataset", "metrics_and_tolerances", "fit_parameters"})
    _exact(gates["status"], "blocked", "fit_and_experiment_gates.status")
    expected = {
        "detailed_mechanism": ("not-selected", ["explicit fit and redistribution terms", "versioned source and independent mechanism review", "declared gas-phase applicability envelope"]),
        "experimental_training_dataset": ("unselected", ["machine-readable raw data with immutable source hash", "explicit permission for fitting and redistribution", "conditions, units, uncertainty and protocol metadata", "predeclared mapping to protocol training cases"]),
        "experimental_holdout_dataset": ("unselected", ["machine-readable raw data with immutable source hash", "explicit permission for validation and redistribution", "conditions, units, uncertainty and protocol metadata", "independent and disjoint partition frozen before fitting"]),
        "metrics_and_tolerances": ("not-selected", ["predeclared ignition-delay metric", "predeclared temperature, major-species and heat-release trajectory metrics", "predeclared objective weights and tolerances", "common rate, thermodynamic and control-volume basis"]),
    }
    for key, (status, requirements) in expected.items():
        item = _object(gates[key], "fit_and_experiment_gates." + key, {"status", "required_before_fit"})
        _exact(item["status"], status, "fit_and_experiment_gates." + key + ".status")
        _exact(item["required_before_fit"], requirements, "fit_and_experiment_gates." + key + ".required_before_fit")
    fit_parameters = _object(gates["fit_parameters"], "fit_and_experiment_gates.fit_parameters", {"status", "allowed_after_all_gates"})
    _exact(fit_parameters["status"], "not-run", "fit_and_experiment_gates.fit_parameters.status")
    _exact(fit_parameters["allowed_after_all_gates"], ["R1 pre-exponential factor", "R1 activation energy", "R1 reaction orders", "R2 pre-exponential factor", "R2 activation energy", "R2 reaction orders", "formation-energy basis"], "fit_and_experiment_gates.fit_parameters.allowed_after_all_gates")


def validate(root: Path, contract_path: Path, verify_reference: bool, cantera_python: str) -> int:
    root = root.resolve()
    contract = _load_json(contract_path)
    document = _object(contract, "root", {"schema", "contract_version", "status", "scope", "approvals", "upstream", "artifacts", "local_reference", "protocol", "fit_and_experiment_gates", "cross_scope", "stop_conditions"})
    _exact(document["schema"], "rfam.r2-0d-calibration-contract/v1", "root.schema")
    if isinstance(document["contract_version"], bool) or document["contract_version"] != 1:
        raise ContractError("invalid-int:root.contract_version")
    _exact(document["status"], "blocked", "root.status")
    _exact(document["scope"], "local-cantera-0d-protocol-only-no-fitting-or-physical-validation", "root.scope")
    approvals = _object(document["approvals"], "approvals", set(EXPECTED_APPROVALS))
    for key, expected in EXPECTED_APPROVALS.items():
        if _bool(approvals[key], "approvals." + key) != expected:
            raise ContractError("approval-must-remain-false:" + key)
    _verify_upstream(root, document["upstream"])
    artifacts = _array(document["artifacts"], "artifacts")
    if len(artifacts) != len(EXPECTED_ARTIFACTS):
        raise ContractError("artifact-count-mismatch")
    artifact_paths = {
        expected["id"]: _verify_file_binding(root, artifacts[index], expected, f"artifacts[{index}]")
        for index, expected in enumerate(EXPECTED_ARTIFACTS)
    }
    _verify_reference(root, artifact_paths)
    _verify_protocol(document)
    _verify_gates(document)
    _exact(document["cross_scope"], {"r2_implementation_equivalence": "required-upstream-not-physical-validation", "kij": "not-a-cantera-calibration-input; sensitive-pair-policy-remains-fail-closed", "hem": "rfam-rho-e-y-is-the-only-two-phase-state-owner", "condensed_carbon": "wiring-or-diagnostic-only-no-kinetics-validation", "fluent": "not-invoked-and-not-approved"}, "cross_scope")
    _exact(document["stop_conditions"], ["any approval flag true", "local reference mechanism marked fit-eligible or selected", "changed local artifact identity or Cantera version", "overlapping, missing or reordered protocol partitions", "missing experimental rights, raw data, uncertainty, metric or holdout gate", "attempt to treat ideal-gas local trajectories as VLE/HEM evidence", "attempt to treat local reference output as physical R2 validation"], "stop_conditions")

    if verify_reference:
        if cantera_python != "/usr/bin/python3":
            raise ContractError("cantera-python-command-contract-mismatch")
        with tempfile.TemporaryDirectory(prefix="rfam-r2-0d-") as directory:
            sandbox = Path(directory)
            generator = sandbox / "generate_reference.py"
            mechanism = sandbox / "cantera_gri30_v3_0.yaml"
            output = sandbox / "regenerated.json"
            # Execute only copies of the bytes just bound above. This prevents
            # a mutable worktree path from changing between identity check and
            # subprocess launch.
            generator.write_bytes(artifact_paths["reference_generator"].read_bytes())
            mechanism.write_bytes(artifact_paths["local_cantera_mechanism"].read_bytes())
            result = subprocess.run(
                [cantera_python, str(generator), "--mechanism", str(mechanism), "--output", str(output)],
                cwd=root,
                text=True,
                capture_output=True,
                check=False,
            )
            if result.returncode != 0:
                raise ContractError("reference-regeneration-command-failed:" + result.stderr.strip())
            if not output.is_file() or output.read_bytes() != artifact_paths["reference_output"].read_bytes():
                raise ContractError("reference-regeneration-byte-mismatch")
    return 1 + len(EXPECTED_ARTIFACTS) + len(EXPECTED_CASES)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--contract", type=Path)
    parser.add_argument("--verify-reference", action="store_true")
    parser.add_argument("--cantera-python", default=sys.executable)
    args = parser.parse_args()
    contract_path = args.contract or args.root / "tests/reference/r2_0d_calibration_contract_v1.json"
    if not contract_path.is_absolute():
        contract_path = args.root / contract_path
    try:
        checks = validate(args.root, contract_path, args.verify_reference, args.cantera_python)
    except ContractError as error:
        print("R2_0D_CALIBRATION_CONTRACT_FAIL " + str(error), file=sys.stderr)
        return 2
    except Exception as error:  # Must never let malformed contract data escape as a traceback.
        print("R2_0D_CALIBRATION_CONTRACT_FAIL unexpected:" + type(error).__name__, file=sys.stderr)
        return 2
    print("R2_0D_CALIBRATION_CONTRACT_BLOCKED checks=" + str(checks))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
