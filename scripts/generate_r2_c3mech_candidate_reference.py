#!/usr/bin/env python3
"""Generate one partition of the bounded C3Mech R2 trajectory screen."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Any

import cantera as ct


class ReferenceError(RuntimeError):
    pass


SPECIES = ["C3H8", "O2", "CO", "CO2", "H2O", "N2"]
RUNTIME = {
    "python_executable": "C:\\Users\\guhai\\AppData\\Local\\Programs\\Python\\Python312\\python.exe",
    "python_sha256": "4d6f5f81a4bca11191c4c7c6b43632694d0a4ce74e068619d8fdc161d469859a",
    "python_size": 104952,
    "cantera_version": "3.2.0",
    "cantera_wheel": "F:\\tmp\\cantera-3.2.0-cp312-cp312-win_amd64.whl",
    "cantera_wheel_sha256": "258b95a4e29f9afd50788eab284e88440174302123c8b3a4a3a493904b40e238",
    "cantera_wheel_size": 5291378,
    "cantera_module": "F:\\tmp\\cantera-3.2.0-python312\\cantera\\__init__.py",
    "cantera_module_sha256": "a793b6d7643cc875ef0e0bbe5e3051f538af08418242aa87eb6cc7f995e34396",
    "cantera_module_size": 1014,
    "cantera_extension": "F:\\tmp\\cantera-3.2.0-python312\\cantera\\_cantera.cp312-win_amd64.pyd",
    "cantera_extension_sha256": "00f00749cc90b76a7266656fd02e6a7d40be82e7d7e1389472ad285123064e9d",
    "cantera_extension_size": 6564352,
    "execution_authority": "bound-installed-module-and-extension",
    "wheel_to_installed_tree_equivalence_proven": False,
}
PARTITIONS = {
    "protocol-training": {
        "manifest_sha256": "2baf34d21180d9b81c2abe841a00ae8bf263f15cd0abbf949ac9ba36f92d6601",
        "case_ids": ["train_phi0p8_diluted", "train_phi1p0", "train_phi1p2"],
    },
    "protocol-holdout": {
        "manifest_sha256": "cb63bfd73fdf7cc56eb179ee95d84be185c7397406f25b93183588bb8b4d4785",
        "case_ids": ["holdout_phi0p6_diluted", "holdout_phi1p1", "holdout_phi1p4"],
    },
}
MODELS = {
    "44I0": {
        "mechanism_sha256": "d5cec179ffbb6221b1b91264e1439ec850658dd97a0298054c8dbc336f3ab148",
        "mechanism_size": 2227356,
        "species_count": 729,
        "header_reaction_count": 3964,
        "yaml_entry_count": 4080,
        "thermo_min_k": 400.0,
        "thermo_max_k": 2000.0,
        "role": "primary-lt-ht-reference-candidate-not-selected",
    },
    "3R00": {
        "mechanism_sha256": "e307306e8fe9497ebc5a88f32545936c6f65a12c131c90ec30cff18fdf3a3e9d",
        "mechanism_size": 1483615,
        "species_count": 378,
        "header_reaction_count": 2521,
        "yaml_entry_count": 2626,
        "thermo_min_k": 300.0,
        "thermo_max_k": 2000.0,
        "role": "bounded-ht-comparator-not-selected",
    },
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ReferenceError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        document = json.loads(
            path.read_text(encoding="utf-8"), object_pairs_hook=strict_object
        )
    except ReferenceError:
        raise
    except Exception as exc:
        raise ReferenceError(f"invalid-json:{label}:{type(exc).__name__}") from exc
    if not isinstance(document, dict):
        raise ReferenceError(f"not-an-object:{label}")
    return document


def finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ReferenceError(f"not-a-number:{label}")
    result = float(value)
    if not math.isfinite(result):
        raise ReferenceError(f"nonfinite:{label}")
    return result


def stable(value: float, label: str) -> float:
    if not math.isfinite(value):
        raise ReferenceError(f"nonfinite:{label}")
    return float(format(float(value), ".12g"))


def stable_event_time(value: float, label: str) -> float:
    if not math.isfinite(value):
        raise ReferenceError(f"nonfinite:{label}")
    return float(format(float(value), ".8g"))


def verify_runtime() -> dict[str, Any]:
    wheel = Path(RUNTIME["cantera_wheel"]).resolve()
    actual = {
        "python_executable": str(Path(sys.executable).resolve()),
        "python_sha256": sha256(Path(sys.executable).resolve()),
        "python_size": Path(sys.executable).resolve().stat().st_size,
        "cantera_version": ct.__version__,
        "cantera_wheel": str(wheel),
        "cantera_wheel_sha256": sha256(wheel),
        "cantera_wheel_size": wheel.stat().st_size,
        "cantera_module": str(Path(ct.__file__).resolve()),
        "cantera_module_sha256": sha256(Path(ct.__file__).resolve()),
        "cantera_module_size": Path(ct.__file__).resolve().stat().st_size,
        "cantera_extension": str(Path(ct._cantera.__file__).resolve()),
        "cantera_extension_sha256": sha256(Path(ct._cantera.__file__).resolve()),
        "cantera_extension_size": Path(ct._cantera.__file__).resolve().stat().st_size,
        "execution_authority": "bound-installed-module-and-extension",
        "wheel_to_installed_tree_equivalence_proven": False,
    }
    if actual != RUNTIME:
        raise ReferenceError("runtime-identity-mismatch")
    return actual


def verify_manifest(path: Path, requested_partition: str) -> dict[str, Any]:
    expected = PARTITIONS.get(requested_partition)
    if expected is None:
        raise ReferenceError("unknown-partition")
    if sha256(path) != expected["manifest_sha256"]:
        raise ReferenceError("case-manifest-sha256-mismatch")
    document = load_json(path, "case-manifest")
    if set(document) != {"schema", "partition", "partition_locked", "sample_times_s", "cases"}:
        raise ReferenceError("case-manifest-schema-open")
    if document["schema"] != "rfam.r2-c3mech-case-partition/v1":
        raise ReferenceError("case-manifest-schema-mismatch")
    if document["partition"] != requested_partition or document["partition_locked"] is not True:
        raise ReferenceError("case-manifest-partition-mismatch")
    times = document["sample_times_s"]
    if not isinstance(times, list) or len(times) != 14:
        raise ReferenceError("sample-time-count-mismatch")
    parsed_times = [finite(value, f"sample-time[{index}]") for index, value in enumerate(times)]
    if parsed_times[0] != 0.0 or any(
        right <= left for left, right in zip(parsed_times, parsed_times[1:])
    ):
        raise ReferenceError("sample-times-not-strictly-increasing")
    cases = document["cases"]
    if not isinstance(cases, list) or [case.get("id") for case in cases if isinstance(case, dict)] != expected["case_ids"]:
        raise ReferenceError("case-order-mismatch")
    allowed = {"id", "temperature_k", "pressure_pa", "equivalence_ratio", "fuel", "oxidizer"}
    for case in cases:
        if not isinstance(case, dict) or set(case) != allowed:
            raise ReferenceError("case-schema-mismatch")
        if finite(case["temperature_k"], f"{case['id']}.temperature") <= 0.0:
            raise ReferenceError("case-temperature-invalid")
        if finite(case["pressure_pa"], f"{case['id']}.pressure") <= 0.0:
            raise ReferenceError("case-pressure-invalid")
        if finite(case["equivalence_ratio"], f"{case['id']}.phi") <= 0.0:
            raise ReferenceError("case-equivalence-ratio-invalid")
        if not isinstance(case["fuel"], str) or not isinstance(case["oxidizer"], str):
            raise ReferenceError("case-composition-invalid")
    return document


def parse_mechanism_header(path: Path) -> tuple[str, int, int]:
    text = path.read_text(encoding="utf-8")
    model_match = re.search(r"model ID \(MID\):\s*([A-Za-z0-9]+)", text)
    count_match = re.search(r"number of reactions:\s*([0-9]+)", text)
    yaml_entries = len(re.findall(r"(?m)^- equation:", text))
    if model_match is None or count_match is None:
        raise ReferenceError("mechanism-header-metadata-missing")
    return model_match.group(1), int(count_match.group(1)), yaml_entries


def load_mechanism(path: Path, model_id: str) -> tuple[ct.Solution, dict[str, Any]]:
    expected = MODELS.get(model_id)
    if expected is None:
        raise ReferenceError("unknown-model-id")
    if not path.is_file() or path.stat().st_size != expected["mechanism_size"]:
        raise ReferenceError("mechanism-size-mismatch")
    if sha256(path) != expected["mechanism_sha256"]:
        raise ReferenceError("mechanism-sha256-mismatch")
    header_model, header_count, yaml_entries = parse_mechanism_header(path)
    if header_model != model_id:
        raise ReferenceError("mechanism-header-model-id-mismatch")
    if header_count != expected["header_reaction_count"]:
        raise ReferenceError("mechanism-header-reaction-count-mismatch")
    if yaml_entries != expected["yaml_entry_count"]:
        raise ReferenceError("mechanism-yaml-entry-count-mismatch")
    gas = ct.Solution(str(path))
    missing = [name for name in SPECIES if name not in gas.species_names]
    if missing:
        raise ReferenceError("mechanism-required-species-missing:" + ",".join(missing))
    actual = {
        "model_id": model_id,
        "role": expected["role"],
        "species_count": gas.n_species,
        "author_header_reaction_count": header_count,
        "yaml_reaction_entry_count": yaml_entries,
        "cantera_loaded_reaction_count": gas.n_reactions,
        "thermo_min_k": stable(float(gas.min_temp), "thermo-min"),
        "thermo_max_k": stable(float(gas.max_temp), "thermo-max"),
        "counting_scope": "header count is provenance text; YAML entries and Cantera loaded count are execution authority",
        "count_discrepancy_resolved_as_equivalence": False,
    }
    required = {
        "species_count": expected["species_count"],
        "cantera_loaded_reaction_count": expected["yaml_entry_count"],
        "thermo_min_k": expected["thermo_min_k"],
        "thermo_max_k": expected["thermo_max_k"],
    }
    for key, value in required.items():
        if actual[key] != value:
            raise ReferenceError(f"mechanism-runtime-mismatch:{key}")
    if gas.thermo_model != "ideal-gas":
        raise ReferenceError("mechanism-not-ideal-gas")
    return gas, actual


def configure(gas: ct.Solution, case: dict[str, Any]) -> None:
    gas.TP = float(case["temperature_k"]), float(case["pressure_pa"])
    gas.set_equivalence_ratio(
        float(case["equivalence_ratio"]), case["fuel"], case["oxidizer"]
    )


def major_mass_fractions(gas: ct.Solution) -> dict[str, float]:
    result = {name: stable(float(gas[name].Y[0]), f"Y.{name}") for name in SPECIES}
    if any(value < 0.0 or value > 1.0 for value in result.values()):
        raise ReferenceError("major-species-mass-fraction-out-of-range")
    return result


def state_record(time_s: float, reactor: ct.IdealGasConstPressureReactor) -> dict[str, Any]:
    phase = reactor.phase
    return {
        "time_s": stable(time_s, "sample-time"),
        "status": "eligible-within-mechanism-thermo-intersection",
        "temperature_k": stable(float(phase.T), "sample-temperature"),
        "pressure_pa": stable(float(phase.P), "sample-pressure"),
        "density_kg_m3": stable(float(phase.density), "sample-density"),
        "heat_release_rate_w_m3": stable(float(phase.heat_release_rate), "sample-qdot"),
        "major_species_mass_fractions": major_mass_fractions(phase),
    }


def trajectory(
    gas: ct.Solution,
    case: dict[str, Any],
    sample_times: list[float],
    thermo_min: float,
    thermo_max: float,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    configure(gas, case)
    reactor = ct.IdealGasConstPressureReactor(gas, energy="on", clone=False)
    network = ct.ReactorNet([reactor])
    network.rtol = 1.0e-10
    network.atol = 1.0e-18
    reactor.set_advance_limit("temperature", 10.0)
    samples = [state_record(0.0, reactor)]
    exit_record: dict[str, Any] | None = None
    for target in sample_times[1:]:
        if exit_record is not None:
            samples.append({
                "time_s": stable(target, "sample-time"),
                "status": "not-evaluated-after-thermo-envelope-exit",
            })
            continue
        while network.time < target:
            previous = float(network.time)
            previous_temperature = float(reactor.T)
            reached = float(network.advance(target, apply_limit=True))
            if not math.isfinite(reached) or reached <= previous or reached > target:
                raise ReferenceError("integrator-time-progress-invalid")
            temperature = float(reactor.T)
            if not math.isfinite(temperature):
                raise ReferenceError("integrator-temperature-nonfinite")
            if temperature < thermo_min or temperature > thermo_max:
                if not (previous_temperature <= thermo_max < temperature):
                    raise ReferenceError("thermo-envelope-exit-bracket-invalid")
                fraction = (thermo_max - previous_temperature) / (
                    temperature - previous_temperature
                )
                crossing_time = previous + fraction * (reached - previous)
                exit_record = {
                    "requested_sample_time_s": stable(target, "exit-request-time"),
                    "thermo_min_k": stable(thermo_min, "thermo-min"),
                    "thermo_max_k": stable(thermo_max, "thermo-max"),
                    "first_crossing_estimate": {
                        "time_s": stable_event_time(crossing_time, "crossing-time"),
                        "temperature_k": stable(thermo_max, "crossing-temperature"),
                        "pressure_pa": stable(
                            float(reactor.phase.P), "crossing-pressure"
                        ),
                        "method": "10-K-advance-limited-linear-temperature-interpolation",
                        "publishable_chemistry_state": False,
                    },
                    "status": "thermo-envelope-exit-detected-no-endpoint-state-published",
                }
                samples.append({
                    "time_s": stable(target, "sample-time"),
                    "status": "outside-thermo-envelope-unscored",
                    "exit_detection": exit_record,
                })
                break
        else:
            samples.append(state_record(target, reactor))
    if exit_record is None:
        exit_record = {"status": "no-exit-detected-within-protocol-window"}
    return samples, exit_record


def ignition_delay(
    gas: ct.Solution,
    case: dict[str, Any],
    end_time: float,
    thermo_min: float,
    thermo_max: float,
) -> dict[str, Any]:
    configure(gas, case)
    initial_temperature = float(gas.T)
    threshold_temperature = initial_temperature + 400.0
    if threshold_temperature >= thermo_max:
        return {
            "status": "threshold-at-or-outside-thermo-upper-bound-unscored",
            "threshold_temperature_k": stable(threshold_temperature, "ignition-threshold"),
            "thermo_max_k": stable(thermo_max, "thermo-max"),
        }
    reactor = ct.IdealGasConstPressureReactor(gas, energy="on", clone=False)
    network = ct.ReactorNet([reactor])
    network.rtol = 1.0e-10
    network.atol = 1.0e-18
    reactor.set_advance_limit("temperature", 1.0)
    previous_time = float(network.time)
    previous_temperature = float(reactor.T)
    while network.time < end_time:
        reached = float(network.advance(end_time, apply_limit=True))
        temperature = float(reactor.T)
        if not math.isfinite(reached) or not math.isfinite(temperature):
            raise ReferenceError("ignition-integrator-nonfinite")
        if temperature < thermo_min or temperature > thermo_max:
            return {
                "status": "outside-thermo-envelope-unscored",
                "threshold_temperature_k": stable(threshold_temperature, "ignition-threshold"),
            }
        if temperature >= threshold_temperature:
            delta_temperature = temperature - previous_temperature
            if delta_temperature <= 0.0:
                raise ReferenceError("ignition-threshold-bracket-invalid")
            fraction = (threshold_temperature - previous_temperature) / delta_temperature
            if fraction < 0.0 or fraction > 1.0:
                raise ReferenceError("ignition-threshold-fraction-invalid")
            interpolated_time = previous_time + fraction * (reached - previous_time)
            return {
                "status": "eligible-within-mechanism-thermo-intersection",
                "ignition_delay_s": stable_event_time(interpolated_time, "ignition-delay"),
                "threshold_temperature_k": stable(threshold_temperature, "ignition-threshold"),
                "method": "1-K-advance-limited-linear-threshold-interpolation",
            }
        previous_time = reached
        previous_temperature = temperature
    return {
        "status": "not-detected-within-protocol-window",
        "end_time_s": stable(end_time, "ignition-end-time"),
        "threshold_temperature_k": stable(threshold_temperature, "ignition-threshold"),
    }


def run_case(
    gas: ct.Solution,
    case: dict[str, Any],
    sample_times: list[float],
    thermo_min: float,
    thermo_max: float,
) -> dict[str, Any]:
    samples, exit_record = trajectory(
        gas, case, sample_times, thermo_min, thermo_max
    )
    idt = ignition_delay(gas, case, sample_times[-1], thermo_min, thermo_max)
    initial = samples[0]["major_species_mass_fractions"]
    return {
        **case,
        "partition": "protocol-training" if case["id"].startswith("train_") else "protocol-holdout",
        "full_protocol_status": (
            "OUT_OF_MECHANISM_THERMO_ENVELOPE"
            if exit_record["status"]
            == "thermo-envelope-exit-detected-no-endpoint-state-published"
            else "WITHIN_MECHANISM_THERMO_ENVELOPE"
        ),
        "initial_major_species_mass_fractions": initial,
        "ignition_delay_definition": "1 K advance-limited linear interpolation to T - T0 = 400 K; threshold at thermo upper bound is unscored",
        "ignition_delay": idt,
        "thermo_envelope_exit": exit_record,
        "samples": samples,
    }


def generate(mechanism: Path, model_id: str, manifest_path: Path, partition: str) -> dict[str, Any]:
    runtime = verify_runtime()
    manifest = verify_manifest(manifest_path, partition)
    gas, inventory = load_mechanism(mechanism, model_id)
    if hasattr(ct, "suppress_thermo_warnings"):
        ct.suppress_thermo_warnings()
    sample_times = [float(value) for value in manifest["sample_times_s"]]
    cases = [
        run_case(
            gas,
            case,
            sample_times,
            float(inventory["thermo_min_k"]),
            float(inventory["thermo_max_k"]),
        )
        for case in manifest["cases"]
    ]
    any_exit = any(
        case["thermo_envelope_exit"]["status"]
        == "thermo-envelope-exit-detected-no-endpoint-state-published"
        for case in cases
    )
    return {
        "schema": "rfam.r2-c3mech-cantera-reference/v2",
        "status": "OUT_OF_MECHANISM_THERMO_ENVELOPE" if any_exit else "WITHIN_MECHANISM_THERMO_ENVELOPE",
        "model": inventory,
        "partition": partition,
        "runtime_identity": runtime,
        "source": {
            "mechanism_path": str(mechanism),
            "mechanism_sha256": sha256(mechanism),
            "mechanism_size": mechanism.stat().st_size,
            "generator_sha256": sha256(Path(__file__).resolve()),
            "case_manifest_path": str(manifest_path),
            "case_manifest_sha256": sha256(manifest_path),
        },
        "reactor": {
            "model": "IdealGasConstPressureReactor",
            "phase": "ideal-gas",
            "energy": "on",
            "rtol": 1.0e-10,
            "atol": 1.0e-18,
            "temperature_advance_limit_k": 10.0,
            "ignition_temperature_advance_limit_k": 1.0,
            "outside_envelope_policy": "stop-at-first-detected-exit; publish-no-state; leave-later-samples-unevaluated",
        },
        "sample_times_s": [stable(value, "sample-time") for value in sample_times],
        "major_species": SPECIES,
        "prefix_diagnostic": {
            "eligible_samples_only": True,
            "outside_endpoint_state_published": False,
            "later_samples_evaluated": False,
        },
        "cases": cases,
        "approvals": {
            "mechanism_selected_for_fit": False,
            "fit_executed": False,
            "holdout_used_for_selection": False,
            "physical_validation_approved": False,
            "production_ready": False,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mechanism", type=Path, required=True)
    parser.add_argument("--model-id", choices=sorted(MODELS), required=True)
    parser.add_argument("--case-manifest", type=Path, required=True)
    parser.add_argument("--partition", choices=sorted(PARTITIONS), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    output = args.output.resolve()
    temporary = output.with_name(output.name + ".tmp")
    try:
        document = generate(
            args.mechanism.resolve(),
            args.model_id,
            args.case_manifest.resolve(),
            args.partition,
        )
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary.write_text(
            json.dumps(document, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, output)
    except Exception as exc:
        try:
            temporary.unlink(missing_ok=True)
        except OSError:
            pass
        print(
            "R2_C3MECH_REFERENCE_FAIL reason=" + str(exc).replace("\n", " | "),
            file=sys.stderr,
        )
        return 2
    print(
        "R2_C3MECH_REFERENCE_BLOCKED "
        f"model={args.model_id} partition={args.partition} cases={len(document['cases'])}"
    )
    print("fit_executed=0 physical_validation_approved=0 production_ready=0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
