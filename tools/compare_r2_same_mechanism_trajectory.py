#!/usr/bin/env python3
"""Compare RFAM and Cantera trajectories for the same TwoStepCO mechanism.

This is an implementation-equivalence gate. It does not fit the mechanism or
validate the physical R2 closure.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

import cantera as ct


class ComparisonError(RuntimeError):
    pass


RUNNER_HEADER_KEYS = {
    "schema",
    "status",
    "case_id",
    "reactor",
    "chemistry",
    "energy_basis",
    "max_internal_step_s",
    "solver_tolerance_profile",
    "solver_relative_tolerance",
    "solver_species_absolute_tolerance",
    "solver_enthalpy_absolute_tolerance",
    "species_order",
    "sample_count",
    "ignition_definition",
    "ignition_detected",
    "ignition_delay_s",
    "positivity_strategy",
    "depletion_event_count",
    "interval_bisection_count",
    "roundoff_clamp_count",
    "accepted_physical_time_s",
    "implicit_time_extension_s",
    "fit_approved",
    "physical_validation_approved",
    "production_ready",
}

RUNNER_SAMPLE_KEYS = {
    "sample",
    "time_s",
    "temperature_k",
    "pressure_pa",
    "density_kg_m3",
    "sensible_enthalpy_j_kg",
    "total_enthalpy_j_kg",
    "qdot_w_m3",
    "temperature_rate_k_s",
    "mass_fraction_rates_s_inv",
    "mass_fractions",
}


def strict_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ComparisonError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"), object_pairs_hook=strict_object)
    except ComparisonError:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ComparisonError(f"json-load-failed:{label}:{error}") from error
    if not isinstance(value, dict):
        raise ComparisonError(f"json-root-not-object:{label}")
    return value


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise ComparisonError(f"artifact-open-failed:{path}:{error}") from error


def finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ComparisonError(f"invalid-number:{label}")
    try:
        result = float(value)
    except ValueError as error:
        raise ComparisonError(f"invalid-number:{label}") from error
    if not math.isfinite(result):
        raise ComparisonError(f"nonfinite-number:{label}")
    return result


def vector(text: str, length: int, label: str) -> list[float]:
    values = [finite(item, label) for item in text.split(",")]
    if len(values) != length:
        raise ComparisonError(f"vector-length-mismatch:{label}")
    return values


def require_close(
    actual: float,
    expected: float,
    relative: float,
    absolute: float,
    label: str,
) -> float:
    error = abs(actual - expected)
    if error > absolute + relative * abs(expected):
        raise ComparisonError(
            f"value-mismatch:{label}:actual={actual:.17g}:"
            f"expected={expected:.17g}:error={error:.17g}"
        )
    return error / max(absolute, abs(expected))


def validate_contract(root: Path, contract: Path) -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(root / "tools/validate_r2_same_mechanism_trajectory_contract.py"),
            "--root",
            str(root),
            "--contract",
            str(contract),
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    combined = result.stdout + result.stderr
    if (
        result.returncode != 0
        or "R2_SAME_MECHANISM_TRAJECTORY_CONTRACT_BLOCKED" not in result.stdout
        or "R2_SAME_MECHANISM_TRAJECTORY_CONTRACT_FAIL" in combined
    ):
        raise ComparisonError(
            "contract-validator-rejected:"
            + combined.strip().replace("\n", "|")
        )


def artifact_paths(root: Path, contract: dict[str, Any]) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for item in contract["artifacts"]:
        path = (root / item["path"]).resolve()
        try:
            path.relative_to(root)
        except ValueError as error:
            raise ComparisonError(f"artifact-path-escapes-root:{item['id']}") from error
        if not path.is_file() or sha256(path) != item["sha256"]:
            raise ComparisonError(f"artifact-identity-mismatch:{item['id']}")
        paths[item["id"]] = path
    return paths


def build_time_grid(contract: dict[str, Any]) -> list[float]:
    grid = contract["time_grid"]
    points = int(grid["log_points_including_endpoints"])
    lower = math.log10(float(grid["minimum_positive_time_s"]))
    upper = math.log10(float(grid["maximum_time_s"]))
    result = [0.0]
    for index in range(points):
        exponent = lower + (upper - lower) * index / (points - 1)
        result.append(10.0 ** exponent)
    if len(result) != int(grid["total_points_per_case"]):
        raise ComparisonError("time-grid-count-mismatch")
    if result[1] != float(grid["minimum_positive_time_s"]):
        raise ComparisonError("time-grid-lower-bound-mismatch")
    if result[-1] != float(grid["maximum_time_s"]):
        raise ComparisonError("time-grid-upper-bound-mismatch")
    return result


def parse_runner_output(
    stdout: str,
    stderr: str,
    case: dict[str, Any],
    sample_times: list[float],
    species: list[str],
    profile: dict[str, Any],
) -> tuple[dict[str, str], list[dict[str, Any]]]:
    if stderr.strip():
        raise ComparisonError(f"runner-stderr:{case['id']}:{stderr.strip()}")
    lines = stdout.splitlines()
    if not lines or lines[-1] != "R2_0D_RFAM_RUNNER_OK":
        raise ComparisonError(f"runner-success-marker-missing:{case['id']}")
    if lines.count("R2_0D_RFAM_RUNNER_OK") != 1:
        raise ComparisonError(f"runner-success-marker-count:{case['id']}")

    headers: dict[str, str] = {}
    sample_lines: list[str] = []
    for line in lines[:-1]:
        if line.startswith("sample="):
            sample_lines.append(line)
            continue
        if "=" not in line:
            raise ComparisonError(f"runner-header-invalid:{case['id']}")
        key, value = line.split("=", 1)
        if not key or key in headers:
            raise ComparisonError(f"runner-header-duplicate:{case['id']}:{key}")
        headers[key] = value
    if set(headers) != RUNNER_HEADER_KEYS:
        raise ComparisonError(f"runner-header-shape-mismatch:{case['id']}")

    expected_text = {
        "schema": "rfam.r2-0d-rfam-trajectory/v1",
        "status": "implementation-protocol-only-not-fit-eligible",
        "case_id": case["id"],
        "reactor": "ideal-gas-constant-pressure-sensible-enthalpy",
        "chemistry": "rfam-twostepco-candidate",
        "energy_basis": "sensible-enthalpy-plus-explicit-formation-heat",
        "solver_tolerance_profile": profile["name"],
        "species_order": ",".join(species),
        "ignition_definition": "first-root-T-minus-T0-equals-400K",
        "positivity_strategy": (
            "fuel-depletion-root-plus-rejected-interval-bisection"
        ),
        "implicit_time_extension_s": "0",
        "fit_approved": "0",
        "physical_validation_approved": "0",
        "production_ready": "0",
    }
    for key, expected in expected_text.items():
        if headers[key] != expected:
            raise ComparisonError(
                f"runner-header-value-mismatch:{case['id']}:{key}")
    numeric_headers = {
        "max_internal_step_s": profile["maximum_internal_step_s"],
        "solver_relative_tolerance": profile["relative_tolerance"],
        "solver_species_absolute_tolerance": profile[
            "species_absolute_tolerance"
        ],
        "solver_enthalpy_absolute_tolerance": profile[
            "sensible_enthalpy_absolute_tolerance_j_per_kg"
        ],
        "accepted_physical_time_s": sample_times[-1],
    }
    for key, expected in numeric_headers.items():
        if finite(headers[key], f"{case['id']}.{key}") != float(expected):
            raise ComparisonError(
                f"runner-header-number-mismatch:{case['id']}:{key}")
    if int(headers["sample_count"]) != len(sample_times):
        raise ComparisonError(f"runner-sample-count-mismatch:{case['id']}")
    if int(headers["ignition_detected"]) != 1:
        raise ComparisonError(f"runner-ignition-not-detected:{case['id']}")
    ignition_delay = finite(headers["ignition_delay_s"], "ignition-delay")
    if not 0.0 < ignition_delay < sample_times[-1]:
        raise ComparisonError(f"runner-ignition-delay-invalid:{case['id']}")
    for key in (
        "depletion_event_count",
        "interval_bisection_count",
        "roundoff_clamp_count",
    ):
        count = int(headers[key])
        if count < 0 or count > 1000000:
            raise ComparisonError(f"runner-counter-invalid:{case['id']}:{key}")

    if len(sample_lines) != len(sample_times):
        raise ComparisonError(f"runner-sample-lines-mismatch:{case['id']}")
    samples: list[dict[str, Any]] = []
    for index, (line, expected_time) in enumerate(zip(sample_lines, sample_times)):
        fields: dict[str, str] = {}
        for token in line.split("|"):
            if "=" not in token:
                raise ComparisonError(f"runner-sample-token-invalid:{case['id']}")
            key, value = token.split("=", 1)
            if not key or key in fields:
                raise ComparisonError(
                    f"runner-sample-key-invalid:{case['id']}:{key}")
            fields[key] = value
        if set(fields) != RUNNER_SAMPLE_KEYS:
            raise ComparisonError(f"runner-sample-shape-mismatch:{case['id']}")
        if int(fields["sample"]) != index:
            raise ComparisonError(f"runner-sample-index-mismatch:{case['id']}")
        time = finite(fields["time_s"], f"{case['id']}.time")
        if time != expected_time:
            raise ComparisonError(f"runner-sample-time-mismatch:{case['id']}")
        mass_fractions = vector(
            fields["mass_fractions"], len(species), f"{case['id']}.Y")
        mass_fraction_rates = vector(
            fields["mass_fraction_rates_s_inv"],
            len(species),
            f"{case['id']}.dYdt",
        )
        sample = {
            "time_s": time,
            "temperature_k": finite(
                fields["temperature_k"], f"{case['id']}.temperature"),
            "pressure_pa": finite(fields["pressure_pa"], f"{case['id']}.pressure"),
            "density_kg_m3": finite(
                fields["density_kg_m3"], f"{case['id']}.density"),
            "sensible_enthalpy_j_kg": finite(
                fields["sensible_enthalpy_j_kg"], f"{case['id']}.hsens"),
            "total_enthalpy_j_kg": finite(
                fields["total_enthalpy_j_kg"], f"{case['id']}.htotal"),
            "qdot_w_m3": finite(fields["qdot_w_m3"], f"{case['id']}.qdot"),
            "temperature_rate_k_s": finite(
                fields["temperature_rate_k_s"], f"{case['id']}.dTdt"),
            "mass_fraction_rates_s_inv": mass_fraction_rates,
            "mass_fractions": mass_fractions,
        }
        if sample["pressure_pa"] != float(case["pressure_pa"]):
            raise ComparisonError(f"runner-pressure-mismatch:{case['id']}")
        if abs(sum(mass_fractions) - 1.0) > 2.0e-10:
            raise ComparisonError(f"runner-composition-sum-mismatch:{case['id']}")
        samples.append(sample)
    return headers, samples


def run_rfam_case(
    runner: Path,
    case: dict[str, Any],
    sample_times: list[float],
    species: list[str],
    profile: dict[str, Any],
) -> tuple[dict[str, str], list[dict[str, Any]]]:
    command = [
        str(runner),
        case["id"],
        format(float(case["temperature_k"]), ".17g"),
        format(float(case["pressure_pa"]), ".17g"),
        str(len(sample_times)),
        format(float(profile["maximum_internal_step_s"]), ".17g"),
        *[format(float(value), ".17g") for value in case["mass_fractions"]],
        *[format(value, ".17g") for value in sample_times],
        "--solver-tolerance-profile=" + profile["name"],
    ]
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    if result.returncode != 0:
        raise ComparisonError(
            f"runner-rejected:{case['id']}:rc={result.returncode}:"
            f"{result.stdout.strip()}:{result.stderr.strip()}"
        )
    return parse_runner_output(
        result.stdout, result.stderr, case, sample_times, species, profile)


def reference_enthalpies(
    gas: ct.Solution,
    species: list[str],
    energy_reference: dict[str, Any],
) -> tuple[list[float], list[dict[str, Any]]]:
    reference_temperature = float(energy_reference["reference_temperature_k"])
    result: list[float] = []
    audit: list[dict[str, Any]] = []
    lower_tail_species: list[str] = []
    for name in species:
        thermo = gas.species(name).thermo
        minimum = finite(
            float(thermo.min_temp), f"reference-minimum-temperature:{name}")
        cantera_polynomial_reference = finite(
            float(thermo.h(reference_temperature)),
            f"cantera-reference-enthalpy:{name}",
        )
        if reference_temperature < minimum:
            lower_tail_species.append(name)
            boundary_enthalpy = finite(
                float(thermo.h(minimum)), f"boundary-enthalpy:{name}")
            boundary_heat_capacity = finite(
                float(thermo.cp(minimum)), f"boundary-heat-capacity:{name}")
            value = boundary_enthalpy + boundary_heat_capacity * (
                reference_temperature - minimum)
            branch = "lower-boundary-cp-linear-tail"
        else:
            value = cantera_polynomial_reference
            branch = "in-range-nasa7"
        if not math.isfinite(value):
            raise ComparisonError(f"reference-enthalpy-nonfinite:{name}")
        result.append(value)
        audit.append(
            {
                "species": name,
                "minimum_temperature_k": minimum,
                "reference_temperature_k": reference_temperature,
                "branch": branch,
                "cantera_polynomial_reference_j_per_kmol": (
                    cantera_polynomial_reference
                ),
                "rfam_reference_j_per_kmol": value,
                "rfam_minus_cantera_reference_j_per_kmol": (
                    value - cantera_polynomial_reference
                ),
            }
        )
    if lower_tail_species != energy_reference["lower_linear_tail_species"]:
        raise ComparisonError("reference-lower-tail-species-mismatch")
    return result, audit


def reaction_heats(instantaneous: dict[str, Any]) -> tuple[float, float]:
    try:
        r1 = float(
            instantaneous["reaction_contract"]["R1"][
                "DH_J_per_mol_progress"
            ]
        )
        r2 = float(
            instantaneous["reaction_contract"]["R2"][
                "DH_J_per_mol_progress"
            ]
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ComparisonError("reaction-heat-contract-invalid") from error
    if not math.isfinite(r1) or not math.isfinite(r2) or r1 >= 0.0 or r2 >= 0.0:
        raise ComparisonError("reaction-heat-contract-invalid")
    return r1, r2


def rfam_qdot(gas: ct.Solution, heats: tuple[float, float]) -> float:
    rates = gas.net_rates_of_progress
    if len(rates) != 2:
        raise ComparisonError("cantera-reaction-count-mismatch")
    rate_r1 = finite(float(rates[0]), "cantera-rate-of-progress:R1")
    rate_r2 = finite(float(rates[1]), "cantera-rate-of-progress:R2")
    return finite(
        rate_r1 * 1000.0 * (-heats[0])
        + rate_r2 * 1000.0 * (-heats[1]),
        "cantera-rfam-qdot",
    )


class RFAMEnergyReferenceReactor(ct.ExtensibleIdealGasConstPressureReactor):
    def after_eval(self, time_s, lhs, rhs):
        del time_s, lhs
        production = self.thermo.net_production_rates
        if len(production) != len(self.rfam_reference_enthalpies):
            raise ComparisonError("cantera-reference-rate-species-count-mismatch")
        reference_rate = finite(
            sum(
                finite(
                    float(rate), f"cantera-reference-production-rate:{index}")
                * finite(reference, f"rfam-reference-enthalpy:{index}")
                for index, (rate, reference) in enumerate(
                    zip(production, self.rfam_reference_enthalpies)
                )
            ),
            "cantera-reference-enthalpy-rate",
        )
        qdot = rfam_qdot(self.thermo, self.rfam_reaction_heats)
        rhs[1] = finite(
            float(rhs[1])
            + finite(float(self.volume), "cantera-reactor-volume")
            * (reference_rate + qdot),
            "cantera-energy-rhs",
        )


def make_reactor(
    gas: ct.Solution,
    references: list[float],
    heats: tuple[float, float],
) -> RFAMEnergyReferenceReactor:
    reactor = RFAMEnergyReferenceReactor(gas, energy="on")
    reactor.rfam_reference_enthalpies = list(references)
    reactor.rfam_reaction_heats = heats
    return reactor


def configure_network(
    reactor: RFAMEnergyReferenceReactor,
    profile: dict[str, Any],
) -> ct.ReactorNet:
    network = ct.ReactorNet([reactor])
    network.rtol = float(profile["relative_tolerance"])
    network.atol = float(profile["absolute_tolerance"])
    network.max_time_step = float(profile["maximum_time_step_s"])
    return network


def cantera_trajectory(
    mechanism: Path,
    case: dict[str, Any],
    sample_times: list[float],
    species: list[str],
    references: list[float],
    heats: tuple[float, float],
    profile: dict[str, Any],
) -> list[dict[str, Any]]:
    gas = ct.Solution(str(mechanism), "rfam-two-step")
    gas.TPY = (
        float(case["temperature_k"]),
        float(case["pressure_pa"]),
        case["mass_fractions"],
    )
    reactor = make_reactor(gas, references, heats)
    network = configure_network(reactor, profile)
    samples: list[dict[str, Any]] = []
    for index, time_s in enumerate(sample_times):
        if time_s > 0.0:
            network.advance(time_s)
        temperature = finite(
            float(reactor.T),
            f"cantera-trajectory:{case['id']}[{index}].temperature",
        )
        mass_fractions = [
            finite(
                float(reactor.thermo[name].Y[0]),
                f"cantera-trajectory:{case['id']}[{index}].Y[{name}]",
            )
            for name in species
        ]
        samples.append(
            {
                "time_s": time_s,
                "temperature_k": temperature,
                "mass_fractions": mass_fractions,
            }
        )
    return samples


def cantera_ignition_delay(
    mechanism: Path,
    case: dict[str, Any],
    references: list[float],
    heats: tuple[float, float],
    profile: dict[str, Any],
    maximum_time: float,
) -> float:
    gas = ct.Solution(str(mechanism), "rfam-two-step")
    gas.TPY = (
        float(case["temperature_k"]),
        float(case["pressure_pa"]),
        case["mass_fractions"],
    )
    reactor = make_reactor(gas, references, heats)
    network = configure_network(reactor, profile)
    threshold = finite(
        float(case["temperature_k"]) + 400.0,
        f"cantera-ignition-threshold:{case['id']}",
    )
    previous_time = 0.0
    previous_temperature = finite(
        float(reactor.T), f"cantera-ignition-initial-temperature:{case['id']}")
    steps = 0
    while True:
        network_time = finite(
            float(network.time), f"cantera-ignition-network-time:{case['id']}")
        if network_time >= maximum_time:
            break
        current_time = finite(
            float(network.step()), f"cantera-ignition-time:{case['id']}")
        current_temperature = finite(
            float(reactor.T), f"cantera-ignition-temperature:{case['id']}")
        if current_time <= previous_time:
            raise ComparisonError(f"cantera-ignition-time-not-increasing:{case['id']}")
        steps += 1
        if current_temperature >= threshold:
            temperature_span = current_temperature - previous_temperature
            if not math.isfinite(temperature_span) or temperature_span <= 0.0:
                raise ComparisonError(f"cantera-ignition-crossing-invalid:{case['id']}")
            fraction = finite(
                (threshold - previous_temperature) / temperature_span,
                f"cantera-ignition-fraction:{case['id']}",
            )
            if fraction < 0.0 or fraction > 1.0:
                raise ComparisonError(f"cantera-ignition-fraction-invalid:{case['id']}")
            result = previous_time + fraction * (current_time - previous_time)
            if not math.isfinite(result) or result <= 0.0:
                raise ComparisonError(f"cantera-ignition-invalid:{case['id']}")
            return result
        previous_time = current_time
        previous_temperature = current_temperature
        if steps > 2000000:
            break
    raise ComparisonError(f"cantera-ignition-not-detected:{case['id']}")


def rfam_derivatives_at_state(
    gas: ct.Solution,
    sample: dict[str, Any],
    species: list[str],
    references: list[float],
    heats: tuple[float, float],
) -> tuple[float, list[float], float]:
    gas.TDY = (
        sample["temperature_k"],
        sample["density_kg_m3"],
        sample["mass_fractions"],
    )
    density = finite(sample["density_kg_m3"], "rfam-derivative-density")
    if density <= 0.0:
        raise ComparisonError("rfam-derivative-density-nonpositive")
    production = gas.net_production_rates
    molecular_weights = gas.molecular_weights
    if len(production) != len(species) or len(molecular_weights) != len(species):
        raise ComparisonError("cantera-derivative-species-count-mismatch")
    mass_fraction_rates = [
        finite(
            finite(
                float(production[index]),
                f"cantera-net-production-rate:{species[index]}",
            )
            * finite(
                float(molecular_weights[index]),
                f"cantera-molecular-weight:{species[index]}",
            )
            / density,
            f"cantera-mass-fraction-rate:{species[index]}",
        )
        for index in range(len(species))
    ]
    qdot = finite(rfam_qdot(gas, heats), "cantera-derivative-qdot")
    sensible_species = [
        finite(
            (
                finite(
                    float(gas.species(name).thermo.h(sample["temperature_k"])),
                    f"cantera-species-enthalpy:{name}",
                )
                - finite(references[index], f"rfam-reference-enthalpy:{name}")
            )
            / finite(
                float(molecular_weights[index]),
                f"cantera-molecular-weight:{name}",
            ),
            f"cantera-sensible-species-enthalpy:{name}",
        )
        for index, name in enumerate(species)
    ]
    heat_capacity = finite(float(gas.cp_mass), "cantera-mixture-cp")
    if heat_capacity <= 0.0:
        raise ComparisonError("cantera-mixture-cp-nonpositive")
    temperature_rate = (
        qdot / density
        - sum(
            sensible * rate
            for sensible, rate in zip(sensible_species, mass_fraction_rates)
        )
    ) / heat_capacity
    temperature_rate = finite(
        temperature_rate, "cantera-temperature-derivative")
    return qdot, mass_fraction_rates, temperature_rate


def compare(root: Path, contract_path: Path, runner: Path) -> dict[str, Any]:
    validate_contract(root, contract_path)
    contract = load_json(contract_path, "contract")
    paths = artifact_paths(root, contract)
    species = list(contract["species_order"])
    mechanism = paths["cantera_same_mechanism"]
    instantaneous = load_json(paths["instantaneous_reference"], "instantaneous")
    heats = reaction_heats(instantaneous)

    probe_gas = ct.Solution(str(mechanism), "rfam-two-step")
    if ct.__version__ != contract["solver_profiles"]["cantera"]["version"]:
        raise ComparisonError("cantera-version-mismatch")
    if list(probe_gas.species_names) != species or probe_gas.n_reactions != 2:
        raise ComparisonError("cantera-mechanism-shape-mismatch")
    references, reference_audit = reference_enthalpies(
        probe_gas, species, contract["energy_reference"])
    shared_minimum = max(
        float(probe_gas.species(name).thermo.min_temp) for name in species)
    shared_maximum = min(
        float(probe_gas.species(name).thermo.max_temp) for name in species)
    envelope = contract["thermo_envelope"]
    if (
        shared_minimum != float(envelope["shared_minimum_temperature_k"])
        or shared_maximum != float(envelope["diagnostic_maximum_temperature_k"])
    ):
        raise ComparisonError("cantera-thermo-envelope-mismatch")

    sample_times = build_time_grid(contract)
    rfam_profile = contract["solver_profiles"]["rfam"]
    cantera_profile = contract["solver_profiles"]["cantera"]
    tolerance = contract["acceptance_tolerances"]
    comparison_maximum = float(envelope["comparison_maximum_temperature_k"])
    diagnostic_maximum = float(envelope["diagnostic_maximum_temperature_k"])

    counts = {"cases": 0, "total_samples": 0, "eligible_samples": 0,
              "guard_samples": 0, "outside_samples": 0}
    maximum_metrics = {
        "trajectory_temperature_absolute_k": 0.0,
        "trajectory_mass_fraction_absolute": 0.0,
        "ignition_delay_absolute_s": 0.0,
        "qdot_scaled_error": 0.0,
        "species_derivative_scaled_error": 0.0,
        "temperature_derivative_scaled_error": 0.0,
    }
    cases_evidence: list[dict[str, Any]] = []

    for case in contract["cases"]:
        headers, rfam_samples = run_rfam_case(
            runner, case, sample_times, species, rfam_profile)
        cantera_samples = cantera_trajectory(
            mechanism,
            case,
            sample_times,
            species,
            references,
            heats,
            cantera_profile,
        )
        cantera_delay = cantera_ignition_delay(
            mechanism,
            case,
            references,
            heats,
            cantera_profile,
            sample_times[-1],
        )
        rfam_delay = finite(headers["ignition_delay_s"], "rfam-ignition")
        ignition_error = abs(rfam_delay - cantera_delay)
        if ignition_error > float(tolerance["ignition_delay_absolute_s"]):
            raise ComparisonError(f"ignition-delay-mismatch:{case['id']}")
        maximum_metrics["ignition_delay_absolute_s"] = max(
            maximum_metrics["ignition_delay_absolute_s"], ignition_error)

        case_counts = {"eligible_samples": 0, "guard_samples": 0,
                       "outside_samples": 0}
        case_metrics = {
            "trajectory_temperature_absolute_k": 0.0,
            "trajectory_mass_fraction_absolute": 0.0,
            "qdot_scaled_error": 0.0,
            "species_derivative_scaled_error": 0.0,
            "temperature_derivative_scaled_error": 0.0,
        }
        audit_samples: list[dict[str, Any]] = []
        for index, (rfam, cantera) in enumerate(zip(rfam_samples, cantera_samples)):
            maximum_temperature = max(
                rfam["temperature_k"], cantera["temperature_k"])
            if maximum_temperature <= comparison_maximum:
                classification = "ELIGIBLE"
                case_counts["eligible_samples"] += 1
                temperature_error = abs(
                    rfam["temperature_k"] - cantera["temperature_k"])
                composition_error = max(
                    abs(actual - expected)
                    for actual, expected in zip(
                        rfam["mass_fractions"], cantera["mass_fractions"])
                )
                if temperature_error > float(
                    tolerance["trajectory_temperature_absolute_k"]
                ):
                    raise ComparisonError(
                        f"trajectory-temperature-mismatch:{case['id']}:{index}")
                if composition_error > float(
                    tolerance["trajectory_mass_fraction_absolute"]
                ):
                    raise ComparisonError(
                        f"trajectory-composition-mismatch:{case['id']}:{index}")

                qdot, species_rates, temperature_rate = rfam_derivatives_at_state(
                    probe_gas, rfam, species, references, heats)
                qdot_error = require_close(
                    rfam["qdot_w_m3"],
                    qdot,
                    float(tolerance["qdot_relative"]),
                    float(tolerance["qdot_absolute_w_per_m3"]),
                    f"{case['id']}[{index}].qdot",
                )
                species_error = 0.0
                for species_index, (actual, expected) in enumerate(
                    zip(rfam["mass_fraction_rates_s_inv"], species_rates)
                ):
                    species_error = max(
                        species_error,
                        require_close(
                            actual,
                            expected,
                            float(tolerance["species_derivative_relative"]),
                            float(tolerance["species_derivative_absolute_per_s"]),
                            f"{case['id']}[{index}].dYdt[{species_index}]",
                        ),
                    )
                temperature_rate_error = require_close(
                    rfam["temperature_rate_k_s"],
                    temperature_rate,
                    float(tolerance["temperature_derivative_relative"]),
                    float(tolerance["temperature_derivative_absolute_k_per_s"]),
                    f"{case['id']}[{index}].dTdt",
                )
                sample_errors: dict[str, Any] | None = {
                    "temperature_absolute_k": temperature_error,
                    "mass_fraction_max_absolute": composition_error,
                    "qdot_scaled": qdot_error,
                    "species_derivative_scaled": species_error,
                    "temperature_derivative_scaled": temperature_rate_error,
                }
                updates = {
                    "trajectory_temperature_absolute_k": temperature_error,
                    "trajectory_mass_fraction_absolute": composition_error,
                    "qdot_scaled_error": qdot_error,
                    "species_derivative_scaled_error": species_error,
                    "temperature_derivative_scaled_error": temperature_rate_error,
                }
                for key, value in updates.items():
                    case_metrics[key] = max(case_metrics[key], value)
                    maximum_metrics[key] = max(maximum_metrics[key], value)
            elif maximum_temperature <= diagnostic_maximum:
                classification = "THERMO_GUARD_DIAGNOSTIC_ONLY"
                case_counts["guard_samples"] += 1
                sample_errors = None
            else:
                classification = "OUTSIDE_THERMO_ENVELOPE"
                case_counts["outside_samples"] += 1
                sample_errors = None

            audit_samples.append(
                {
                    "index": index,
                    "time_s": rfam["time_s"],
                    "classification": classification,
                    "scored": classification == "ELIGIBLE",
                    "rfam_temperature_k": rfam["temperature_k"],
                    "cantera_temperature_k": cantera["temperature_k"],
                    "rfam_mass_fractions": rfam["mass_fractions"],
                    "cantera_mass_fractions": cantera["mass_fractions"],
                    "errors": sample_errors,
                }
            )

        counts["cases"] += 1
        counts["total_samples"] += len(audit_samples)
        for key in ("eligible_samples", "guard_samples", "outside_samples"):
            counts[key] += case_counts[key]
        cases_evidence.append(
            {
                "id": case["id"],
                "partition": case["partition"],
                "counts": case_counts,
                "metrics": case_metrics,
                "rfam_ignition_delay_s": rfam_delay,
                "cantera_ignition_delay_s": cantera_delay,
                "ignition_delay_absolute_error_s": ignition_error,
                "samples": audit_samples,
            }
        )

    if counts != contract["expected_counts"]:
        raise ComparisonError(
            "evidence-count-mismatch:actual=" + json.dumps(counts, sort_keys=True))

    return {
        "schema": contract["output_contract"]["schema"],
        "status": "implementation-equivalence-pass-physical-validation-blocked",
        "contract_sha256": sha256(contract_path),
        "runner_executable_sha256": sha256(runner),
        "cantera_version": ct.__version__,
        "species_order": species,
        "energy_reference": contract["energy_reference"],
        "energy_reference_construction": reference_audit,
        "solver_profiles": contract["solver_profiles"],
        "counts": counts,
        "maximum_metrics": maximum_metrics,
        "cases": cases_evidence,
        "approvals": {
            "fit_executed": False,
            "physical_validation_claimed": False,
            "physical_validation_approved": False,
            "production_ready": False,
            "reacting_hot_start_approved": False,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        root = args.root.resolve()
        contract = args.contract.resolve()
        runner = args.runner.resolve()
        output = args.output.resolve()
        if not root.is_dir():
            raise ComparisonError("root-not-directory")
        if not contract.is_file():
            raise ComparisonError("contract-missing")
        if not runner.is_file():
            raise ComparisonError("runner-missing")
        if output.exists():
            raise ComparisonError("output-already-exists")
        evidence = compare(root, contract, runner)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(
                evidence,
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=True,
                allow_nan=False,
            )
            + "\n",
            encoding="ascii",
        )
    except ComparisonError as error:
        print(f"R2_SAME_MECHANISM_TRAJECTORY_FAIL {error}", file=sys.stderr)
        return 2
    except Exception as error:
        print(
            "R2_SAME_MECHANISM_TRAJECTORY_FAIL unexpected:"
            + type(error).__name__,
            file=sys.stderr,
        )
        return 2

    counts = evidence["counts"]
    metrics = evidence["maximum_metrics"]
    print(
        "OUTSIDE_THERMO_ENVELOPE "
        f"samples={counts['outside_samples']} scored=0 "
        f"guard_samples={counts['guard_samples']}"
    )
    print(
        "R2_SAME_MECHANISM_TRAJECTORY_PASS "
        f"cases={counts['cases']} samples={counts['total_samples']} "
        f"eligible={counts['eligible_samples']} "
        f"max_temperature_error_k="
        f"{metrics['trajectory_temperature_absolute_k']:.8e} "
        f"max_mass_fraction_error="
        f"{metrics['trajectory_mass_fraction_absolute']:.8e}"
    )
    print(
        "fit_executed=0 physical_validation_approved=0 production_ready=0"
    )
    print("R2_PHYSICS_GATE_BLOCKED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
