#!/usr/bin/env python3
"""Validate RFAM/Cantera same-mechanism instantaneous rate equivalence."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any


SCHEMA = "rfam.r2-instantaneous-reference/v1"
PROBE_SCHEMA = "rfam.r2-instantaneous-probe/v1"
SPECIES = ["C3H8", "O2", "CO", "CO2", "H2O", "N2"]
MW = [0.04409562, 0.0319988, 0.0280101, 0.0440095, 0.01801528, 0.0280134]
HF = [-103850.0, 0.0, -110530.0, -393510.0, -241820.0, 0.0]
ELEMENTS = [
    [3, 8, 0, 0],
    [0, 0, 2, 0],
    [1, 0, 1, 0],
    [1, 0, 2, 0],
    [0, 2, 1, 0],
    [0, 0, 0, 2],
]
CASE_IDS = [
    "r1_cold",
    "r1_hot",
    "r1_order_base",
    "r1_order_fuel_x2",
    "r1_order_o2_x2",
    "r2_forward",
    "r2_order_base",
    "r2_order_co_x2",
    "r2_order_o2_x2",
    "r2_reversible_hot",
    "r2_dissociation_hot",
    "mixed",
]
TOLERANCES = {
    "relative": 2.0e-8,
    "absolute_rate_mol_m3_s": 1.0e-8,
    "absolute_source_kg_m3_s": 1.0e-10,
    "absolute_qdot_w_m3": 1.0e-3,
    "closure_relative": 2.0e-12,
    "closure_absolute_mol_m3_s": 1.0e-9,
}
ROOT_KEYS = {
    "schema",
    "status",
    "cantera_version",
    "mechanism_sha256",
    "thermo_sha256",
    "species_registry_sha256",
    "generator_sha256",
    "thermo_import",
    "species_order",
    "rate_units",
    "mass_source_units",
    "energy_units",
    "energy_basis",
    "comparison_tolerances",
    "reaction_contract",
    "cases",
}
CASE_KEYS = {
    "id",
    "temperature_k",
    "pressure_pa",
    "density_kg_m3",
    "mass_fractions",
    "concentrations_mol_m3",
    "forward_rates_mol_m3_s",
    "reverse_rates_mol_m3_s",
    "net_rates_mol_m3_s",
    "net_production_rates_mol_m3_s",
    "mass_source_kg_m3_s",
    "qdot_w_m3",
    "element_residuals_mol_m3_s",
}


class GateError(ValueError):
    pass


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as exc:
        raise GateError(f"artifact-open-failed:{path}:{exc}") from exc


def finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise GateError(f"invalid-number:{label}")
    result = float(value)
    if not math.isfinite(result):
        raise GateError(f"nonfinite-number:{label}")
    return result


def finite_vector(value: Any, length: int, label: str) -> list[float]:
    if not isinstance(value, list) or len(value) != length:
        raise GateError(f"invalid-vector:{label}")
    return [finite_number(item, f"{label}[{index}]") for index, item in enumerate(value)]


def close(actual: float, expected: float, relative: float, absolute: float) -> bool:
    return abs(actual - expected) <= absolute + relative * abs(expected)


def require_close(
    actual: float,
    expected: float,
    relative: float,
    absolute: float,
    label: str,
) -> float:
    if not close(actual, expected, relative, absolute):
        raise GateError(
            f"value-mismatch:{label}:actual={actual:.17g}:expected={expected:.17g}"
        )
    scale = max(abs(expected), absolute)
    return abs(actual - expected) / scale


def validate_reference(
    reference: Path,
    generator: Path,
    mechanism: Path,
    thermo: Path,
    species_registry: Path,
) -> dict[str, Any]:
    try:
        document = json.loads(reference.read_text(encoding="utf-8"))
    except OSError as exc:
        raise GateError(f"reference-open-failed:{exc}") from exc
    except json.JSONDecodeError as exc:
        raise GateError(f"reference-json-invalid:{exc.msg}") from exc
    if not isinstance(document, dict):
        raise GateError("reference-not-object")
    if set(document) != ROOT_KEYS:
        raise GateError("reference-root-keys-mismatch")
    if document["schema"] != SCHEMA:
        raise GateError("reference-schema-mismatch")
    if document["status"] != "diagnostic-equivalence-only":
        raise GateError("reference-status-mismatch")
    if document["cantera_version"] != "3.1.0":
        raise GateError("reference-cantera-version-mismatch")
    if document["generator_sha256"] != sha256(generator):
        raise GateError("reference-generator-sha256-mismatch")
    if document["mechanism_sha256"] != sha256(mechanism):
        raise GateError("reference-mechanism-sha256-mismatch")
    if document["thermo_sha256"] != sha256(thermo):
        raise GateError("reference-thermo-sha256-mismatch")
    if document["species_registry_sha256"] != sha256(species_registry):
        raise GateError("reference-species-registry-sha256-mismatch")
    if document["thermo_import"] != f"{thermo.name}/species":
        raise GateError("reference-thermo-import-mismatch")
    if document["species_order"] != SPECIES:
        raise GateError("reference-species-order-mismatch")
    if document["rate_units"] != "mol/m3/s":
        raise GateError("reference-rate-units-mismatch")
    if document["mass_source_units"] != "kg/m3/s":
        raise GateError("reference-mass-source-units-mismatch")
    if document["energy_units"] != "W/m3":
        raise GateError("reference-energy-units-mismatch")
    if document["energy_basis"] != (
        "RFAM standard formation enthalpy at 298.15 K; sensible flow energy excluded"
    ):
        raise GateError("reference-energy-basis-mismatch")
    if document["comparison_tolerances"] != TOLERANCES:
        raise GateError("reference-tolerance-contract-mismatch")

    reaction_contract = document["reaction_contract"]
    if not isinstance(reaction_contract, dict) or set(reaction_contract) != {"R1", "R2"}:
        raise GateError("reference-reaction-contract-mismatch")
    if reaction_contract["R1"] != {
        "stoich": {"C3H8": -1.0, "O2": -3.5, "CO": 3.0, "H2O": 4.0},
        "concentration_orders": {"C3H8": 0.1, "O2": 1.65},
        "A": 1.0e10,
        "Ea_J_per_mol": 125520.0,
        "DH_J_per_mol_progress": -1.19502e6,
    }:
        raise GateError("reference-r1-contract-mismatch")
    if reaction_contract["R2"] != {
        "stoich": {"CO": -1.0, "O2": -0.5, "CO2": 1.0},
        "concentration_orders": {"CO": 1.0, "O2": 0.5},
        "A_forward": 2.0e12,
        "Ea_forward_J_per_mol": 167360.0,
        "reverse": "detailed_balance_from_same_thermo",
        "DH_J_per_mol_progress": -2.8298e5,
    }:
        raise GateError("reference-r2-contract-mismatch")

    cases = document["cases"]
    if not isinstance(cases, list) or [case.get("id") for case in cases] != CASE_IDS:
        raise GateError("reference-case-identity-mismatch")
    for case in cases:
        if not isinstance(case, dict) or set(case) != CASE_KEYS:
            raise GateError(f"reference-case-keys-mismatch:{case.get('id', 'unknown')}")
        case_id = case["id"]
        temperature = finite_number(case["temperature_k"], f"{case_id}.temperature")
        pressure = finite_number(case["pressure_pa"], f"{case_id}.pressure")
        density = finite_number(case["density_kg_m3"], f"{case_id}.density")
        if temperature <= 0.0 or pressure <= 0.0 or density <= 0.0:
            raise GateError(f"reference-state-not-positive:{case_id}")
        mass_fractions = finite_vector(case["mass_fractions"], 6, f"{case_id}.Y")
        concentrations = finite_vector(
            case["concentrations_mol_m3"], 6, f"{case_id}.C"
        )
        if any(value < 0.0 for value in mass_fractions + concentrations):
            raise GateError(f"reference-negative-composition:{case_id}")
        require_close(sum(mass_fractions), 1.0, 0.0, 2.0e-12, f"{case_id}.sumY")
        require_close(
            sum(concentration * mw for concentration, mw in zip(concentrations, MW)),
            density,
            2.0e-12,
            1.0e-14,
            f"{case_id}.density-from-concentrations",
        )
        for index, (concentration, mw, fraction) in enumerate(
            zip(concentrations, MW, mass_fractions)
        ):
            require_close(
                concentration * mw / density,
                fraction,
                2.0e-12,
                1.0e-14,
                f"{case_id}.Y[{index}]",
            )

        forward = finite_vector(case["forward_rates_mol_m3_s"], 2, f"{case_id}.forward")
        reverse = finite_vector(case["reverse_rates_mol_m3_s"], 2, f"{case_id}.reverse")
        net = finite_vector(case["net_rates_mol_m3_s"], 2, f"{case_id}.net")
        production = finite_vector(
            case["net_production_rates_mol_m3_s"], 6, f"{case_id}.production"
        )
        mass_source = finite_vector(
            case["mass_source_kg_m3_s"], 6, f"{case_id}.mass-source"
        )
        for reaction in range(2):
            require_close(
                forward[reaction] - reverse[reaction],
                net[reaction],
                2.0e-12,
                1.0e-10,
                f"{case_id}.net[{reaction}]",
            )
        expected_production = [
            -net[0],
            -3.5 * net[0] - 0.5 * net[1],
            3.0 * net[0] - net[1],
            net[1],
            4.0 * net[0],
            0.0,
        ]
        for index, expected in enumerate(expected_production):
            require_close(
                production[index], expected, 2.0e-12, 1.0e-8,
                f"{case_id}.stoich[{index}]",
            )
            require_close(
                mass_source[index], production[index] * MW[index],
                2.0e-12, 1.0e-10, f"{case_id}.mass-source[{index}]",
            )
        qdot = finite_number(case["qdot_w_m3"], f"{case_id}.qdot")
        require_close(
            qdot,
            -sum(rate * hf for rate, hf in zip(production, HF)),
            2.0e-12,
            1.0e-3,
            f"{case_id}.formation-qdot",
        )
        residuals = case["element_residuals_mol_m3_s"]
        if not isinstance(residuals, dict) or set(residuals) != {"C", "H", "O", "N"}:
            raise GateError(f"reference-element-residual-shape:{case_id}")
        for element_index, element in enumerate(("C", "H", "O", "N")):
            element_terms = [
                production[index] * ELEMENTS[index][element_index]
                for index in range(6)
            ]
            calculated = sum(element_terms)
            element_scale = sum(abs(term) for term in element_terms)
            stored = finite_number(residuals[element], f"{case_id}.element.{element}")
            # The stored residual is generated by a separate floating-point
            # reduction.  Compare it with a scale-aware tolerance so a
            # different summation order cannot reject an otherwise identical
            # source vector, while the subsequent element-zero check keeps the
            # physical closure threshold unchanged.
            require_close(
                stored,
                calculated,
                0.0,
                TOLERANCES["closure_absolute_mol_m3_s"]
                + TOLERANCES["closure_relative"] * element_scale,
                f"{case_id}.element.{element}",
            )
            require_close(
                stored,
                0.0,
                0.0,
                TOLERANCES["closure_absolute_mol_m3_s"]
                + TOLERANCES["closure_relative"] * element_scale,
                f"{case_id}.element-zero.{element}",
            )

    by_id = {case["id"]: case for case in cases}
    order_checks = [
        ("r1_order_fuel_x2", "r1_order_base", 0, 2.0 ** 0.1),
        ("r1_order_o2_x2", "r1_order_base", 0, 2.0 ** 1.65),
        ("r2_order_co_x2", "r2_order_base", 1, 2.0),
        ("r2_order_o2_x2", "r2_order_base", 1, 2.0 ** 0.5),
    ]
    for numerator_id, denominator_id, reaction, expected_ratio in order_checks:
        denominator = by_id[denominator_id]["forward_rates_mol_m3_s"][reaction]
        numerator = by_id[numerator_id]["forward_rates_mol_m3_s"][reaction]
        if denominator <= 0.0:
            raise GateError(f"reference-order-denominator-not-positive:{denominator_id}")
        require_close(
            numerator / denominator,
            expected_ratio,
            2.0e-12,
            1.0e-14,
            f"reaction-order:{numerator_id}",
        )
    return document


def parse_probe(stdout: str) -> dict[str, str]:
    values: dict[str, str] = {}
    marker = False
    for line in stdout.splitlines():
        if line == "R2_INSTANTANEOUS_PROBE_OK":
            marker = True
            continue
        if not line or "=" not in line:
            raise GateError(f"probe-output-invalid:{line}")
        key, value = line.split("=", 1)
        if not key or key in values:
            raise GateError(f"probe-output-key-invalid:{key}")
        values[key] = value
    if not marker:
        raise GateError("probe-success-marker-missing")
    if values.get("schema") != PROBE_SCHEMA:
        raise GateError("probe-schema-mismatch")
    return values


def probe_number(values: dict[str, str], key: str) -> float:
    if key not in values:
        raise GateError(f"probe-key-missing:{key}")
    try:
        result = float(values[key])
    except ValueError as exc:
        raise GateError(f"probe-number-invalid:{key}") from exc
    if not math.isfinite(result):
        raise GateError(f"probe-number-nonfinite:{key}")
    return result


def probe_vector(values: dict[str, str], key: str) -> list[float]:
    if key not in values:
        raise GateError(f"probe-key-missing:{key}")
    try:
        result = [float(item) for item in values[key].split(",")]
    except ValueError as exc:
        raise GateError(f"probe-vector-invalid:{key}") from exc
    if len(result) != 6 or any(not math.isfinite(item) for item in result):
        raise GateError(f"probe-vector-invalid:{key}")
    return result


def validate_probe_case(
    args: argparse.Namespace,
    case: dict[str, Any],
) -> tuple[int, float]:
    command = [
        str(args.probe),
        str(args.generator),
        str(args.reference),
        str(args.mechanism),
        str(args.thermo),
        str(args.rfam_mechanism),
        str(args.species_registry),
        str(args.mass_action_thermo),
        str(args.nasa7_thermo),
        str(args.mechanism_reader),
        str(args.thermo_provenance_header),
        str(args.thermo_provenance_source),
        format(case["temperature_k"], ".17g"),
        format(case["pressure_pa"], ".17g"),
        format(case["density_kg_m3"], ".17g"),
        *[format(value, ".17g") for value in case["mass_fractions"]],
    ]
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    if result.returncode != 0:
        raise GateError(
            f"probe-rejected:{case['id']}:rc={result.returncode}:"
            f"{result.stdout.strip()}:{result.stderr.strip()}"
        )
    values = parse_probe(result.stdout)
    relative = TOLERANCES["relative"]
    rate_absolute = TOLERANCES["absolute_rate_mol_m3_s"]
    source_absolute = TOLERANCES["absolute_source_kg_m3_s"]
    qdot_absolute = TOLERANCES["absolute_qdot_w_m3"]
    checks = 0
    errors: list[float] = []

    rate_pairs = [
        ("r1_forward_mol_m3_s", case["forward_rates_mol_m3_s"][0]),
        ("r2_forward_mol_m3_s", case["forward_rates_mol_m3_s"][1]),
        ("r2_reverse_mol_m3_s", case["reverse_rates_mol_m3_s"][1]),
        ("r2_net_mol_m3_s", case["net_rates_mol_m3_s"][1]),
    ]
    for key, expected in rate_pairs:
        errors.append(require_close(
            probe_number(values, key), expected, relative, rate_absolute,
            f"{case['id']}.{key}",
        ))
        checks += 1

    expected_source = case["mass_source_kg_m3_s"]
    two_step_source = probe_vector(values, "twostep_mass_source_kg_m3_s")
    dynamic_source = probe_vector(values, "dynamic_mass_source_kg_m3_s")
    for index in range(6):
        errors.append(require_close(
            two_step_source[index], expected_source[index], relative,
            source_absolute, f"{case['id']}.twostep-source[{index}]",
        ))
        errors.append(require_close(
            dynamic_source[index], expected_source[index], relative,
            source_absolute, f"{case['id']}.dynamic-source[{index}]",
        ))
        checks += 2
    for key in ("twostep_qdot_w_m3", "dynamic_qdot_w_m3"):
        errors.append(require_close(
            probe_number(values, key), case["qdot_w_m3"], relative,
            qdot_absolute, f"{case['id']}.{key}",
        ))
        checks += 1

    source_scale = sum(abs(value) for value in expected_source)
    formation_scale = abs(case["qdot_w_m3"]) + sum(
        abs(case["net_production_rates_mol_m3_s"][index] * HF[index])
        for index in range(6)
    )
    closure_checks = {
        "mass_residual_kg_m3_s": 1.0e-10
        + TOLERANCES["closure_relative"] * source_scale,
        "formation_energy_residual_w_m3": 1.0e-3
        + TOLERANCES["closure_relative"] * formation_scale,
    }
    for element_index, element in enumerate(("C", "H", "O", "N")):
        element_scale = sum(
            abs(case["net_production_rates_mol_m3_s"][index]
                * ELEMENTS[index][element_index])
            for index in range(6)
        )
        closure_checks[f"element_residual_{element}_mol_m3_s"] = (
            TOLERANCES["closure_absolute_mol_m3_s"]
            + TOLERANCES["closure_relative"] * element_scale
        )
    for key, absolute in closure_checks.items():
        require_close(probe_number(values, key), 0.0, 0.0, absolute, f"{case['id']}.{key}")
        checks += 1
    return checks, max(errors, default=0.0)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--generator", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--mechanism", type=Path, required=True)
    parser.add_argument("--thermo", type=Path, required=True)
    parser.add_argument("--rfam-mechanism", type=Path, required=True)
    parser.add_argument("--species-registry", type=Path, required=True)
    parser.add_argument("--mass-action-thermo", type=Path, required=True)
    parser.add_argument("--nasa7-thermo", type=Path, required=True)
    parser.add_argument("--mechanism-reader", type=Path, required=True)
    parser.add_argument("--thermo-provenance-header", type=Path, required=True)
    parser.add_argument("--thermo-provenance-source", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if not args.probe.is_file():
            raise GateError(f"probe-open-failed:{args.probe}")
        document = validate_reference(
            args.reference,
            args.generator,
            args.mechanism,
            args.thermo,
            args.species_registry,
        )
        total_checks = 0
        maximum_error = 0.0
        for case in document["cases"]:
            checks, error = validate_probe_case(args, case)
            total_checks += checks
            maximum_error = max(maximum_error, error)
    except (GateError, OSError) as error:
        print(f"R2_INSTANTANEOUS_EQUIVALENCE_FAIL {error}", file=sys.stderr)
        return 2

    print(
        "R2_INSTANTANEOUS_EQUIVALENCE_PASS "
        f"cases={len(document['cases'])} checks={total_checks} "
        f"order_checks=4 max_scaled_error={maximum_error:.8e}"
    )
    print("physical_holdout_validated=0 production_ready=0")
    print("R2_PHYSICS_GATE_BLOCKED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
