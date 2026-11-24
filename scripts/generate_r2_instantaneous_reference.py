#!/usr/bin/env python3
"""Generate the pinned Cantera reference for the RFAM TwoStepCO rate gate.

This generator compares prescribed T/rho/Y states, not an evolving mixed-EOS
trajectory.  That isolates reaction orders, rate units, detailed-balance
direction, and the RFAM formation-energy qdot convention from the separate
real-fluid closure.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import cantera as ct


SPECIES = ["C3H8", "O2", "CO", "CO2", "H2O", "N2"]
RFAM_MW = {
    "C3H8": 0.04409562,
    "O2": 0.0319988,
    "CO": 0.0280101,
    "CO2": 0.0440095,
    "H2O": 0.01801528,
    "N2": 0.0280134,
}
RFAM_HF = {
    "C3H8": -103850.0,
    "O2": 0.0,
    "CO": -110530.0,
    "CO2": -393510.0,
    "H2O": -241820.0,
    "N2": 0.0,
}
RFAM_ELEMENTS = {
    "C3H8": (3, 8, 0, 0),
    "O2": (0, 0, 2, 0),
    "CO": (1, 0, 1, 0),
    "CO2": (1, 0, 2, 0),
    "H2O": (0, 2, 1, 0),
    "N2": (0, 0, 0, 2),
}
DH1 = -1.19502e6
DH2 = -2.8298e5

# Concentrations are mol/m3 and deliberately stay above the RFAM
# regularization transition so this gate tests the declared mass-action law.
CASES = [
    {"id": "r1_cold", "temperature_k": 900.0, "concentrations_mol_m3": {"C3H8": 2.0, "O2": 8.0}},
    {"id": "r1_hot", "temperature_k": 1400.0, "concentrations_mol_m3": {"C3H8": 20.0, "O2": 70.0}},
    {"id": "r1_order_base", "temperature_k": 1200.0, "concentrations_mol_m3": {"C3H8": 2.0, "O2": 8.0}},
    {"id": "r1_order_fuel_x2", "temperature_k": 1200.0, "concentrations_mol_m3": {"C3H8": 4.0, "O2": 8.0}},
    {"id": "r1_order_o2_x2", "temperature_k": 1200.0, "concentrations_mol_m3": {"C3H8": 2.0, "O2": 16.0}},
    {"id": "r2_forward", "temperature_k": 1200.0, "concentrations_mol_m3": {"CO": 10.0, "O2": 20.0, "CO2": 1.0}},
    {"id": "r2_order_base", "temperature_k": 1600.0, "concentrations_mol_m3": {"CO": 2.0, "O2": 4.0}},
    {"id": "r2_order_co_x2", "temperature_k": 1600.0, "concentrations_mol_m3": {"CO": 4.0, "O2": 4.0}},
    {"id": "r2_order_o2_x2", "temperature_k": 1600.0, "concentrations_mol_m3": {"CO": 2.0, "O2": 8.0}},
    {"id": "r2_reversible_hot", "temperature_k": 2500.0, "concentrations_mol_m3": {"CO": 5.0, "O2": 3.0, "CO2": 30.0}},
    {"id": "r2_dissociation_hot", "temperature_k": 3500.0, "concentrations_mol_m3": {"CO": 1.0, "O2": 0.5, "CO2": 40.0}},
    {"id": "mixed", "temperature_k": 1800.0, "concentrations_mol_m3": {"C3H8": 1.0, "O2": 4.0, "CO": 2.0, "CO2": 0.1, "H2O": 0.5}},
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def configure_from_concentrations(
    gas: ct.Solution,
    temperature: float,
    concentrations: dict[str, float],
) -> tuple[float, list[float], float, list[float]]:
    vector = [float(concentrations.get(name, 0.0)) for name in SPECIES]
    total = sum(vector)
    if total <= 0.0:
        raise ValueError("concentration vector is empty")
    density, mass_fractions = mass_fractions_from_concentrations(vector)
    # Both implementations receive this exact T/rho/Y state.  The custom
    # element weights in the Cantera mechanism make its species MW identical
    # to the RFAM registry, so concentration reconstruction is also identical.
    gas.TDY = temperature, density, mass_fractions
    actual = [value * 1000.0 for value in gas.concentrations]
    for expected, observed in zip(vector, actual):
        if abs(expected - observed) > 2.0e-12 * max(1.0, abs(expected)):
            raise ValueError(f"Cantera concentration mismatch {expected} != {observed}")
    return float(gas.P), actual, density, mass_fractions


def mass_fractions_from_concentrations(concentrations: list[float]) -> tuple[float, list[float]]:
    mass_density = [
        concentration * RFAM_MW[name]
        for name, concentration in zip(SPECIES, concentrations)
    ]
    density = sum(mass_density)
    if density <= 0.0:
        raise ValueError("mass density is not positive")
    return density, [value / density for value in mass_density]


def formation_qdot(net_production: list[float]) -> float:
    return -sum(
        rate * RFAM_HF[name]
        for name, rate in zip(SPECIES, net_production)
    )


def element_residual(net_production: list[float], element: int) -> float:
    return sum(
        rate * RFAM_ELEMENTS[name][element]
        for name, rate in zip(SPECIES, net_production)
    )


def generate(mechanism: Path, thermo: Path, species_registry: Path) -> dict:
    if mechanism.parent != thermo.parent:
        raise ValueError("mechanism and thermo artifacts must share a directory")
    thermo_import = f"{thermo.name}/species"
    if thermo_import not in mechanism.read_text(encoding="utf-8"):
        raise ValueError("mechanism does not reference the pinned thermo artifact")
    gas = ct.Solution(str(mechanism), "rfam-two-step")
    if list(gas.species_names) != SPECIES:
        raise ValueError("Cantera species order does not match the RFAM contract")
    for name, cantera_mw in zip(SPECIES, gas.molecular_weights):
        if abs(cantera_mw / 1000.0 - RFAM_MW[name]) > 2.0e-12 * RFAM_MW[name]:
            raise ValueError(f"molecular-weight mismatch for {name}")
    cases = []
    for item in CASES:
        temperature = float(item["temperature_k"])
        pressure, actual, density, mass_fractions = configure_from_concentrations(
            gas, temperature, item["concentrations_mol_m3"])
        forward = [float(value * 1000.0) for value in gas.forward_rates_of_progress]
        reverse = [float(value * 1000.0) for value in gas.reverse_rates_of_progress]
        net = [float(value * 1000.0) for value in gas.net_rates_of_progress]
        production = [float(value * 1000.0) for value in gas.net_production_rates]
        r1, r2 = net
        qdot = formation_qdot(production)
        cases.append(
            {
                "id": item["id"],
                "temperature_k": temperature,
                "pressure_pa": pressure,
                "density_kg_m3": density,
                "mass_fractions": mass_fractions,
                "concentrations_mol_m3": actual,
                "forward_rates_mol_m3_s": forward,
                "reverse_rates_mol_m3_s": reverse,
                "net_rates_mol_m3_s": net,
                "net_production_rates_mol_m3_s": production,
                "mass_source_kg_m3_s": [
                    rate * RFAM_MW[name]
                    for name, rate in zip(SPECIES, production)
                ],
                "qdot_w_m3": qdot,
                "element_residuals_mol_m3_s": {
                    element: element_residual(production, index)
                    for index, element in enumerate(("C", "H", "O", "N"))
                },
            }
        )
    return {
        "schema": "rfam.r2-instantaneous-reference/v1",
        "status": "diagnostic-equivalence-only",
        "cantera_version": ct.__version__,
        "mechanism_sha256": sha256(mechanism),
        "thermo_sha256": sha256(thermo),
        "species_registry_sha256": sha256(species_registry),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "thermo_import": thermo_import,
        "species_order": SPECIES,
        "rate_units": "mol/m3/s",
        "mass_source_units": "kg/m3/s",
        "energy_units": "W/m3",
        "energy_basis": "RFAM standard formation enthalpy at 298.15 K; sensible flow energy excluded",
        "comparison_tolerances": {
            "relative": 2.0e-8,
            "absolute_rate_mol_m3_s": 1.0e-8,
            "absolute_source_kg_m3_s": 1.0e-10,
            "absolute_qdot_w_m3": 1.0e-3,
            "closure_relative": 2.0e-12,
            "closure_absolute_mol_m3_s": 1.0e-9,
        },
        "reaction_contract": {
            "R1": {
                "stoich": {"C3H8": -1.0, "O2": -3.5, "CO": 3.0, "H2O": 4.0},
                "concentration_orders": {"C3H8": 0.1, "O2": 1.65},
                "A": 1.0e10,
                "Ea_J_per_mol": 125520.0,
                "DH_J_per_mol_progress": DH1,
            },
            "R2": {
                "stoich": {"CO": -1.0, "O2": -0.5, "CO2": 1.0},
                "concentration_orders": {"CO": 1.0, "O2": 0.5},
                "A_forward": 2.0e12,
                "Ea_forward_J_per_mol": 167360.0,
                "reverse": "detailed_balance_from_same_thermo",
                "DH_J_per_mol_progress": DH2,
            },
        },
        "cases": cases,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mechanism", type=Path, required=True)
    parser.add_argument("--thermo", type=Path, required=True)
    parser.add_argument("--species-registry", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    document = generate(
        args.mechanism.resolve(),
        args.thermo.resolve(),
        args.species_registry.resolve(),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(f"R2_REFERENCE_WRITTEN cases={len(document['cases'])}")
    print(f"mechanism_sha256={document['mechanism_sha256']}")
    print(f"thermo_sha256={document['thermo_sha256']}")
    print(f"species_registry_sha256={document['species_registry_sha256']}")
    print(f"generator_sha256={document['generator_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
