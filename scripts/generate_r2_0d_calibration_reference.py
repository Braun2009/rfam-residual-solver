#!/usr/bin/env python3
"""Generate a local Cantera 0D protocol reference for the R2 pre-fit gate.

The output is deliberately a *local reference trajectory*, not an approved
kinetics fit target.  Its only inputs are the checked-in GRI30 YAML and the
predeclared high-temperature, gas-only condition matrix below.  The companion
contract keeps mechanism selection, experimental-data selection, fitting, and
physical validation blocked.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import cantera as ct


CANONICAL_CANTERA_VERSION = "3.1.0"
CANONICAL_PYTHON_EXECUTABLE = "/usr/bin/python3.8"
CANONICAL_PYTHON_SHA256 = "298a9e830ed52f36c299427565485d717d1ce0179c0597cc16560513eb780b06"
CANONICAL_CANTERA_MODULE = "/home/nvoyager/.local/lib/python3.8/site-packages/cantera/__init__.py"
CANONICAL_CANTERA_MODULE_SHA256 = "9a7ceae3bae1a7e41dc581e353d44e04af8447114748f5b7c75995e26f8819b9"
CANONICAL_CANTERA_EXTENSION = "/home/nvoyager/.local/lib/python3.8/site-packages/cantera/_cantera.cpython-38-x86_64-linux-gnu.so"
CANONICAL_CANTERA_EXTENSION_SHA256 = "c57a47c201b0cb8f965b38acc7ebfafbbe95db7464d80359499c64d50555fae8"
MAJOR_SPECIES = ("C3H8", "O2", "CO", "CO2", "H2O", "N2")
SAMPLE_TIMES_S = (
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
)

# These are high-temperature, ideal-gas, constant-pressure protocol points.
# They intentionally exclude the cryogenic/VLE domain that belongs to RFAM
# rho-e-Y HEM, not to a Cantera ideal-gas reference reactor.
CASES = (
    {
        "id": "train_phi0p8_diluted",
        "partition": "protocol-training",
        "temperature_k": 1200.0,
        "pressure_pa": 4.0e6,
        "equivalence_ratio": 0.8,
        "fuel": "C3H8:1.0",
        "oxidizer": "O2:1.0,N2:0.25",
    },
    {
        "id": "train_phi1p0",
        "partition": "protocol-training",
        "temperature_k": 1300.0,
        "pressure_pa": 4.0e6,
        "equivalence_ratio": 1.0,
        "fuel": "C3H8:1.0",
        "oxidizer": "O2:1.0",
    },
    {
        "id": "train_phi1p2",
        "partition": "protocol-training",
        "temperature_k": 1500.0,
        "pressure_pa": 2.0e6,
        "equivalence_ratio": 1.2,
        "fuel": "C3H8:1.0",
        "oxidizer": "O2:1.0",
    },
    {
        "id": "holdout_phi0p6_diluted",
        "partition": "protocol-holdout",
        "temperature_k": 1250.0,
        "pressure_pa": 2.0e6,
        "equivalence_ratio": 0.6,
        "fuel": "C3H8:1.0",
        "oxidizer": "O2:1.0,N2:0.5",
    },
    {
        "id": "holdout_phi1p1",
        "partition": "protocol-holdout",
        "temperature_k": 1400.0,
        "pressure_pa": 6.0e6,
        "equivalence_ratio": 1.1,
        "fuel": "C3H8:1.0",
        "oxidizer": "O2:1.0",
    },
    {
        "id": "holdout_phi1p4",
        "partition": "protocol-holdout",
        "temperature_k": 1600.0,
        "pressure_pa": 1.0e6,
        "equivalence_ratio": 1.4,
        "fuel": "C3H8:1.0",
        "oxidizer": "O2:1.0",
    },
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_float(value: float) -> float:
    """Keep the frozen JSON portable across same-version local reruns."""
    if not math.isfinite(value):
        raise ValueError("Cantera produced a non-finite trajectory value")
    return float(format(float(value), ".12g"))


def runtime_identity() -> dict[str, str]:
    python_executable = Path(sys.executable).resolve()
    module = Path(ct.__file__).resolve()
    extension = Path(ct._cantera.__file__).resolve()
    identity = {
        "python_executable": str(python_executable),
        "python_sha256": sha256(python_executable),
        "cantera_module": str(module),
        "cantera_module_sha256": sha256(module),
        "cantera_extension": str(extension),
        "cantera_extension_sha256": sha256(extension),
    }
    expected = {
        "python_executable": CANONICAL_PYTHON_EXECUTABLE,
        "python_sha256": CANONICAL_PYTHON_SHA256,
        "cantera_module": CANONICAL_CANTERA_MODULE,
        "cantera_module_sha256": CANONICAL_CANTERA_MODULE_SHA256,
        "cantera_extension": CANONICAL_CANTERA_EXTENSION,
        "cantera_extension_sha256": CANONICAL_CANTERA_EXTENSION_SHA256,
    }
    if identity != expected:
        raise ValueError("local Cantera runtime identity does not match the frozen protocol")
    return identity


def major_mass_fractions(phase: ct.Solution) -> dict[str, float]:
    return {
        name: stable_float(float(phase[name].Y[0]))
        for name in MAJOR_SPECIES
    }


def configure_case(mechanism: Path, definition: dict[str, Any]) -> ct.Solution:
    gas = ct.Solution(str(mechanism))
    if gas.thermo_model != "ideal-gas":
        raise ValueError("local Cantera mechanism is not an ideal-gas phase")
    missing = [name for name in MAJOR_SPECIES if name not in gas.species_names]
    if missing:
        raise ValueError("local Cantera mechanism misses required species: " + ",".join(missing))
    gas.TP = definition["temperature_k"], definition["pressure_pa"]
    gas.set_equivalence_ratio(
        definition["equivalence_ratio"], definition["fuel"], definition["oxidizer"])
    return gas


def record_state(time_s: float, reactor: ct.IdealGasConstPressureReactor) -> dict[str, Any]:
    phase = reactor.thermo
    return {
        "time_s": stable_float(time_s),
        "temperature_k": stable_float(float(phase.T)),
        "pressure_pa": stable_float(float(phase.P)),
        "density_kg_m3": stable_float(float(phase.density)),
        "heat_release_rate_w_m3": stable_float(float(phase.heat_release_rate)),
        "major_species_mass_fractions": major_mass_fractions(phase),
    }


def ignition_delay_s(mechanism: Path, definition: dict[str, Any]) -> float:
    gas = configure_case(mechanism, definition)
    reactor = ct.IdealGasConstPressureReactor(gas, energy="on")
    network = ct.ReactorNet([reactor])
    network.rtol = 1.0e-10
    network.atol = 1.0e-18
    initial_temperature = float(reactor.T)
    while network.time < SAMPLE_TIMES_S[-1]:
        time_s = float(network.step())
        if reactor.T - initial_temperature >= 400.0:
            return stable_float(time_s)
    raise ValueError("protocol case did not reach the 400 K ignition threshold")


def run_case(mechanism: Path, definition: dict[str, Any]) -> dict[str, Any]:
    gas = configure_case(mechanism, definition)
    reactor = ct.IdealGasConstPressureReactor(gas, energy="on")
    network = ct.ReactorNet([reactor])
    network.rtol = 1.0e-10
    network.atol = 1.0e-18
    samples = [record_state(0.0, reactor)]
    for time_s in SAMPLE_TIMES_S[1:]:
        network.advance(time_s)
        samples.append(record_state(time_s, reactor))
    return {
        **definition,
        "ignition_delay_definition": "first adaptive integrator time where T - T0 >= 400 K",
        "ignition_delay_s": ignition_delay_s(mechanism, definition),
        "initial_major_species_mass_fractions": samples[0]["major_species_mass_fractions"],
        "samples": samples,
    }


def generate(mechanism: Path) -> dict[str, Any]:
    if ct.__version__ != CANONICAL_CANTERA_VERSION:
        raise ValueError(
            "Cantera version mismatch: expected "
            f"{CANONICAL_CANTERA_VERSION}, got {ct.__version__}")
    if not mechanism.is_file():
        raise ValueError(f"Cantera mechanism does not exist: {mechanism}")
    identity = runtime_identity()
    return {
        "schema": "rfam.r2-0d-cantera-reference/v1",
        "status": "local-cantera-reference-not-fit-eligible",
        "cantera_version": ct.__version__,
        "runtime_identity": identity,
        "generator_sha256": sha256(Path(__file__).resolve()),
        "mechanism_sha256": sha256(mechanism),
        "reactor": {
            "model": "IdealGasConstPressureReactor",
            "phase": "ideal-gas",
            "energy": "on",
            "rtol": 1.0e-10,
            "atol": 1.0e-18,
            "phase_scope": "high-temperature gas-only; VLE, near-critical and solid-risk states excluded",
        },
        "sample_times_s": [stable_float(value) for value in SAMPLE_TIMES_S],
        "major_species": list(MAJOR_SPECIES),
        "units": {
            "temperature": "K",
            "pressure": "Pa",
            "density": "kg/m3",
            "heat_release_rate": "W/m3",
            "mass_fraction": "1",
            "time": "s",
        },
        "cases": [run_case(mechanism, definition) for definition in CASES],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mechanism", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    document = generate(args.mechanism.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(document, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n",
        encoding="utf-8",
    )
    print("R2_0D_REFERENCE_WRITTEN cases=" + str(len(document["cases"])))
    print("generator_sha256=" + document["generator_sha256"])
    print("mechanism_sha256=" + document["mechanism_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
