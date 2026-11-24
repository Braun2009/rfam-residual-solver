#!/usr/bin/env python3
"""Screen the thermodynamic limit of the six-species R2 mechanism."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import cantera as ct


EXPECTED_CANTERA_VERSION = "3.1.0"
MAJOR_SPECIES = ("C3H8", "O2", "CO", "CO2", "H2O", "N2")
REQUIRED_DISSOCIATION_SPECIES = ("H2", "OH", "O", "H")


class ProbeError(RuntimeError):
    pass


def strict_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ProbeError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def load_json(path: Path) -> Any:
    try:
        return json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=strict_object,
        )
    except ProbeError:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ProbeError(f"json-load-failed:{path}:{error}") from error


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ProbeError(message)


def finite(value: Any, location: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ProbeError(f"invalid-number:{location}")
    result = float(value)
    if not math.isfinite(result):
        raise ProbeError(f"nonfinite-number:{location}")
    return result


def equilibrium_state(
    phase: ct.Solution,
    temperature: float,
    pressure: float,
    mass_fractions: dict[str, float],
) -> dict[str, float]:
    phase.TPY = temperature, pressure, mass_fractions
    phase.equilibrate("HP", solver="auto", rtol=1.0e-10, max_steps=20000)
    return {
        "temperature_k": float(phase.T),
        **{
            f"{species}_mass_fraction": float(phase[species].Y[0])
            for species in MAJOR_SPECIES
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    args = parser.parse_args()

    try:
        root = args.root.resolve()
        detailed_path = root / "tests/reference/cantera_gri30_v3_0.yaml"
        reduced_path = root / "tests/reference/rfam_twostep_cantera_v1.yaml"
        trajectory_path = root / "tests/reference/r2_0d_cantera_reference_v1.json"
        instantaneous_path = root / "tests/reference/r2_instantaneous_reference_v1.json"
        for path in (
            detailed_path,
            reduced_path,
            trajectory_path,
            instantaneous_path,
        ):
            require(path.is_file(), f"missing-artifact:{path}")

        require(
            ct.__version__ == EXPECTED_CANTERA_VERSION,
            f"cantera-version-mismatch:{ct.__version__}",
        )
        trajectory = load_json(trajectory_path)
        instantaneous = load_json(instantaneous_path)
        require(
            trajectory.get("mechanism_sha256") == sha256(detailed_path),
            "detailed-mechanism-hash-mismatch",
        )
        require(
            instantaneous.get("mechanism_sha256") == sha256(reduced_path),
            "reduced-mechanism-hash-mismatch",
        )
        require(
            trajectory.get("major_species") == list(MAJOR_SPECIES),
            "major-species-contract-mismatch",
        )
        cases = trajectory.get("cases")
        require(isinstance(cases, list) and len(cases) == 6, "case-matrix-mismatch")

        detailed = ct.Solution(str(detailed_path))
        reduced = ct.Solution(str(reduced_path), "rfam-two-step")
        require(
            tuple(reduced.species_names) == MAJOR_SPECIES,
            "reduced-species-set-mismatch",
        )
        require(
            all(name in detailed.species_names for name in REQUIRED_DISSOCIATION_SPECIES),
            "detailed-mechanism-misses-dissociation-species",
        )

        fieldnames = [
            "case_id",
            "partition",
            "detailed_equilibrium_temperature_k",
            "reduced_equilibrium_temperature_k",
            "temperature_gap_k",
            "reference_terminal_temperature_residual_k",
            "maximum_major_species_mass_fraction_gap",
            "fit_approved",
            "physical_validation_approved",
            "production_ready",
        ]
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, lineterminator="\n")
        print("# rfam.r2-six-species-equilibrium-screen/v1")
        writer.writeheader()

        temperature_gaps: list[float] = []
        reference_residuals: list[float] = []
        species_gaps: list[float] = []
        for index, case in enumerate(cases):
            require(isinstance(case, dict), f"invalid-case:{index}")
            case_id = case.get("id")
            partition = case.get("partition")
            require(isinstance(case_id, str) and case_id, f"invalid-case-id:{index}")
            require(
                partition in ("protocol-training", "protocol-holdout"),
                f"invalid-partition:{case_id}",
            )
            initial = case.get("initial_major_species_mass_fractions")
            require(
                isinstance(initial, dict) and set(initial) == set(MAJOR_SPECIES),
                f"invalid-initial-composition:{case_id}",
            )
            mass_fractions = {
                name: finite(initial[name], f"{case_id}.{name}")
                for name in MAJOR_SPECIES
            }
            require(
                abs(sum(mass_fractions.values()) - 1.0) <= 2.0e-10,
                f"initial-composition-sum-mismatch:{case_id}",
            )
            temperature = finite(case.get("temperature_k"), f"{case_id}.temperature")
            pressure = finite(case.get("pressure_pa"), f"{case_id}.pressure")
            detailed_state = equilibrium_state(
                detailed, temperature, pressure, mass_fractions
            )
            reduced_state = equilibrium_state(
                reduced, temperature, pressure, mass_fractions
            )
            samples = case.get("samples")
            require(isinstance(samples, list) and samples, f"missing-samples:{case_id}")
            terminal = samples[-1]
            require(isinstance(terminal, dict), f"invalid-terminal-sample:{case_id}")
            reference_temperature = finite(
                terminal.get("temperature_k"), f"{case_id}.terminal-temperature"
            )
            reference_composition = terminal.get("major_species_mass_fractions")
            require(
                isinstance(reference_composition, dict)
                and set(reference_composition) == set(MAJOR_SPECIES),
                f"invalid-terminal-composition:{case_id}",
            )

            temperature_gap = (
                reduced_state["temperature_k"] - detailed_state["temperature_k"]
            )
            reference_residual = abs(
                detailed_state["temperature_k"] - reference_temperature
            )
            maximum_species_gap = max(
                abs(
                    reduced_state[f"{name}_mass_fraction"]
                    - finite(
                        reference_composition[name],
                        f"{case_id}.terminal.{name}",
                    )
                )
                for name in MAJOR_SPECIES
            )
            require(
                temperature_gap > 0.0,
                f"six-species-temperature-gap-not-positive:{case_id}",
            )
            require(
                reference_residual <= 1.0e-3,
                f"reference-terminal-not-equilibrated:{case_id}",
            )
            temperature_gaps.append(temperature_gap)
            reference_residuals.append(reference_residual)
            species_gaps.append(maximum_species_gap)
            writer.writerow(
                {
                    "case_id": case_id,
                    "partition": partition,
                    "detailed_equilibrium_temperature_k": format(
                        detailed_state["temperature_k"], ".17g"
                    ),
                    "reduced_equilibrium_temperature_k": format(
                        reduced_state["temperature_k"], ".17g"
                    ),
                    "temperature_gap_k": format(temperature_gap, ".17g"),
                    "reference_terminal_temperature_residual_k": format(
                        reference_residual, ".17g"
                    ),
                    "maximum_major_species_mass_fraction_gap": format(
                        maximum_species_gap, ".17g"
                    ),
                    "fit_approved": 0,
                    "physical_validation_approved": 0,
                    "production_ready": 0,
                }
            )

        print(
            "# R2_SIX_SPECIES_EQUILIBRIUM_BLOCKED"
            f" cases={len(cases)}"
            f" min_temperature_gap_k={min(temperature_gaps):.17g}"
            f" max_temperature_gap_k={max(temperature_gaps):.17g}"
            f" max_reference_terminal_residual_k={max(reference_residuals):.17g}"
            f" max_major_species_gap={max(species_gaps):.17g}"
            " fit_approved=0 physical_validation_approved=0 production_ready=0"
        )
        return 0
    except (ProbeError, OSError, RuntimeError, ValueError, ct.CanteraError) as error:
        print(f"R2_SIX_SPECIES_EQUILIBRIUM_FAIL reason={error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
