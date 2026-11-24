#!/usr/bin/env python3
"""Generate the frozen CoolProp 7.2.0 PR/TP-flash reference matrix."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import math
from pathlib import Path
from typing import Any

import CoolProp
from CoolProp import CoolProp as CP


EXPECTED_COOLPROP_VERSION = "7.2.0"
SCHEMA_VERSION = "rfam.pr_tp_flash_reference/v1"
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "tests"
    / "reference"
    / "pr_tp_flash_coolprop72.csv"
)

FIELDS = (
    "schema_version",
    "case_id",
    "record_kind",
    "backend",
    "coolprop_version",
    "components",
    "z_o2_mole",
    "z_c3h8_mole",
    "pressure_pa",
    "temperature_k",
    "phase",
    "quality_mole",
    "x_o2_mole",
    "x_c3h8_mole",
    "y_o2_mole",
    "y_c3h8_mole",
    "rho_mass_kg_m3",
    "h_mass_j_kg",
    "rho_liquid_mass_kg_m3",
    "rho_vapor_mass_kg_m3",
    "h_liquid_mass_j_kg",
    "h_vapor_mass_j_kg",
    "latent_heat_mass_j_kg",
    "notes",
)

FLOAT_FIELDS = frozenset(
    {
        "z_o2_mole",
        "z_c3h8_mole",
        "pressure_pa",
        "temperature_k",
        "quality_mole",
        "x_o2_mole",
        "x_c3h8_mole",
        "y_o2_mole",
        "y_c3h8_mole",
        "rho_mass_kg_m3",
        "h_mass_j_kg",
        "rho_liquid_mass_kg_m3",
        "rho_vapor_mass_kg_m3",
        "h_liquid_mass_j_kg",
        "h_vapor_mass_j_kg",
        "latent_heat_mass_j_kg",
    }
)

PHASE_NAMES = {
    CP.iphase_liquid: "liquid",
    CP.iphase_gas: "gas",
    CP.iphase_twophase: "two_phase",
    CP.iphase_supercritical: "supercritical",
    CP.iphase_supercritical_gas: "supercritical_gas",
    CP.iphase_supercritical_liquid: "supercritical_liquid",
    CP.iphase_critical_point: "critical_point",
    CP.iphase_unknown: "unknown",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def close(a: float, b: float, *, rtol: float = 1.0e-10, atol: float = 1.0e-10) -> bool:
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def base_row(case_id: str, record_kind: str, components: str) -> dict[str, Any]:
    row: dict[str, Any] = {field: None for field in FIELDS}
    row.update(
        {
            "schema_version": SCHEMA_VERSION,
            "case_id": case_id,
            "record_kind": record_kind,
            "backend": "HEOS",
            "coolprop_version": EXPECTED_COOLPROP_VERSION,
            "components": components,
        }
    )
    return row


def pure_o2_saturation(pressure_pa: float, case_id: str) -> dict[str, Any]:
    state = CoolProp.AbstractState("HEOS", "Oxygen")

    state.update(CP.PQ_INPUTS, pressure_pa, 0.0)
    temperature_liquid = state.T()
    rho_liquid = state.rhomass()
    h_liquid = state.hmass()

    state.update(CP.PQ_INPUTS, pressure_pa, 1.0)
    temperature_vapor = state.T()
    rho_vapor = state.rhomass()
    h_vapor = state.hmass()

    require(
        close(temperature_liquid, temperature_vapor, rtol=1.0e-12),
        f"{case_id}: liquid/vapor saturation temperatures disagree",
    )
    require(h_vapor > h_liquid, f"{case_id}: latent heat is not positive")

    row = base_row(case_id, "pure_saturation", "Oxygen")
    row.update(
        {
            "z_o2_mole": 1.0,
            "z_c3h8_mole": 0.0,
            "pressure_pa": pressure_pa,
            "temperature_k": temperature_liquid,
            "phase": "two_phase",
            "rho_liquid_mass_kg_m3": rho_liquid,
            "rho_vapor_mass_kg_m3": rho_vapor,
            "h_liquid_mass_j_kg": h_liquid,
            "h_vapor_mass_j_kg": h_vapor,
            "latent_heat_mass_j_kg": h_vapor - h_liquid,
            "notes": "fixed-pressure pure-fluid saturation",
        }
    )
    return row


def pure_o2_supercritical_sample(
    pressure_pa: float, temperature_k: float, case_id: str
) -> dict[str, Any]:
    state = CoolProp.AbstractState("HEOS", "Oxygen")
    require(
        pressure_pa > state.p_critical(),
        f"{case_id}: pressure must exceed the Oxygen critical pressure",
    )
    state.update(CP.PT_INPUTS, pressure_pa, temperature_k)
    phase = PHASE_NAMES.get(state.phase())
    require(phase is not None, f"{case_id}: unmapped CoolProp phase index")
    require(
        phase in {"supercritical", "supercritical_gas", "supercritical_liquid"},
        f"{case_id}: expected a supercritical phase label, got {phase}",
    )

    row = base_row(case_id, "supercritical_isobar_sample", "Oxygen")
    row.update(
        {
            "z_o2_mole": 1.0,
            "z_c3h8_mole": 0.0,
            "pressure_pa": pressure_pa,
            "temperature_k": temperature_k,
            "phase": phase,
            "rho_mass_kg_m3": state.rhomass(),
            "h_mass_j_kg": state.hmass(),
            "notes": "pressure_above_critical;not_a_saturation_record",
        }
    )
    return row


def binary_pair_parameters() -> dict[str, float]:
    oxygen_cas = CP.get_fluid_param_string("Oxygen", "CAS")
    propane_cas = CP.get_fluid_param_string("Propane", "CAS")

    # CoolProp stores this pair in Propane/Oxygen CAS order.
    parameters = {
        name: float(CP.get_mixture_binary_pair_data(propane_cas, oxygen_cas, name))
        for name in ("betaT", "gammaT", "betaV", "gammaV", "F")
    }
    expected = {"betaT": 1.0, "gammaT": 1.0, "betaV": 1.0, "gammaV": 1.0, "F": 0.0}
    require(
        parameters == expected,
        "CoolProp Oxygen/Propane binary parameters changed; review the frozen reference",
    )
    return parameters


def o2_propane_tp_flash() -> dict[str, Any]:
    case_id = "08_o2_c3h8_z050_050_p2mpa_t200k"
    pressure_pa = 2.0e6
    temperature_k = 200.0
    z = (0.5, 0.5)

    parameters = binary_pair_parameters()
    state = CoolProp.AbstractState("HEOS", "Oxygen&Propane")
    state.set_mole_fractions(z)
    state.update(CP.PT_INPUTS, pressure_pa, temperature_k)

    require(state.phase() == CP.iphase_twophase, f"{case_id}: expected a two-phase flash")
    quality = state.Q()
    x = tuple(state.mole_fractions_liquid())
    y = tuple(state.mole_fractions_vapor())

    require(0.0 < quality < 1.0, f"{case_id}: quality is outside (0,1)")
    require(len(x) == len(y) == len(z) == 2, f"{case_id}: unexpected composition size")
    require(close(sum(x), 1.0), f"{case_id}: liquid composition does not sum to one")
    require(close(sum(y), 1.0), f"{case_id}: vapor composition does not sum to one")
    for i, component in enumerate(("Oxygen", "Propane")):
        reconstructed = (1.0 - quality) * x[i] + quality * y[i]
        require(
            close(reconstructed, z[i], rtol=2.0e-10),
            f"{case_id}: molar lever rule fails for {component}",
        )

    parameter_note = "_".join(f"{name}={parameters[name]:g}" for name in parameters)
    row = base_row(case_id, "mixture_tp_flash", "Oxygen&Propane")
    row.update(
        {
            "z_o2_mole": z[0],
            "z_c3h8_mole": z[1],
            "pressure_pa": pressure_pa,
            "temperature_k": temperature_k,
            "phase": "two_phase",
            "quality_mole": quality,
            "x_o2_mole": x[0],
            "x_c3h8_mole": x[1],
            "y_o2_mole": y[0],
            "y_c3h8_mole": y[1],
            "rho_mass_kg_m3": state.rhomass(),
            "h_mass_j_kg": state.hmass(),
            "notes": (
                "Q_is_molar_vapor_fraction;"
                f"default_binary_pair_{parameter_note};"
                "no_fitted_departure_term;implementation_reference_only"
            ),
        }
    )
    return row


def build_rows() -> list[dict[str, Any]]:
    package_version = CoolProp.__version__
    library_version = CP.get_global_param_string("version")
    require(
        package_version == EXPECTED_COOLPROP_VERSION,
        f"expected CoolProp {EXPECTED_COOLPROP_VERSION}, got Python package {package_version}",
    )
    require(
        library_version == EXPECTED_COOLPROP_VERSION,
        f"expected CoolProp {EXPECTED_COOLPROP_VERSION}, got library {library_version}",
    )

    rows = [
        pure_o2_saturation(2.0e6, "01_o2_saturation_p2mpa"),
        pure_o2_saturation(4.0e6, "02_o2_saturation_p4mpa"),
    ]
    for index, temperature_k in enumerate((130.0, 150.0, 170.0, 200.0, 250.0), start=3):
        rows.append(
            pure_o2_supercritical_sample(
                6.0e6,
                temperature_k,
                f"{index:02d}_o2_supercritical_p6mpa_t{int(temperature_k):03d}k",
            )
        )
    rows.append(o2_propane_tp_flash())
    rows.sort(key=lambda row: row["case_id"])
    validate_rows(rows)
    return rows


def validate_rows(rows: list[dict[str, Any]]) -> None:
    require(len(rows) == 8, f"expected 8 reference rows, got {len(rows)}")
    case_ids = [row["case_id"] for row in rows]
    require(case_ids == sorted(case_ids), "rows are not in stable case_id order")
    require(len(case_ids) == len(set(case_ids)), "case_id values are not unique")

    for row in rows:
        require(tuple(row.keys()) == FIELDS, f"{row['case_id']}: CSV schema mismatch")
        require(row["schema_version"] == SCHEMA_VERSION, f"{row['case_id']}: schema mismatch")
        require(row["backend"] == "HEOS", f"{row['case_id']}: backend mismatch")
        require(
            row["coolprop_version"] == EXPECTED_COOLPROP_VERSION,
            f"{row['case_id']}: version mismatch",
        )
        for field in FLOAT_FIELDS:
            value = row[field]
            if value is not None:
                require(
                    isinstance(value, (int, float)) and math.isfinite(float(value)),
                    f"{row['case_id']}: {field} is not finite",
                )

        if row["record_kind"] == "pure_saturation":
            latent = row["h_vapor_mass_j_kg"] - row["h_liquid_mass_j_kg"]
            require(
                close(latent, row["latent_heat_mass_j_kg"], rtol=1.0e-13),
                f"{row['case_id']}: latent heat closure failed",
            )
        elif row["record_kind"] == "supercritical_isobar_sample":
            require(
                row["latent_heat_mass_j_kg"] is None,
                f"{row['case_id']}: supercritical sample must not contain latent heat",
            )
        elif row["record_kind"] == "mixture_tp_flash":
            require(row["phase"] == "two_phase", f"{row['case_id']}: mixture is not two-phase")
        else:
            raise RuntimeError(f"{row['case_id']}: unsupported record_kind")


def format_float(value: float) -> str:
    if value == 0.0:
        value = 0.0
    return f"{value:.17e}"


def render_csv(rows: list[dict[str, Any]]) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        rendered = {
            field: (
                ""
                if row[field] is None
                else format_float(float(row[field]))
                if field in FLOAT_FIELDS
                else str(row[field])
            )
            for field in FIELDS
        }
        writer.writerow(rendered)
    return output.getvalue().encode("ascii")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"output CSV (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify that the existing output is bitwise identical without rewriting it",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    payload = render_csv(build_rows())
    digest = hashlib.sha256(payload).hexdigest()
    output = args.output.resolve()

    if args.check:
        require(output.is_file(), f"reference file does not exist: {output}")
        require(output.read_bytes() == payload, f"reference file differs: {output}")
        print(f"REFERENCE_OK path={output} sha256={digest}")
        return 0

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(payload)
    print(f"REFERENCE_WRITTEN path={output} sha256={digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
