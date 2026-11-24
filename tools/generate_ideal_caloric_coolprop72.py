#!/usr/bin/env python3
"""Generate frozen low-temperature ideal-caloric data from CoolProp 7.2.0."""

from __future__ import annotations

import argparse
import csv
import io
import math
from pathlib import Path
from typing import Any, Iterable

import CoolProp
from CoolProp import CoolProp as CP


EXPECTED_COOLPROP_VERSION = "7.2.0"
SCHEMA_VERSION = "rfam.ideal_caloric_coolprop72/v1"
REFERENCE_TEMPERATURE_K = 298.15
QUERY_PRESSURE_PA = 100000.0
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "tests"
    / "reference"
    / "ideal_caloric_coolprop72.csv"
)

SPECIES = (
    ("O2", "Oxygen", ((80.0, 300.0, 5.0),)),
    (
        "C3H8",
        "Propane",
        ((80.0, 150.0, 2.5), (150.0, 350.0, 5.0)),
    ),
)

FIELDS = (
    "schema_version",
    "species",
    "coolprop_fluid",
    "backend",
    "coolprop_version",
    "source_api",
    "enthalpy_property",
    "cp_property",
    "reference_temperature_k",
    "query_pressure_pa",
    "temperature_k",
    "hsens_ideal_j_kg",
    "cp0_j_kg_k",
    "notes",
)

FLOAT_FIELDS = frozenset(
    {
        "reference_temperature_k",
        "query_pressure_pa",
        "temperature_k",
        "hsens_ideal_j_kg",
        "cp0_j_kg_k",
    }
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative_error(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1.0)


def temperatures(
    segments: tuple[tuple[float, float, float], ...]
) -> list[float]:
    values: list[float] = []
    for lower, upper, spacing in segments:
        require(
            lower > 0.0 and upper > lower and spacing > 0.0,
            "invalid temperature segment",
        )
        count = int(round((upper - lower) / spacing))
        require(
            abs(lower + spacing * count - upper) < 1.0e-12,
            "temperature segment is not exactly divisible by its spacing",
        )
        values.extend(lower + spacing * i for i in range(count + 1))
    lower = min(segment[0] for segment in segments)
    upper = max(segment[1] for segment in segments)
    if lower < REFERENCE_TEMPERATURE_K < upper:
        values.append(REFERENCE_TEMPERATURE_K)
    return sorted(set(values))


def query(fluid: str, temperature: float, pressure: float) -> tuple[float, float]:
    h = CP.PropsSI(
        "HMASS_IDEALGAS", "T", temperature, "P", pressure, fluid
    )
    cp = CP.PropsSI("CP0MASS", "T", temperature, "P", pressure, fluid)
    require(math.isfinite(h), f"{fluid} {temperature:g} K: non-finite ideal enthalpy")
    require(
        math.isfinite(cp) and cp > 0.0,
        f"{fluid} {temperature:g} K: non-positive CP0MASS",
    )
    return h, cp


def hermite(
    t0: float,
    h0: float,
    cp0: float,
    t1: float,
    h1: float,
    cp1: float,
    temperature: float,
) -> tuple[float, float]:
    dt = t1 - t0
    u = (temperature - t0) / dt
    u2 = u * u
    u3 = u2 * u
    h00 = 2.0 * u3 - 3.0 * u2 + 1.0
    h10 = u3 - 2.0 * u2 + u
    h01 = -2.0 * u3 + 3.0 * u2
    h11 = u3 - u2
    h = h00 * h0 + h10 * dt * cp0 + h01 * h1 + h11 * dt * cp1

    dh00 = 6.0 * u2 - 6.0 * u
    dh10 = 3.0 * u2 - 4.0 * u + 1.0
    dh01 = -6.0 * u2 + 6.0 * u
    dh11 = 3.0 * u2 - 2.0 * u
    cp = (
        (dh00 * h0 + dh01 * h1) / dt
        + dh10 * cp0
        + dh11 * cp1
    )
    return h, cp


def audit_source_derivative(fluid: str, temperature: float, cp: float) -> None:
    delta = 1.0e-3
    h_upper, _ = query(fluid, temperature + delta, QUERY_PRESSURE_PA)
    h_lower, _ = query(fluid, temperature - delta, QUERY_PRESSURE_PA)
    derivative = (h_upper - h_lower) / (2.0 * delta)
    require(
        relative_error(derivative, cp) < 2.0e-8,
        f"{fluid} {temperature:g} K: d(HMASS_IDEALGAS)/dT != CP0MASS",
    )


def build_species_rows(
    species: str,
    fluid: str,
    segments: tuple[tuple[float, float, float], ...],
) -> list[dict[str, Any]]:
    reference_h, _ = query(
        fluid, REFERENCE_TEMPERATURE_K, QUERY_PRESSURE_PA
    )
    rows: list[dict[str, Any]] = []
    for temperature in temperatures(segments):
        h, cp = query(fluid, temperature, QUERY_PRESSURE_PA)
        h_second_pressure, cp_second_pressure = query(
            fluid, temperature, 2.0 * QUERY_PRESSURE_PA
        )
        require(
            relative_error(h, h_second_pressure) < 2.0e-13,
            f"{fluid} {temperature:g} K: ideal enthalpy depends on query pressure",
        )
        require(
            relative_error(cp, cp_second_pressure) < 2.0e-13,
            f"{fluid} {temperature:g} K: CP0MASS depends on query pressure",
        )
        audit_source_derivative(fluid, temperature, cp)
        hsens = h - reference_h
        if temperature == REFERENCE_TEMPERATURE_K:
            require(hsens == 0.0, f"{fluid}: reference subtraction is not exact")
        rows.append(
            {
                "schema_version": SCHEMA_VERSION,
                "species": species,
                "coolprop_fluid": fluid,
                "backend": "HEOS",
                "coolprop_version": EXPECTED_COOLPROP_VERSION,
                "source_api": "PropsSI",
                "enthalpy_property": "HMASS_IDEALGAS",
                "cp_property": "CP0MASS",
                "reference_temperature_k": REFERENCE_TEMPERATURE_K,
                "query_pressure_pa": QUERY_PRESSURE_PA,
                "temperature_k": temperature,
                "hsens_ideal_j_kg": hsens,
                "cp0_j_kg_k": cp,
                "notes": (
                    "ideal_caloric_only;not_liquid_or_two_phase_transport;"
                    "not_HEOS_residual"
                ),
            }
        )

    for left, right in zip(rows, rows[1:]):
        t0 = float(left["temperature_k"])
        t1 = float(right["temperature_k"])
        h0 = float(left["hsens_ideal_j_kg"])
        h1 = float(right["hsens_ideal_j_kg"])
        cp0 = float(left["cp0_j_kg_k"])
        cp1 = float(right["cp0_j_kg_k"])
        require(h1 > h0, f"{fluid}: frozen ideal enthalpy is not monotone")
        for fraction in (0.125, 0.25, 0.5, 0.75, 0.875):
            temperature = t0 + fraction * (t1 - t0)
            interpolated_h, interpolated_cp = hermite(
                t0, h0, cp0, t1, h1, cp1, temperature
            )
            source_h, source_cp = query(fluid, temperature, QUERY_PRESSURE_PA)
            source_h -= reference_h
            require(
                interpolated_cp > 0.0,
                f"{fluid} {temperature:g} K: Hermite cp is non-positive",
            )
            require(
                relative_error(interpolated_h, source_h) < 1.0e-7,
                f"{fluid} {temperature:g} K: Hermite enthalpy audit failed",
            )
            require(
                relative_error(interpolated_cp, source_cp) < 1.0e-6,
                f"{fluid} {temperature:g} K: Hermite cp audit failed",
            )
    return rows


def build_rows() -> list[dict[str, Any]]:
    require(
        CoolProp.__version__ == EXPECTED_COOLPROP_VERSION,
        "expected CoolProp "
        f"{EXPECTED_COOLPROP_VERSION}, got Python package {CoolProp.__version__}",
    )
    library_version = CP.get_global_param_string("version")
    require(
        library_version == EXPECTED_COOLPROP_VERSION,
        "expected CoolProp "
        f"{EXPECTED_COOLPROP_VERSION}, got library {library_version}",
    )
    rows: list[dict[str, Any]] = []
    for species, fluid, segments in SPECIES:
        rows.extend(build_species_rows(species, fluid, segments))
    return rows


def format_float(value: float) -> str:
    require(math.isfinite(value), "attempted to serialize a non-finite value")
    return format(value, ".17e")


def render_csv(rows: Iterable[dict[str, Any]]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        serialized = dict(row)
        for field in FLOAT_FIELDS:
            serialized[field] = format_float(float(serialized[field]))
        writer.writerow(serialized)
    return stream.getvalue()


def render_cpp(rows: Iterable[dict[str, Any]]) -> str:
    grouped: dict[str, list[dict[str, Any]]] = {name: [] for name, *_ in SPECIES}
    for row in rows:
        grouped[str(row["species"])].append(row)
    output: list[str] = []
    for species, *_ in SPECIES:
        identifier = species.lower().replace("3h8", "3h8")
        output.append(f"constexpr IdealCaloricNode {identifier}Nodes[] =")
        output.append("{")
        for row in grouped[species]:
            output.append(
                "    {"
                f"{format_float(float(row['temperature_k']))}, "
                f"{format_float(float(row['hsens_ideal_j_kg']))}, "
                f"{format_float(float(row['cp0_j_kg_k']))}"
                "},"
            )
        output.append("};")
        output.append("")
    return "\n".join(output)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail unless the checked-in CSV is byte-for-byte current",
    )
    parser.add_argument(
        "--emit-cpp",
        action="store_true",
        help="print C++ initializer rows for the audited built-in tables",
    )
    args = parser.parse_args()

    rows = build_rows()
    if args.emit_cpp:
        print(render_cpp(rows), end="")
        return 0

    rendered = render_csv(rows).encode("ascii")
    if args.check:
        require(args.output.exists(), f"missing frozen reference: {args.output}")
        require(
            args.output.read_bytes() == rendered,
            f"frozen reference is stale: {args.output}",
        )
        print(f"IDEAL_CALORIC_REFERENCE_CURRENT {args.output}")
        return 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(rendered)
    print(f"WROTE {args.output} rows={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
