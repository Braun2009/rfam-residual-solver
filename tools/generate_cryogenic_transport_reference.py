#!/usr/bin/env python3
"""Generate the frozen CoolProp 7.2.0 pure-fluid transport matrix."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
from typing import Any, Iterable

import CoolProp
from CoolProp import CoolProp as CP


EXPECTED_COOLPROP_VERSION = "7.2.0"
EXPECTED_COOLPROP_GIT_REVISION = "98b3523d5daa98454618d381d2ae53f7471d216b"
SCHEMA_VERSION = "rfam.cryogenic-transport-pure-reference/v1"
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "tests"
    / "reference"
    / "cryogenic_transport_pure_coolprop72_v1.csv"
)

FLUIDS = {
    "O2": {
        "coolprop_name": "Oxygen",
        "cas": "7782-44-7",
        "viscosity_reference": "Lemmon-IJT-2004",
        "conductivity_reference": "Lemmon-IJT-2004",
        "transport_json_sha256": (
            "17d25c4a10f8ccd5ec69d14f0fb1d332"
            "7a687423b6bdeea77ee1f79bff449b8a"
        ),
    },
    "C3H8": {
        "coolprop_name": "Propane",
        "cas": "74-98-6",
        "viscosity_reference": "Vogel-JPCRD-1998",
        "conductivity_reference": "Marsh-JCED-2002",
        "transport_json_sha256": (
            "c3a43f7b18a507f67f04e08fc28f353"
            "c0ae5a8370503cb21485b2298782afbfb"
        ),
    },
}

FIELDS = (
    "schema_version",
    "case_id",
    "species",
    "coolprop_fluid",
    "cas",
    "backend",
    "coolprop_version",
    "coolprop_git_revision",
    "input_mode",
    "temperature_k",
    "pressure_pa",
    "quality_mole",
    "phase",
    "rho_mass_kg_m3",
    "viscosity_pa_s",
    "conductivity_w_m_k",
    "viscosity_reference",
    "conductivity_reference",
    "transport_json_sha256",
    "source_api",
    "production_ready",
    "notes",
)

FLOAT_FIELDS = frozenset(
    {
        "temperature_k",
        "pressure_pa",
        "quality_mole",
        "rho_mass_kg_m3",
        "viscosity_pa_s",
        "conductivity_w_m_k",
    }
)

VALID_PT_PHASES = frozenset(
    {"gas", "liquid", "supercritical", "supercritical_gas", "supercritical_liquid"}
)

NOTES = (
    "pure-fluid-reference-only;not-mixture-or-effective-two-phase-transport;"
    "not-runtime-table;production-ready=0"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def format_float(value: float) -> str:
    require(math.isfinite(value), "attempted to serialize a non-finite value")
    return format(value, ".17e")


def transport_identity(fluid: str) -> dict[str, str]:
    payload = json.loads(CP.get_fluid_param_string(fluid, "JSON"))
    require(isinstance(payload, list) and len(payload) == 1, f"{fluid}: invalid JSON")
    transport = payload[0].get("TRANSPORT")
    require(isinstance(transport, dict), f"{fluid}: missing TRANSPORT JSON")
    canonical = json.dumps(
        transport, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    viscosity = transport.get("viscosity")
    conductivity = transport.get("conductivity")
    require(isinstance(viscosity, dict), f"{fluid}: missing viscosity JSON")
    require(isinstance(conductivity, dict), f"{fluid}: missing conductivity JSON")
    return {
        "viscosity_reference": str(viscosity.get("BibTeX", "")),
        "conductivity_reference": str(conductivity.get("BibTeX", "")),
        "transport_json_sha256": hashlib.sha256(canonical).hexdigest(),
    }


def point_specs() -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []

    for temperature in (110.0, 130.0):
        for pressure_mpa in (2.0, 4.0, 6.0):
            specs.append(
                {
                    "case_id": f"o2_pt_{int(temperature)}k_{int(pressure_mpa)}mpa",
                    "species": "O2",
                    "mode": "PT",
                    "temperature": temperature,
                    "pressure": pressure_mpa * 1.0e6,
                }
            )
    for temperature in (180.0, 300.0):
        for pressure_mpa in (2.0, 4.0, 6.0):
            specs.append(
                {
                    "case_id": f"o2_pt_{int(temperature)}k_{int(pressure_mpa)}mpa",
                    "species": "O2",
                    "mode": "PT",
                    "temperature": temperature,
                    "pressure": pressure_mpa * 1.0e6,
                }
            )

    for pressure_pa, pressure_label in ((2.0e6, "2mpa"), (4.0e6, "4mpa")):
        specs.extend(
            saturation_specs("O2", "o2", pressure_pa, pressure_label)
        )

    for scale, label in ((0.99, "below"), (1.01, "above")):
        specs.append(
            {
                "case_id": f"o2_near_critical_{label}",
                "species": "O2",
                "mode": "PT_CRITICAL_OFFSET",
                "temperature_scale": scale,
                "pressure_scale": 1.02,
            }
        )

    for pressure_pa, pressure_label in (
        (1.0e5, "0p1mpa"),
        (1.0e6, "1mpa"),
        (4.0e6, "4mpa"),
        (6.0e6, "6mpa"),
    ):
        specs.append(
            {
                "case_id": f"c3h8_pt_300k_{pressure_label}",
                "species": "C3H8",
                "mode": "PT",
                "temperature": 300.0,
                "pressure": pressure_pa,
            }
        )
    for pressure_pa, pressure_label in (
        (3.5e6, "3p5mpa"),
        (4.0e6, "4mpa"),
        (4.5e6, "4p5mpa"),
    ):
        specs.append(
            {
                "case_id": f"c3h8_pt_450k_{pressure_label}",
                "species": "C3H8",
                "mode": "PT",
                "temperature": 450.0,
                "pressure": pressure_pa,
            }
        )

    for pressure_pa, pressure_label in ((1.0e6, "1mpa"), (4.0e6, "4mpa")):
        specs.extend(
            saturation_specs("C3H8", "c3h8", pressure_pa, pressure_label)
        )

    for scale, label in ((0.99, "below"), (1.01, "above")):
        specs.append(
            {
                "case_id": f"c3h8_near_critical_{label}",
                "species": "C3H8",
                "mode": "PT_CRITICAL_OFFSET",
                "temperature_scale": scale,
                "pressure_scale": 1.02,
            }
        )

    return sorted(specs, key=lambda item: str(item["case_id"]))


def saturation_specs(
    species: str, prefix: str, pressure_pa: float, pressure_label: str
) -> list[dict[str, Any]]:
    return [
        {
            "case_id": f"{prefix}_sat_liquid_{pressure_label}",
            "species": species,
            "mode": "PQ",
            "pressure": pressure_pa,
            "quality": 0.0,
        },
        {
            "case_id": f"{prefix}_sat_vapor_{pressure_label}",
            "species": species,
            "mode": "PQ",
            "pressure": pressure_pa,
            "quality": 1.0,
        },
        {
            "case_id": f"{prefix}_sat_subcooled_0p25k_{pressure_label}",
            "species": species,
            "mode": "PT_SAT_OFFSET",
            "pressure": pressure_pa,
            "temperature_offset": -0.25,
        },
        {
            "case_id": f"{prefix}_sat_superheated_0p25k_{pressure_label}",
            "species": species,
            "mode": "PT_SAT_OFFSET",
            "pressure": pressure_pa,
            "temperature_offset": 0.25,
        },
    ]


def evaluate(spec: dict[str, Any]) -> dict[str, Any]:
    species = str(spec["species"])
    metadata = FLUIDS[species]
    fluid = str(metadata["coolprop_name"])
    mode = str(spec["mode"])
    quality: float | None = None

    if mode == "PQ":
        pressure = float(spec["pressure"])
        quality = float(spec["quality"])
        temperature = CP.PropsSI("T", "P", pressure, "Q", quality, fluid)
        rho = CP.PropsSI("DMASS", "P", pressure, "Q", quality, fluid)
        viscosity = CP.PropsSI("VISCOSITY", "P", pressure, "Q", quality, fluid)
        conductivity = CP.PropsSI(
            "CONDUCTIVITY", "P", pressure, "Q", quality, fluid
        )
        phase = "saturated_liquid" if quality == 0.0 else "saturated_vapor"
        input_mode = "PQ"
    else:
        if mode == "PT_SAT_OFFSET":
            pressure = float(spec["pressure"])
            saturation_temperature = CP.PropsSI("T", "P", pressure, "Q", 0.0, fluid)
            temperature = saturation_temperature + float(spec["temperature_offset"])
        elif mode == "PT_CRITICAL_OFFSET":
            critical_temperature = CP.PropsSI("TCRIT", fluid)
            critical_pressure = CP.PropsSI("PCRIT", fluid)
            temperature = critical_temperature * float(spec["temperature_scale"])
            pressure = critical_pressure * float(spec["pressure_scale"])
        else:
            temperature = float(spec["temperature"])
            pressure = float(spec["pressure"])
        rho = CP.PropsSI("DMASS", "T", temperature, "P", pressure, fluid)
        viscosity = CP.PropsSI("VISCOSITY", "T", temperature, "P", pressure, fluid)
        conductivity = CP.PropsSI(
            "CONDUCTIVITY", "T", temperature, "P", pressure, fluid
        )
        phase = CP.PhaseSI("T", temperature, "P", pressure, fluid)
        require(phase in VALID_PT_PHASES, f"{spec['case_id']}: unsupported phase {phase}")
        input_mode = "PT"

    for label, value in (
        ("temperature", temperature),
        ("pressure", pressure),
        ("density", rho),
        ("viscosity", viscosity),
        ("conductivity", conductivity),
    ):
        require(math.isfinite(value) and value > 0.0, f"{spec['case_id']}: invalid {label}")

    return {
        "schema_version": SCHEMA_VERSION,
        "case_id": spec["case_id"],
        "species": species,
        "coolprop_fluid": fluid,
        "cas": metadata["cas"],
        "backend": "HEOS",
        "coolprop_version": EXPECTED_COOLPROP_VERSION,
        "coolprop_git_revision": EXPECTED_COOLPROP_GIT_REVISION,
        "input_mode": input_mode,
        "temperature_k": temperature,
        "pressure_pa": pressure,
        "quality_mole": quality,
        "phase": phase,
        "rho_mass_kg_m3": rho,
        "viscosity_pa_s": viscosity,
        "conductivity_w_m_k": conductivity,
        "viscosity_reference": metadata["viscosity_reference"],
        "conductivity_reference": metadata["conductivity_reference"],
        "transport_json_sha256": metadata["transport_json_sha256"],
        "source_api": "CoolProp.PropsSI+PhaseSI+fluid-JSON",
        "production_ready": 0,
        "notes": NOTES,
    }


def build_rows() -> list[dict[str, Any]]:
    require(
        CoolProp.__version__ == EXPECTED_COOLPROP_VERSION,
        f"expected CoolProp {EXPECTED_COOLPROP_VERSION}, got {CoolProp.__version__}",
    )
    require(
        CP.get_global_param_string("version") == EXPECTED_COOLPROP_VERSION,
        "CoolProp library version mismatch",
    )
    require(
        CoolProp.__gitrevision__ == EXPECTED_COOLPROP_GIT_REVISION,
        "CoolProp git revision mismatch",
    )
    for metadata in FLUIDS.values():
        fluid = str(metadata["coolprop_name"])
        require(
            CP.get_fluid_param_string(fluid, "CAS") == metadata["cas"],
            f"{fluid}: CAS mismatch",
        )
        identity = transport_identity(fluid)
        for field in (
            "viscosity_reference",
            "conductivity_reference",
            "transport_json_sha256",
        ):
            require(
                identity[field] == metadata[field],
                f"{fluid}: transport identity mismatch for {field}",
            )
    rows = [evaluate(spec) for spec in point_specs()]
    require(len(rows) == 39, "reference point count changed")
    return rows


def render_csv(rows: Iterable[dict[str, Any]]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        serialized = dict(row)
        for field in FLOAT_FIELDS:
            value = serialized[field]
            serialized[field] = "" if value is None else format_float(float(value))
        writer.writerow(serialized)
    return stream.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail unless the frozen CSV is byte-for-byte current",
    )
    args = parser.parse_args()

    rendered = render_csv(build_rows()).encode("ascii")
    if args.check:
        require(args.output.exists(), f"missing frozen reference: {args.output}")
        require(
            args.output.read_bytes() == rendered,
            f"frozen reference is stale: {args.output}",
        )
        print(
            "CRYOGENIC_TRANSPORT_REFERENCE_CURRENT "
            f"path={args.output} rows=39 production_ready=0"
        )
        return 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(rendered)
    print(
        "CRYOGENIC_TRANSPORT_REFERENCE_WROTE "
        f"path={args.output} rows=39 production_ready=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
