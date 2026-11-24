#!/usr/bin/env python3
"""Create an independent Vortex cold-air configuration with an open pressure outlet.

The immutable r2 configured case uses a pressure-far-field boundary.  Fluent
requires its initialization Mach number to be strictly positive, which is not
an appropriate way to represent a quiescent external ambient in this cold-air
commissioning test.  This utility reads r2, changes only the ``farfield``
boundary type to ``pressure-outlet`` and supplies its ambient backflow state.

It never initializes, iterates, writes the source case, or modifies the final
mesh.  The resulting case is an explicitly separate *ideal-air commissioning*
branch, not a real-fluid/reacting production configuration.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sys
import traceback
from typing import Any

os.environ.setdefault("AWP_ROOT261", "/ansys_inc/v261")
os.environ.setdefault("ANSYSLMD_LICENSE_FILE", "1055@localhost")
os.environ.setdefault("ANSYSCL_FNP_PATH", os.environ["ANSYSLMD_LICENSE_FILE"])

import ansys.fluent.core as pyfluent
from ansys.fluent.core import FluentMode, UIMode

SOURCE = Path(
    "/root/vortex_fluent/gmsh_fixed_v2/solver_cold_air_sst_commission_20260815_r2/"
    "cold_air_sst_configured.cas.h5"
)
EXPECTED_FLUID = {
    "vortex_full360_v2_zoned_ascii-prism-cells",
    "vortex_full360_v2_zoned_ascii-tet-cells",
}


def stamp() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def log(message: str) -> None:
    print(message, flush=True)


def dump(path: Path, document: dict[str, Any]) -> None:
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def verify(settings: Any, temperature_k: float) -> dict[str, Any]:
    setup = settings.setup
    bc = setup.boundary_conditions
    expected_groups = {
        "mass_flow_inlet": {"inlet_spray", "inlet_torch"},
        "pressure_far_field": set(),
        "pressure_outlet": {"farfield"},
        "wall": {"wall_chamber", "wall_plate"},
    }
    result: dict[str, Any] = {}
    for group_name, expected in expected_groups.items():
        actual = set(getattr(bc, group_name).get_object_names())
        result[group_name] = sorted(actual)
        if actual != expected:
            raise RuntimeError(
                f"{group_name}: expected {sorted(expected)}, got {sorted(actual)}"
            )
    fluid = set(setup.cell_zone_conditions.fluid.get_object_names())
    result["fluid"] = sorted(fluid)
    if fluid != EXPECTED_FLUID:
        raise RuntimeError(f"fluid zones: expected {sorted(EXPECTED_FLUID)}, got {sorted(fluid)}")

    result["solver"] = setup.general.solver.get_state()
    result["air"] = setup.materials.fluid["air"].get_state()
    result["energy"] = setup.models.energy.get_state()
    result["viscous"] = setup.models.viscous.get_state()
    result["outlet"] = bc.pressure_outlet["farfield"].get_state()
    result["initialization_type"] = settings.solution.initialization.initialization_type()
    if result["solver"]["type"] != "density-based-implicit":
        raise RuntimeError(f"solver mismatch: {result['solver']}")
    if result["air"]["density"]["option"] != "ideal-gas":
        raise RuntimeError(f"air-density mismatch: {result['air']}")
    if not result["energy"]["enabled"]:
        raise RuntimeError("energy equation disabled")
    if result["viscous"]["model"] != "k-omega" or result["viscous"]["k_omega_model"] != "sst":
        raise RuntimeError(f"SST mismatch: {result['viscous']}")
    if result["initialization_type"] != "standard":
        raise RuntimeError(f"initialization type mismatch: {result['initialization_type']}")
    outlet = result["outlet"]
    if abs(outlet["momentum"]["gauge_pressure"]["value"]) > 1.0e-12:
        raise RuntimeError(f"farfield outlet gauge pressure changed: {outlet}")
    if abs(outlet["thermal"]["backflow_total_temperature"]["value"] - temperature_k) > 1.0e-12:
        raise RuntimeError(f"farfield backflow temperature mismatch: {outlet}")
    if outlet["momentum"]["prevent_reverse_flow"]:
        raise RuntimeError(f"reverse flow must remain allowed at open farfield: {outlet}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--processors", type=int, default=4)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    args = parser.parse_args()
    if args.processors < 1:
        raise ValueError("--processors must be positive")
    if args.temperature_k <= 0.0:
        raise ValueError("--temperature-k must be positive")
    if not SOURCE.is_file():
        raise FileNotFoundError(SOURCE)

    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    destination = work / "cold_air_sst_pressure_outlet_configured.cas.h5"
    document: dict[str, Any] = {
        "schema": "vortex.fluent.farfield-pressure-outlet-configuration/v1",
        "started_utc": stamp(),
        "status": "starting",
        "source_case": str(SOURCE),
        "source_is_immutable": True,
        "processors": args.processors,
        "boundary_change": {
            "zone": "farfield",
            "old_type": "pressure-far-field",
            "new_type": "pressure-outlet",
            "gauge_pressure_pa": 0.0,
            "backflow_total_temperature_k": args.temperature_k,
            "backflow_turbulence": "Intensity and Viscosity Ratio; I=0.05, mu_t/mu=10",
            "prevent_reverse_flow": False,
        },
        "scope": "ideal-air density-based SST cold-flow commissioning only; no GERG, reaction, EDC, or LOX/C3H8 performance claim",
        "prohibited_actions": ["initialization", "iteration", "source-case write", "mesh modification"],
    }
    dump(work / "configuration.json", document)
    session = None
    try:
        session = pyfluent.launch_fluent(
            precision="double", processor_count=args.processors,
            mode=FluentMode.SOLVER, ui_mode=UIMode.NO_GUI, start_timeout=1200,
        )
        log("PRESSURE_OUTLET_CONFIG_SOLVER_LAUNCHED")
        settings = session.settings
        settings.file.read_case(file_name=str(SOURCE))
        log("PRESSURE_OUTLET_CONFIG_SOURCE_CASE_READ_OK")
        settings.mesh.modify_zones.zone_type(
            zone_names=["farfield"], new_type="pressure-outlet"
        )
        outlet = settings.setup.boundary_conditions.pressure_outlet["farfield"]
        outlet.momentum.gauge_pressure.value = 0.0
        outlet.momentum.prevent_reverse_flow = False
        outlet.thermal.backflow_total_temperature.value = args.temperature_k
        outlet.turbulence.turbulence_specification = "Intensity and Viscosity Ratio"
        outlet.turbulence.backflow_turbulent_intensity = 0.05
        outlet.turbulence.backflow_turbulent_viscosity_ratio = 10.0
        settings.solution.initialization.initialization_type = "standard"

        document["readback"] = verify(settings, args.temperature_k)
        document["output_case"] = str(destination)
        document["status"] = "configured"
        dump(work / "configuration.json", document)
        settings.file.write_case(file_name=str(destination))
        document["status"] = "ok"
        document["finished_utc"] = stamp()
        dump(work / "configuration.json", document)
        log(f"PRESSURE_OUTLET_CONFIGURED_CASE_WRITTEN path={destination}")
        log("VORTEX_FARFIELD_PRESSURE_OUTLET_CONFIG_OK")
    except Exception:
        document["status"] = "failed"
        document["finished_utc"] = stamp()
        document["traceback"] = traceback.format_exc()
        dump(work / "configuration.json", document)
        raise
    finally:
        if session is not None:
            try:
                session.exit()
                log("PRESSURE_OUTLET_CONFIG_SOLVER_EXIT_OK")
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
