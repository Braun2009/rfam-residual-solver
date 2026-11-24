#!/usr/bin/env python3
"""Create an immutable-source Vortex restart using density-based global pseudo time.

This utility changes only Fluent 26.1's density-based pseudo-time formulation
from ``off`` to ``global-time-step``.  It is intended as an isolated low-Mach
commissioning diagnostic for an already configured ideal-air SST case.  It
never initializes, iterates, changes the mesh, writes data, or overwrites its
source case.
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

BASE_EXPECTED = {
    "mass_flow_inlet": {"inlet_spray", "inlet_torch"},
    "wall": {"wall_chamber", "wall_plate"},
}
FARFIELD_GROUP_BY_TYPE = {
    "pressure-far-field": "pressure_far_field",
    "pressure-outlet": "pressure_outlet",
}
EXPECTED_FLUID = {
    "vortex_full360_v2_zoned_ascii-prism-cells",
    "vortex_full360_v2_zoned_ascii-tet-cells",
}


def stamp() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def emit(message: str) -> None:
    print(message, flush=True)


def dump(path: Path, document: dict[str, Any]) -> None:
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def verify(settings: Any, farfield_type: str) -> dict[str, Any]:
    setup = settings.setup
    bc = setup.boundary_conditions
    farfield_group = FARFIELD_GROUP_BY_TYPE[farfield_type]
    expected_boundaries = {
        **BASE_EXPECTED,
        "pressure_far_field": {"farfield"} if farfield_group == "pressure_far_field" else set(),
        "pressure_outlet": {"farfield"} if farfield_group == "pressure_outlet" else set(),
    }
    result: dict[str, Any] = {}
    for group, expected_names in expected_boundaries.items():
        actual_names = set(getattr(bc, group).get_object_names())
        result[group] = sorted(actual_names)
        if actual_names != expected_names:
            raise RuntimeError(f"{group}: expected {sorted(expected_names)}, got {sorted(actual_names)}")

    fluid = set(setup.cell_zone_conditions.fluid.get_object_names())
    result["fluid"] = sorted(fluid)
    if fluid != EXPECTED_FLUID:
        raise RuntimeError(f"fluid zones: expected {sorted(EXPECTED_FLUID)}, got {sorted(fluid)}")

    result["solver"] = setup.general.solver.get_state()
    result["energy"] = setup.models.energy.get_state()
    result["viscous"] = setup.models.viscous.get_state()
    result["air"] = setup.materials.fluid["air"].get_state()
    result["farfield"] = getattr(bc, farfield_group)["farfield"].get_state()
    result["initialization"] = settings.solution.initialization.get_state()
    result["solution_methods"] = settings.solution.methods.get_state()
    result["solution_controls"] = settings.solution.controls.get_state()
    result["wall_distance_method"] = settings.mesh.wall_distance_method()

    if result["solver"]["type"] != "density-based-implicit":
        raise RuntimeError(f"solver mismatch: {result['solver']}")
    if result["solver"]["time"] != "steady":
        raise RuntimeError(f"time formulation mismatch: {result['solver']}")
    if not result["energy"]["enabled"]:
        raise RuntimeError("energy is disabled")
    if result["air"]["density"]["option"] != "ideal-gas":
        raise RuntimeError(f"air density mismatch: {result['air']}")
    if result["viscous"]["model"] != "k-omega" or result["viscous"]["k_omega_model"] != "sst":
        raise RuntimeError(f"SST mismatch: {result['viscous']}")
    pseudo_time = result["solution_methods"]["pseudo_time_method"]["formulation"]
    if pseudo_time["density_based_solver"] != "global-time-step":
        raise RuntimeError(f"global pseudo time setting lost: {pseudo_time}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--processors", type=int, default=4)
    parser.add_argument(
        "--farfield-type", choices=tuple(FARFIELD_GROUP_BY_TYPE),
        required=True,
        help="actual Fluent type of the retained farfield physical zone",
    )
    args = parser.parse_args()
    if args.processors < 1:
        raise ValueError("--processors must be positive")
    source = args.case.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    output_case = work / "cold_air_sst_global_timestep_configured.cas.h5"
    document: dict[str, Any] = {
        "schema": "vortex.fluent.density-based-global-time-step-configuration/v1",
        "started_utc": stamp(),
        "status": "starting",
        "source_case": str(source),
        "source_is_immutable": True,
        "processors": args.processors,
        "farfield_boundary_type": args.farfield_type,
        "old_density_based_pseudo_time": "read after source case load",
        "new_density_based_pseudo_time": "global-time-step",
        "actions": [
            "read case",
            "verify topology/models",
            "set density-based pseudo time",
            "write independent case",
        ],
        "prohibited_actions": [
            "initialization",
            "iteration",
            "mesh modification",
            "write_data",
            "source-case write",
        ],
    }
    dump(work / "configuration.json", document)
    session = None
    try:
        session = pyfluent.launch_fluent(
            precision="double", processor_count=args.processors, mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI, start_timeout=1200,
        )
        emit("GLOBAL_TIMESTEP_CONFIG_SOLVER_LAUNCHED")
        settings = session.settings
        settings.file.read_case(file_name=str(source))
        emit("GLOBAL_TIMESTEP_CONFIG_SOURCE_CASE_READ_OK")
        before = settings.solution.methods.pseudo_time_method.formulation.density_based_solver()
        settings.solution.methods.pseudo_time_method.formulation.density_based_solver = "global-time-step"
        document["old_density_based_pseudo_time"] = before
        document["readback"] = verify(settings, args.farfield_type)
        document["output_case"] = str(output_case)
        document["status"] = "configured"
        dump(work / "configuration.json", document)
        settings.file.write_case(file_name=str(output_case))
        document["status"] = "ok"
        document["finished_utc"] = stamp()
        dump(work / "configuration.json", document)
        emit(f"GLOBAL_TIMESTEP_CONFIG_CASE_WRITTEN path={output_case}")
        emit("VORTEX_GLOBAL_TIMESTEP_CONFIG_OK")
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
                emit("GLOBAL_TIMESTEP_CONFIG_SOLVER_EXIT_OK")
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
