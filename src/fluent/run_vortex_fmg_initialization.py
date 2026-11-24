#!/usr/bin/env python3
"""Run an isolated Fluent FMG initialization for the Vortex cold-air case.

Reads one immutable ideal-air SST configuration case, verifies topology and
models, changes only the FMG initialization turbulent-viscosity ratio in the
new session, executes Fluent's ``fmg_initialize`` command, and writes fresh
case/data into ``--work``.  It does not iterate, modify the mesh, or overwrite
the source case.

This is a numerical initialization diagnostic only.  It makes no real-fluid,
LOX/C3H8, chemistry, EDC, thrust, or heat-transfer claim.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
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


def verify(settings: Any, farfield_type: str, requested_fmg_ratio: float) -> dict[str, Any]:
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
    result["solution_methods"] = settings.solution.methods.get_state()
    result["solution_controls"] = settings.solution.controls.get_state()
    result["initialization"] = settings.solution.initialization.get_state()
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
    actual_ratio = float(result["initialization"]["fmg"]["turbulent_viscosity_ratio"])
    if abs(actual_ratio - requested_fmg_ratio) > max(1e-12, abs(requested_fmg_ratio) * 1e-12):
        raise RuntimeError(f"FMG ratio readback mismatch: {actual_ratio} != {requested_fmg_ratio}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--processors", type=int, default=36)
    parser.add_argument(
        "--farfield-type", choices=tuple(FARFIELD_GROUP_BY_TYPE), required=True,
        help="actual Fluent type of the retained farfield physical zone",
    )
    parser.add_argument("--fmg-turbulent-viscosity-ratio", type=float, default=10.0)
    args = parser.parse_args()
    if args.processors < 1:
        raise ValueError("--processors must be positive")
    if not math.isfinite(args.fmg_turbulent_viscosity_ratio) or args.fmg_turbulent_viscosity_ratio <= 0.0:
        raise ValueError("--fmg-turbulent-viscosity-ratio must be finite and positive")
    source = args.case.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    output_case = work / "cold_air_sst_fmg_initialized.cas.h5"
    output_data = work / "cold_air_sst_fmg_initialized.dat.h5"
    document: dict[str, Any] = {
        "schema": "vortex.fluent.fmg-initialization/v1",
        "started_utc": stamp(),
        "status": "starting",
        "source_case": str(source),
        "source_is_immutable": True,
        "processors": args.processors,
        "farfield_boundary_type": args.farfield_type,
        "requested_fmg_turbulent_viscosity_ratio": args.fmg_turbulent_viscosity_ratio,
        "actions": [
            "read case", "verify topology/models", "set FMG turbulence ratio",
            "FMG initialize", "write independent case/data",
        ],
        "prohibited_actions": ["iteration", "mesh modification", "source-case write"],
        "scope": "ideal-air density-based implicit SST cold-flow initialization diagnostic only",
    }
    dump(work / "manifest.json", document)
    session = None
    try:
        session = pyfluent.launch_fluent(
            precision="double", processor_count=args.processors, mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI, start_timeout=1200,
        )
        emit("FMG_INIT_SOLVER_LAUNCHED")
        settings = session.settings
        settings.file.read_case(file_name=str(source))
        emit("FMG_INIT_SOURCE_CASE_READ_OK")
        initialization = settings.solution.initialization
        document["fmg_before"] = initialization.fmg.get_state()
        initialization.fmg.turbulent_viscosity_ratio = args.fmg_turbulent_viscosity_ratio
        document["readback_before_fmg"] = verify(
            settings, args.farfield_type, args.fmg_turbulent_viscosity_ratio
        )
        document["status"] = "configured"
        dump(work / "manifest.json", document)
        initialization.fmg.fmg_initialize()
        emit("FMG_INIT_EXECUTED_OK")
        document["readback_after_fmg"] = verify(
            settings, args.farfield_type, args.fmg_turbulent_viscosity_ratio
        )
        document["output_case"] = str(output_case)
        document["output_data"] = str(output_data)
        document["status"] = "initialized"
        dump(work / "manifest.json", document)
        settings.file.write_case(file_name=str(output_case))
        settings.file.write_data(file_name=str(output_data))
        document["status"] = "ok"
        document["finished_utc"] = stamp()
        dump(work / "manifest.json", document)
        emit(f"FMG_INIT_CASE_DATA_WRITTEN case={output_case} data={output_data}")
        emit("VORTEX_FMG_INITIALIZATION_OK")
    except Exception:
        document["status"] = "failed"
        document["finished_utc"] = stamp()
        document["traceback"] = traceback.format_exc()
        dump(work / "manifest.json", document)
        raise
    finally:
        if session is not None:
            try:
                session.exit()
                emit("FMG_INIT_SOLVER_EXIT_OK")
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
