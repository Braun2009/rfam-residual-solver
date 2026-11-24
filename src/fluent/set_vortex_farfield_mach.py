#!/usr/bin/env python3
"""Write an independent Vortex cold-air configuration with valid far-field Mach.

Fluent rejects pressure-far-field Mach=0 during initialization.  This utility
reads the immutable r2 configured case, changes only the far-field Mach and
sets standard initialization explicitly, then writes a new case.  The final
mesh and source case are never modified.
"""
from __future__ import annotations

import argparse
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
EXPECTED_BOUNDARIES = {
    "mass_flow_inlet": {"inlet_spray", "inlet_torch"},
    "pressure_far_field": {"farfield"},
    "wall": {"wall_chamber", "wall_plate"},
}
EXPECTED_FLUID = {
    "vortex_full360_v2_zoned_ascii-prism-cells",
    "vortex_full360_v2_zoned_ascii-tet-cells",
}


def log(message: str) -> None:
    print(message, flush=True)


def verify(settings: Any, mach: float) -> dict[str, Any]:
    setup = settings.setup
    bc = setup.boundary_conditions
    result: dict[str, Any] = {}
    for group_name, expected in EXPECTED_BOUNDARIES.items():
        actual = set(getattr(bc, group_name).get_object_names())
        result[group_name] = sorted(actual)
        if actual != expected:
            raise RuntimeError(f"{group_name}: expected {sorted(expected)}, got {sorted(actual)}")
    fluid = set(setup.cell_zone_conditions.fluid.get_object_names())
    result["fluid"] = sorted(fluid)
    if fluid != EXPECTED_FLUID:
        raise RuntimeError(f"fluid zones: expected {sorted(EXPECTED_FLUID)}, got {sorted(fluid)}")

    result["solver"] = setup.general.solver.get_state()
    result["air"] = setup.materials.fluid["air"].get_state()
    result["energy"] = setup.models.energy.get_state()
    result["viscous"] = setup.models.viscous.get_state()
    result["farfield"] = bc.pressure_far_field["farfield"].get_state()
    result["initialization_type"] = settings.solution.initialization.initialization_type()

    if result["solver"]["type"] != "density-based-implicit":
        raise RuntimeError(f"solver mismatch: {result['solver']}")
    if result["air"]["density"]["option"] != "ideal-gas":
        raise RuntimeError(f"air density mismatch: {result['air']}")
    if not result["energy"]["enabled"]:
        raise RuntimeError("energy equation disabled")
    if result["viscous"]["model"] != "k-omega" or result["viscous"]["k_omega_model"] != "sst":
        raise RuntimeError(f"SST mismatch: {result['viscous']}")
    actual_mach = result["farfield"]["momentum"]["mach_number"]["value"]
    if abs(actual_mach - mach) > 1.0e-12:
        raise RuntimeError(f"farfield Mach readback mismatch: {actual_mach} != {mach}")
    if result["initialization_type"] != "standard":
        raise RuntimeError(f"initialization type mismatch: {result['initialization_type']}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--processors", type=int, default=4)
    parser.add_argument("--farfield-mach", type=float, default=0.01)
    args = parser.parse_args()
    if args.processors < 1:
        raise ValueError("processors must be positive")
    if not 0.0 < args.farfield_mach < 1.0:
        raise ValueError("farfield Mach must be in (0, 1)")
    if not SOURCE.is_file():
        raise FileNotFoundError(SOURCE)

    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    session = None
    try:
        session = pyfluent.launch_fluent(
            precision="double", processor_count=args.processors,
            mode=FluentMode.SOLVER, ui_mode=UIMode.NO_GUI, start_timeout=1200,
        )
        log("FARFIELD_MACH_CONFIG_SOLVER_LAUNCHED")
        settings = session.settings
        settings.file.read_case(file_name=str(SOURCE))
        log("FARFIELD_MACH_CONFIG_SOURCE_CASE_READ_OK")
        farfield = settings.setup.boundary_conditions.pressure_far_field["farfield"]
        old_mach = farfield.momentum.mach_number.value()
        farfield.momentum.mach_number.value = args.farfield_mach
        settings.solution.initialization.initialization_type = "standard"
        readback = verify(settings, args.farfield_mach)
        document = {
            "source_case": str(SOURCE),
            "processors": args.processors,
            "old_farfield_mach": old_mach,
            "new_farfield_mach": args.farfield_mach,
            "note": "Fluent pressure-far-field requires strictly positive Mach for initialization.",
            "readback": readback,
        }
        (work / "configuration.json").write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        # Keep the artifact name consistent with the actual strictly-positive
        # Mach value; do not reuse the legacy mach001 filename for a different
        # far-field condition.
        mach_token = format(args.farfield_mach, ".12g")
        configured_case = work / f"cold_air_sst_mach{mach_token}_configured.cas.h5"
        settings.file.write_case(file_name=str(configured_case))
        document["configured_case"] = str(configured_case)
        (work / "configuration.json").write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        log(f"FARFIELD_MACH_CONFIGURED_CASE_WRITTEN path={configured_case}")
        log("VORTEX_FARFIELD_MACH_CONFIG_OK")
    finally:
        if session is not None:
            session.exit()
            log("FARFIELD_MACH_CONFIG_SOLVER_EXIT_OK")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
