#!/usr/bin/env python3
"""Resume the configured Vortex 3D cold-air SST case on a new MPI allocation.

Uses the preserved ``cold_air_sst_configured.cas.h5`` from commissioning r2;
therefore it does not reread or modify the final mesh file.  Standard
initialization is intentional: the prior hybrid initialization entered Fluent
26.1's interrupted-compute state before producing initialized data.
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

SOURCE_CASE = Path(
    "/root/vortex_fluent/gmsh_fixed_v2/solver_cold_air_sst_commission_20260815_r2/"
    "cold_air_sst_configured.cas.h5"
)
REQUIRED_BOUNDARIES = {
    "mass_flow_inlet": {"inlet_spray", "inlet_torch"},
    "pressure_far_field": {"farfield"},
    "wall": {"wall_chamber", "wall_plate"},
}
REQUIRED_FLUID_ZONES = {
    "vortex_full360_v2_zoned_ascii-prism-cells",
    "vortex_full360_v2_zoned_ascii-tet-cells",
}


def stamp() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def log(message: str) -> None:
    print(message, flush=True)


def dump(path: Path, doc: dict[str, Any]) -> None:
    path.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def check(settings: Any) -> dict[str, Any]:
    setup = settings.setup
    bc = setup.boundary_conditions
    result: dict[str, Any] = {}
    for group_name, expected in REQUIRED_BOUNDARIES.items():
        actual = set(getattr(bc, group_name).get_object_names())
        result[group_name] = sorted(actual)
        if actual != expected:
            raise RuntimeError(f"{group_name}: expected {sorted(expected)}, got {sorted(actual)}")
    fluid = set(setup.cell_zone_conditions.fluid.get_object_names())
    result["fluid"] = sorted(fluid)
    if fluid != REQUIRED_FLUID_ZONES:
        raise RuntimeError(f"fluid zones: expected {sorted(REQUIRED_FLUID_ZONES)}, got {sorted(fluid)}")
    solver = setup.general.solver.get_state()
    air = setup.materials.fluid["air"].get_state()
    viscous = setup.models.viscous.get_state()
    result["solver"] = solver
    result["energy"] = setup.models.energy.get_state()
    result["air"] = air
    result["viscous"] = viscous
    result["inlets"] = {
        name: bc.mass_flow_inlet[name].get_state()
        for name in ("inlet_spray", "inlet_torch")
    }
    result["farfield"] = bc.pressure_far_field["farfield"].get_state()
    if solver["type"] != "density-based-implicit":
        raise RuntimeError(f"wrong solver: {solver}")
    if not result["energy"]["enabled"]:
        raise RuntimeError("energy disabled")
    if air["density"]["option"] != "ideal-gas":
        raise RuntimeError(f"air density wrong: {air}")
    if viscous["model"] != "k-omega" or viscous["k_omega_model"] != "sst":
        raise RuntimeError(f"SST missing: {viscous}")
    return result


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--processors", type=int, default=36)
    p.add_argument("--iterations", type=int, default=5)
    p.add_argument("--courant", type=float, default=0.5)
    return p.parse_args()


def main() -> None:
    a = args()
    if a.processors < 1 or a.iterations < 1 or a.courant <= 0:
        raise ValueError("processors/iterations/courant must be positive")
    if not SOURCE_CASE.is_file():
        raise FileNotFoundError(SOURCE_CASE)
    work = a.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, Any] = {
        "schema": "vortex.fluent.cold-air-sst-resume/v1",
        "started_utc": stamp(),
        "status": "starting",
        "source_configured_case": str(SOURCE_CASE),
        "intent": "ideal-air cold-flow commissioning; no species, reactions, GERG, or EDC",
        "processors": a.processors,
        "iterations": a.iterations,
        "courant": a.courant,
        "initialization": "standard",
    }
    dump(work / "manifest.json", manifest)
    s = None
    try:
        s = pyfluent.launch_fluent(
            precision="double", processor_count=a.processors, mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI, start_timeout=1200,
        )
        log("COLD_AIR_RESUME_SOLVER_LAUNCHED")
        s.settings.file.read_case(file_name=str(SOURCE_CASE))
        log("COLD_AIR_CONFIGURED_CASE_READ_OK")
        s.settings.solution.controls.courant_number = a.courant
        manifest["readback_before_initialization"] = check(s.settings)
        manifest["status"] = "case_read"
        dump(work / "manifest.json", manifest)

        s.settings.solution.initialization.initialization_type = "standard"
        s.settings.solution.initialization.standard_initialize()
        log("COLD_AIR_STANDARD_INITIALIZATION_OK")
        s.settings.file.write_case(file_name=str(work / "cold_air_sst_initialized.cas.h5"))
        s.settings.file.write_data(file_name=str(work / "cold_air_sst_initialized.dat.h5"))
        log("COLD_AIR_INITIALIZED_CASE_DATA_WRITTEN")

        s.settings.solution.run_calculation.iterate(iter_count=a.iterations)
        log(f"COLD_AIR_ITERATION_OK count={a.iterations}")
        manifest["readback_after_iterations"] = check(s.settings)
        manifest["status"] = "iterated"
        manifest["finished_utc"] = stamp()
        dump(work / "manifest.json", manifest)
        s.settings.file.write_case(file_name=str(work / "cold_air_sst_commissioning.cas.h5"))
        s.settings.file.write_data(file_name=str(work / "cold_air_sst_commissioning.dat.h5"))
        log("COLD_AIR_FINAL_CASE_DATA_WRITTEN")
        log("VORTEX_3D_COLD_AIR_SST_RESUME_OK")
    except Exception:
        manifest["status"] = "failed"
        manifest["finished_utc"] = stamp()
        manifest["traceback"] = traceback.format_exc()
        dump(work / "manifest.json", manifest)
        raise
    finally:
        if s is not None:
            try:
                s.exit()
                log("COLD_AIR_RESUME_SOLVER_EXIT_OK")
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
