#!/usr/bin/env python3
"""Create an immutable-source Vortex cold-air restart with a new Courant number.

The utility reads a solved case/data pair, validates that it remains the
ideal-air density-based implicit SST commissioning setup, changes *only* the
solution Courant number, and writes a fresh case/data pair in ``--work``.
It never initializes or iterates and never overwrites the source checkpoint.
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

BASE_EXPECTED_BOUNDARIES = {
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


def dump(path: Path, doc: dict[str, Any]) -> None:
    path.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def verify(settings: Any, requested_courant: float, farfield_type: str) -> dict[str, Any]:
    setup = settings.setup
    bc = setup.boundary_conditions
    farfield_group = FARFIELD_GROUP_BY_TYPE[farfield_type]
    expected_boundaries = {
        **BASE_EXPECTED_BOUNDARIES,
        "pressure_far_field": {"farfield"} if farfield_group == "pressure_far_field" else set(),
        "pressure_outlet": {"farfield"} if farfield_group == "pressure_outlet" else set(),
    }
    out: dict[str, Any] = {}
    for group, expected in expected_boundaries.items():
        actual = set(getattr(bc, group).get_object_names())
        out[group] = sorted(actual)
        if actual != expected:
            raise RuntimeError(f"{group}: expected {sorted(expected)}, got {sorted(actual)}")
    fluid = set(setup.cell_zone_conditions.fluid.get_object_names())
    out["fluid"] = sorted(fluid)
    if fluid != EXPECTED_FLUID:
        raise RuntimeError(f"fluid zones: expected {sorted(EXPECTED_FLUID)}, got {sorted(fluid)}")

    out["solver"] = setup.general.solver.get_state()
    out["energy"] = setup.models.energy.get_state()
    out["viscous"] = setup.models.viscous.get_state()
    out["air"] = setup.materials.fluid["air"].get_state()
    out["farfield"] = getattr(bc, farfield_group)["farfield"].get_state()
    out["farfield_boundary_type"] = farfield_type
    out["solution_controls"] = settings.solution.controls.get_state()
    if out["solver"]["type"] != "density-based-implicit":
        raise RuntimeError(f"solver mismatch: {out['solver']}")
    if not out["energy"]["enabled"]:
        raise RuntimeError("energy is disabled")
    if out["viscous"]["model"] != "k-omega" or out["viscous"]["k_omega_model"] != "sst":
        raise RuntimeError(f"SST mismatch: {out['viscous']}")
    if out["air"]["density"]["option"] != "ideal-gas":
        raise RuntimeError(f"air-density mismatch: {out['air']}")
    actual = float(out["solution_controls"]["courant_number"])
    if abs(actual - requested_courant) > max(1.0e-12, abs(requested_courant) * 1.0e-12):
        raise RuntimeError(f"Courant readback mismatch: {actual} != {requested_courant}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--processors", type=int, default=4)
    parser.add_argument("--courant", type=float, required=True)
    parser.add_argument(
        "--farfield-type", choices=tuple(FARFIELD_GROUP_BY_TYPE),
        default="pressure-far-field",
        help="actual Fluent type of the retained farfield physical zone",
    )
    args = parser.parse_args()
    if args.processors < 1:
        raise ValueError("--processors must be positive")
    if not 0.0 < args.courant <= 100.0:
        raise ValueError("--courant must be in (0, 100]")
    source_case, source_data = args.case.resolve(), args.data.resolve()
    if not source_case.is_file() or not source_data.is_file():
        raise FileNotFoundError(f"missing source case/data: {source_case} / {source_data}")
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    token = format(args.courant, ".12g").replace(".", "p")
    result_case = work / f"cold_air_sst_courant{token}_configured.cas.h5"
    result_data = work / f"cold_air_sst_courant{token}_configured.dat.h5"
    document: dict[str, Any] = {
        "schema": "vortex.fluent.courant-restart-configuration/v1",
        "started_utc": stamp(),
        "status": "starting",
        "source_case": str(source_case),
        "source_data": str(source_data),
        "source_is_immutable": True,
        "processors": args.processors,
        "requested_courant": args.courant,
        "farfield_boundary_type": args.farfield_type,
        "actions": ["read case", "read data", "set Courant", "write independent case/data"],
        "prohibited_actions": ["initialization", "iteration", "mesh modification"],
    }
    dump(work / "configuration.json", document)
    session = None
    try:
        session = pyfluent.launch_fluent(
            precision="double", processor_count=args.processors, mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI, start_timeout=1200,
        )
        emit("COURANT_CONFIG_SOLVER_LAUNCHED")
        settings = session.settings
        settings.file.read_case(file_name=str(source_case))
        settings.file.read_data(file_name=str(source_data))
        emit("COURANT_CONFIG_SOURCE_CASE_DATA_READ_OK")
        before = settings.solution.controls.get_state()
        settings.solution.controls.courant_number = args.courant
        document["old_courant"] = before["courant_number"]
        document["readback"] = verify(settings, args.courant, args.farfield_type)
        document["output_case"] = str(result_case)
        document["output_data"] = str(result_data)
        document["status"] = "configured"
        dump(work / "configuration.json", document)
        settings.file.write_case(file_name=str(result_case))
        settings.file.write_data(file_name=str(result_data))
        emit(f"COURANT_CONFIG_CASE_DATA_WRITTEN case={result_case} data={result_data}")
        document["status"] = "ok"
        document["finished_utc"] = stamp()
        dump(work / "configuration.json", document)
        emit("VORTEX_COURANT_CONFIG_OK")
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
                emit("COURANT_CONFIG_SOLVER_EXIT_OK")
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
