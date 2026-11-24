#!/usr/bin/env python3
"""Create an immutable-source Vortex case using transport-equation wall distance.

Fluent's geometric wall-distance initialization previously deadlocked in
parallel on the large wedge/tet mesh.  This utility changes only
``mesh.wall_distance_method`` on a supplied ideal-air SST checkpoint, verifies
all five physical boundaries (including either supported farfield type), and
writes a new case.  It does not alter the source case, data, or mesh.
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

DEFAULT_SOURCE = Path(
    "/root/vortex_fluent/gmsh_fixed_v2/solver_cold_air_sst_commission_20260815_r2/"
    "cold_air_sst_configured.cas.h5"
)
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


def log(message: str) -> None:
    print(message, flush=True)


def dump(path: Path, document: dict[str, Any]) -> None:
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def verify_topology(settings: Any, farfield_type: str) -> dict[str, Any]:
    bc = settings.setup.boundary_conditions
    farfield_group = FARFIELD_GROUP_BY_TYPE[farfield_type]
    expected = {
        **BASE_EXPECTED,
        "pressure_far_field": {"farfield"} if farfield_group == "pressure_far_field" else set(),
        "pressure_outlet": {"farfield"} if farfield_group == "pressure_outlet" else set(),
    }
    result: dict[str, Any] = {}
    for group_name, required_names in expected.items():
        actual = set(getattr(bc, group_name).get_object_names())
        result[group_name] = sorted(actual)
        if actual != required_names:
            raise RuntimeError(f"{group_name}: expected {sorted(required_names)}, got {sorted(actual)}")
    fluid = set(settings.setup.cell_zone_conditions.fluid.get_object_names())
    result["fluid"] = sorted(fluid)
    if fluid != EXPECTED_FLUID:
        raise RuntimeError(f"fluid zones: expected {sorted(EXPECTED_FLUID)}, got {sorted(fluid)}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--processors", type=int, default=4)
    parser.add_argument(
        "--farfield-type", choices=tuple(FARFIELD_GROUP_BY_TYPE),
        default="pressure-far-field",
        help="actual Fluent type of the retained farfield physical zone",
    )
    args = parser.parse_args()
    if args.processors < 1:
        raise ValueError("processors must be positive")
    source = args.case.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    destination = work / "cold_air_sst_transport_wall_distance_configured.cas.h5"
    session = None
    result: dict[str, Any] = {
        "schema": "vortex.fluent.transport-wall-distance-configuration/v2",
        "started_utc": stamp(),
        "status": "starting",
        "source_case": str(source),
        "source_is_immutable": True,
        "processors": args.processors,
        "farfield_boundary_type": args.farfield_type,
        "actions": ["read case", "verify topology", "set wall distance", "write independent case"],
        "prohibited_actions": ["initialization", "iteration", "mesh modification", "source-case write"],
        "old_wall_distance_method": None,
        "new_wall_distance_method": None,
    }
    dump(work / "configuration.json", result)
    try:
        session = pyfluent.launch_fluent(
            precision="double", processor_count=args.processors, mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI, start_timeout=1200,
        )
        log("WALL_DISTANCE_CONFIG_SOLVER_LAUNCHED")
        st = session.settings
        st.file.read_case(file_name=str(source))
        log("WALL_DISTANCE_CONFIG_SOURCE_CASE_READ_OK")
        result["topology"] = verify_topology(st, args.farfield_type)
        result["old_wall_distance_method"] = st.mesh.wall_distance_method()
        st.mesh.wall_distance_method = "transport-eqn"
        result["new_wall_distance_method"] = st.mesh.wall_distance_method()
        if result["new_wall_distance_method"] != "transport-eqn":
            raise RuntimeError(f"wall-distance update lost: {result}")
        st.solution.initialization.initialization_type = "standard"
        result["initialization_type"] = st.solution.initialization.initialization_type()
        if result["initialization_type"] != "standard":
            raise RuntimeError(f"standard initialization setting lost: {result}")
        result["solver"] = st.setup.general.solver.get_state()
        result["air"] = st.setup.materials.fluid["air"].get_state()
        result["viscous"] = st.setup.models.viscous.get_state()
        result["farfield"] = getattr(
            st.setup.boundary_conditions, FARFIELD_GROUP_BY_TYPE[args.farfield_type]
        )["farfield"].get_state()
        result["output_case"] = str(destination)
        result["status"] = "configured"
        dump(work / "configuration.json", result)
        st.file.write_case(file_name=str(destination))
        result["status"] = "ok"
        result["finished_utc"] = stamp()
        dump(work / "configuration.json", result)
        log("WALL_DISTANCE_TRANSPORT_CONFIGURED_CASE_WRITTEN")
        log("VORTEX_WALL_DISTANCE_TRANSPORT_CONFIG_OK")
    except Exception:
        result["status"] = "failed"
        result["finished_utc"] = stamp()
        result["traceback"] = traceback.format_exc()
        dump(work / "configuration.json", result)
        raise
    finally:
        if session is not None:
            try:
                session.exit()
                log("WALL_DISTANCE_CONFIG_SOLVER_EXIT_OK")
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
