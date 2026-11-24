#!/usr/bin/env python3
"""只读核验 Fluent Solver 对最终 Vortex 3D 网格的导入。

本脚本不写 case/data，也不修改网格。它用于确认 HDF5 ``units=mm``
在 Solver 侧的缩放、体/面区名称与类型，以及 SST/能量模型 API 的可用性。
"""
from __future__ import annotations

import json
import os
import sys
import traceback

os.environ.setdefault("AWP_ROOT261", "/ansys_inc/v261")
os.environ.setdefault("ANSYSLMD_LICENSE_FILE", "1055@localhost")
os.environ.setdefault("ANSYSCL_FNP_PATH", os.environ["ANSYSLMD_LICENSE_FILE"])

import ansys.fluent.core as pyfluent
from ansys.fluent.core import FluentMode, UIMode

MESH = "/root/vortex_fluent/gmsh_fixed_v2/strict_tui_prism_final.msh.h5"


def emit(label: str, value: object) -> None:
    try:
        rendered = json.dumps(value, ensure_ascii=True, default=str, sort_keys=True)
    except Exception:
        rendered = repr(value)
    print(f"{label}={rendered}", flush=True)


def inspect_object(label: str, obj: object) -> None:
    record: dict[str, object] = {"type": str(type(obj))}
    for method_name in (
        "get_state",
        "get_active_child_names",
        "get_active_command_names",
        "get_active_query_names",
    ):
        method = getattr(obj, method_name, None)
        if method is None:
            continue
        try:
            record[method_name] = method()
        except Exception as exc:  # API capability inventory only
            record[f"{method_name}_error"] = repr(exc)
    emit(label, record)


def boundary_inventory(settings: object) -> None:
    bc = settings.setup.boundary_conditions
    groups = (
        "mass_flow_inlet",
        "pressure_far_field",
        "pressure_outlet",
        "wall",
        "interior",
        "symmetry",
        "velocity_inlet",
    )
    result: dict[str, object] = {}
    for group_name in groups:
        try:
            group = getattr(bc, group_name)
            result[group_name] = list(group.get_object_names())
        except Exception as exc:
            result[group_name] = {"error": repr(exc)}
    emit("BOUNDARY_INVENTORY", result)


def main() -> None:
    session = None
    try:
        session = pyfluent.launch_fluent(
            precision="double",
            processor_count=4,
            mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI,
            start_timeout=900,
        )
        print("SOLVER_LAUNCHED", flush=True)
        session.settings.file.read_mesh(file_name=MESH)
        print("FINAL_MESH_READ_OK", flush=True)

        settings = session.settings
        inspect_object("SETTINGS_MESH", settings.mesh)
        for attribute in ("units", "scale", "length_unit", "mesh_units"):
            try:
                inspect_object(f"SETTINGS_MESH_{attribute.upper()}", getattr(settings.mesh, attribute))
            except Exception as exc:
                emit(f"SETTINGS_MESH_{attribute.upper()}_UNAVAILABLE", repr(exc))

        boundary_inventory(settings)
        cell_zones = settings.setup.cell_zone_conditions
        emit("CELL_ZONE_NAMES", list(cell_zones.get_object_names()))
        for group_name in ("fluid", "solid"):
            try:
                emit(f"CELL_ZONE_GROUP_{group_name.upper()}", list(getattr(cell_zones, group_name).get_object_names()))
            except Exception as exc:
                emit(f"CELL_ZONE_GROUP_{group_name.upper()}_ERROR", repr(exc))

        inspect_object("ENERGY_MODEL", settings.setup.models.energy)
        inspect_object("VISCOUS_MODEL", settings.setup.models.viscous)
        emit("DIMENSION", getattr(session, "dimension", "unavailable"))
        print("FINAL_SOLVER_PROBE_OK", flush=True)
    finally:
        if session is not None:
            session.exit()
            print("SOLVER_EXIT_OK", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
