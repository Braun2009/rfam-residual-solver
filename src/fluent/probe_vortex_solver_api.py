#!/usr/bin/env python3
"""枚举 Fluent 26.1 中最终 Vortex 网格的 Solver 设置接口。

仅 read_mesh；不改模型、不写 case/data、不迭代。输出用于让后续冷流脚本使用
当前版本的真实 PyFluent 字段，而不是沿用旧 2D 脚本的猜测路径。
"""
from __future__ import annotations

import json
import os
import sys
import traceback
from typing import Any

os.environ.setdefault("AWP_ROOT261", "/ansys_inc/v261")
os.environ.setdefault("ANSYSLMD_LICENSE_FILE", "1055@localhost")
os.environ.setdefault("ANSYSCL_FNP_PATH", os.environ["ANSYSLMD_LICENSE_FILE"])

import ansys.fluent.core as pyfluent
from ansys.fluent.core import FluentMode, UIMode

MESH = "/root/vortex_fluent/gmsh_fixed_v2/strict_tui_prism_final.msh.h5"


def emit(label: str, value: Any) -> None:
    print(label + "=" + json.dumps(value, default=str, ensure_ascii=True, sort_keys=True), flush=True)


def describe(label: str, obj: Any, *, state: bool = True) -> None:
    record: dict[str, Any] = {"type": str(type(obj))}
    for name in ("get_active_child_names", "get_active_command_names", "get_active_query_names"):
        try:
            record[name] = getattr(obj, name)()
        except Exception as exc:
            record[name + "_error"] = repr(exc)
    if state:
        try:
            record["state"] = obj.get_state()
        except Exception as exc:
            record["state_error"] = repr(exc)
    emit(label, record)


def optional(label: str, parent: Any, name: str, *, state: bool = True) -> Any | None:
    try:
        obj = getattr(parent, name)
    except Exception as exc:
        emit(label + "_UNAVAILABLE", repr(exc))
        return None
    describe(label, obj, state=state)
    return obj


def allowed(label: str, obj: Any) -> None:
    try:
        emit(label, list(obj.allowed_values()))
    except Exception as exc:
        emit(label + "_ERROR", repr(exc))


def main() -> None:
    s = None
    try:
        s = pyfluent.launch_fluent(
            precision="double", processor_count=4, mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI, start_timeout=900,
        )
        print("SOLVER_LAUNCHED", flush=True)
        s.settings.file.read_mesh(file_name=MESH)
        print("MESH_READ_OK", flush=True)
        st = s.settings

        setup = optional("SETUP", st, "setup", state=False)
        solution = optional("SOLUTION", st, "solution", state=False)
        if setup is not None:
            general = optional("SETUP_GENERAL", setup, "general")
            if general is not None:
                for n in ("solver", "solver_type", "time", "gravity", "operating_conditions"):
                    child = optional("SETUP_GENERAL_" + n.upper(), general, n)
                    if child is not None:
                        allowed("SETUP_GENERAL_" + n.upper() + "_ALLOWED", child)
            models = optional("SETUP_MODELS", setup, "models", state=False)
            if models is not None:
                energy = optional("MODEL_ENERGY", models, "energy")
                viscous = optional("MODEL_VISCOUS", models, "viscous")
                if viscous is not None:
                    optional("MODEL_VISCOUS_KOMEGA", viscous, "k_omega")
                    optional("MODEL_VISCOUS_OPTIONS", viscous, "options")
                    allowed("MODEL_VISCOUS_MODEL_ALLOWED", viscous.model)
                    allowed("MODEL_VISCOUS_KOMEGA_ALLOWED", viscous.k_omega_model)
                if energy is not None:
                    optional("MODEL_ENERGY_ENABLED", energy, "enabled")

            mats = optional("MATERIALS", setup, "materials", state=False)
            if mats is not None:
                for n in ("fluid", "mixture", "solid"):
                    child = optional("MATERIALS_" + n.upper(), mats, n, state=False)
                    if child is not None:
                        try:
                            emit("MATERIALS_" + n.upper() + "_NAMES", list(child.get_object_names()))
                        except Exception as exc:
                            emit("MATERIALS_" + n.upper() + "_NAMES_ERROR", repr(exc))
                try:
                    air = mats.fluid["air"]
                    describe("MATERIAL_AIR", air)
                    describe("MATERIAL_AIR_DENSITY", air.density)
                    allowed("MATERIAL_AIR_DENSITY_OPTIONS", air.density.option)
                except Exception as exc:
                    emit("MATERIAL_AIR_ERROR", repr(exc))

            bc = optional("BOUNDARY_CONDITIONS", setup, "boundary_conditions", state=False)
            if bc is not None:
                for group_name, zone_names in {
                    "mass_flow_inlet": ("inlet_spray", "inlet_torch"),
                    "pressure_far_field": ("farfield",),
                    "wall": ("wall_chamber", "wall_plate"),
                }.items():
                    group = optional("BC_GROUP_" + group_name.upper(), bc, group_name, state=False)
                    if group is None:
                        continue
                    try:
                        emit("BC_GROUP_" + group_name.upper() + "_NAMES", list(group.get_object_names()))
                    except Exception as exc:
                        emit("BC_GROUP_" + group_name.upper() + "_NAMES_ERROR", repr(exc))
                    for zone_name in zone_names:
                        try:
                            zone = group[zone_name]
                            describe("BC_" + zone_name.upper(), zone, state=False)
                            for child_name in ("momentum", "thermal", "turbulence", "geometry"):
                                child = optional("BC_" + zone_name.upper() + "_" + child_name.upper(), zone, child_name)
                                if child is not None:
                                    # Relevant enum leaves, if this boundary exposes them.
                                    for leaf in ("mass_flow_specification", "direction_specification", "direction_specification_method", "coordinate_system", "thermal_condition", "turbulence_specification"):
                                        try:
                                            allowed("BC_" + zone_name.upper() + "_" + child_name.upper() + "_" + leaf.upper() + "_ALLOWED", getattr(child, leaf))
                                        except Exception:
                                            pass
                        except Exception as exc:
                            emit("BC_" + zone_name.upper() + "_ERROR", repr(exc))

            czc = optional("CELL_ZONE_CONDITIONS", setup, "cell_zone_conditions", state=False)
            if czc is not None:
                try:
                    emit("CELL_ZONE_NAMES", list(czc))
                except Exception as exc:
                    emit("CELL_ZONE_NAMES_ERROR", repr(exc))
                for group_name in ("fluid", "solid"):
                    group = optional("CELL_ZONE_" + group_name.upper(), czc, group_name, state=False)
                    if group is not None:
                        try:
                            emit("CELL_ZONE_" + group_name.upper() + "_NAMES", list(group.get_object_names()))
                        except Exception as exc:
                            emit("CELL_ZONE_" + group_name.upper() + "_NAMES_ERROR", repr(exc))
                        for zone_name in list(group.get_object_names()) if hasattr(group, "get_object_names") else []:
                            try:
                                describe("CELL_ZONE_" + zone_name.upper(), group[zone_name])
                            except Exception as exc:
                                emit("CELL_ZONE_" + zone_name.upper() + "_ERROR", repr(exc))

        if solution is not None:
            for n in ("methods", "controls", "initialization", "run_calculation", "report_definitions", "monitors"):
                obj = optional("SOLUTION_" + n.upper(), solution, n, state=False)
                if obj is not None:
                    describe("SOLUTION_" + n.upper() + "_DETAIL", obj)

        print("SOLVER_API_PROBE_OK", flush=True)
    finally:
        if s is not None:
            s.exit()
            print("SOLVER_EXIT_OK", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
