#!/usr/bin/env python3
"""在临时 Fluent 会话中验证最终 Vortex 网格的冷流配置接口。

会读取网格并临时启用 ideal-gas / energy / SST / density-based；不迭代，
不写 case/data，因此用于验证 26.1 API 和边界字段。
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
    print(f"{label}=" + json.dumps(value, default=str, ensure_ascii=True, sort_keys=True), flush=True)


def info(label: str, obj: Any) -> None:
    rec: dict[str, Any] = {"type": str(type(obj))}
    for key in ("get_active_child_names", "get_active_command_names", "get_active_query_names", "get_state"):
        try:
            rec[key] = getattr(obj, key)()
        except Exception as exc:
            rec[key + "_error"] = repr(exc)
    emit(label, rec)


def try_info(label: str, obj: Any, attr: str) -> Any | None:
    try:
        child = getattr(obj, attr)
        info(label, child)
        return child
    except Exception as exc:
        emit(label + "_ERROR", repr(exc))
        return None


def allow(label: str, leaf: Any) -> None:
    try:
        emit(label, list(leaf.allowed_values()))
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
        st = s.settings
        st.file.read_mesh(file_name=MESH)
        print("MESH_READ_OK", flush=True)

        # Required to make thermal conditions active and satisfy pressure-far-field.
        st.setup.materials.fluid["air"].density.option = "ideal-gas"
        st.setup.models.energy.enabled = True
        st.setup.models.viscous.model = "k-omega"
        st.setup.models.viscous.k_omega_model = "sst"
        st.setup.general.solver.time = "steady"
        allow("SOLVER_TYPE_ALLOWED", st.setup.general.solver.type)
        st.setup.general.solver.type = "density-based-implicit"
        print("MODELS_SET_OK", flush=True)

        info("GENERAL_SOLVER", st.setup.general.solver)
        info("AIR", st.setup.materials.fluid["air"])
        info("ENERGY", st.setup.models.energy)
        info("VISCOUS", st.setup.models.viscous)
        info("SOLUTION_METHODS", st.solution.methods)
        info("SOLUTION_CONTROLS", st.solution.controls)

        bc = st.setup.boundary_conditions
        for name in ("inlet_spray", "inlet_torch"):
            z = bc.mass_flow_inlet[name]
            print(f"MFI_ZONE={name}", flush=True)
            info(f"MFI_{name}_BEFORE", z)
            info(f"MFI_{name}_MOMENTUM_BEFORE", z.momentum)
            info(f"MFI_{name}_THERMAL_BEFORE", z.thermal)
            info(f"MFI_{name}_TURB_BEFORE", z.turbulence)
            try_info(f"MFI_{name}_THERMAL_CONDITION", z.thermal, "thermal_condition")
            # Probe direction-vector activated fields on temporary session.
            z.momentum.direction_specification = "Direction Vector"
            info(f"MFI_{name}_MOMENTUM_DIRVEC", z.momentum)
            for attr in ("coordinate_system", "flow_direction", "direction_vector", "axis_direction", "axis_origin"):
                child = try_info(f"MFI_{name}_{attr.upper()}", z.momentum, attr)
                if child is not None:
                    allow(f"MFI_{name}_{attr.upper()}_ALLOWED", child)
            # Return to the intended robust first baseline setting.
            z.momentum.direction_specification = "Normal to Boundary"

        ff = bc.pressure_far_field["farfield"]
        info("FARFIELD", ff)
        info("FARFIELD_MOMENTUM", ff.momentum)
        info("FARFIELD_THERMAL", ff.thermal)
        info("FARFIELD_TURB", ff.turbulence)
        try_info("FARFIELD_THERMAL_CONDITION", ff.thermal, "thermal_condition")
        for name in ("wall_chamber", "wall_plate"):
            z = bc.wall[name]
            info("WALL_" + name.upper(), z)
            info("WALL_" + name.upper() + "_THERMAL", z.thermal)
            try_info("WALL_" + name.upper() + "_THERMAL_CONDITION", z.thermal, "thermal_condition")

        print("COLD_SETUP_API_PROBE_OK", flush=True)
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
