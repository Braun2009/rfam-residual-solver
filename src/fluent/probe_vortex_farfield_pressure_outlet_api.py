#!/usr/bin/env python3
"""只读探测 Vortex 配置 case 中 farfield -> pressure-outlet 的 Fluent API。

在临时 Solver 会话中读取 immutable r2 configured case，将 farfield 的类型
仅在内存中从 pressure-far-field 改成 pressure-outlet，然后枚举可设置字段。
不会初始化、迭代或写出 case/data。该探测用于未来非理想气体远场策略，不是
生产物理计算。
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


def safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def describe(obj: Any) -> dict[str, Any]:
    result: dict[str, Any] = {"type": str(type(obj))}
    for name in ("get_active_child_names", "get_active_command_names", "get_active_query_names", "get_state"):
        try:
            result[name] = getattr(obj, name)()
        except Exception as exc:
            result[name + "_error"] = f"{type(exc).__name__}: {exc}"
    return result


def dump(path: Path, document: dict[str, Any]) -> None:
    path.write_text(json.dumps(safe(document), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--processors", type=int, default=4)
    args = parser.parse_args()
    if args.processors < 1:
        raise ValueError("--processors must be positive")
    if not SOURCE.is_file():
        raise FileNotFoundError(SOURCE)
    work = args.work.resolve(); work.mkdir(parents=True, exist_ok=True)
    result: dict[str, Any] = {
        "schema": "vortex.fluent.farfield-pressure-outlet-api-probe/v1",
        "source_case": str(SOURCE),
        "source_is_read_only": True,
        "processors": args.processors,
        "prohibited_actions": ["initialization", "iteration", "write_case", "write_data", "mesh modification on disk"],
        "status": "starting",
    }
    dump(work / "result.json", result)
    session = None
    try:
        session = pyfluent.launch_fluent(
            precision="double", processor_count=args.processors, mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI, start_timeout=1200,
        )
        print("VORTEX_PRESSURE_OUTLET_PROBE_SOLVER_LAUNCHED", flush=True)
        settings = session.settings
        settings.file.read_case(file_name=str(SOURCE))
        print("VORTEX_PRESSURE_OUTLET_PROBE_CASE_READ_OK", flush=True)
        bc = settings.setup.boundary_conditions
        result["before"] = {
            "pressure_far_field": list(bc.pressure_far_field.get_object_names()),
            "pressure_outlet": list(bc.pressure_outlet.get_object_names()),
        }
        settings.mesh.modify_zones.zone_type(
            zone_names=["farfield"], new_type="pressure-outlet"
        )
        result["after"] = {
            "pressure_far_field": list(bc.pressure_far_field.get_object_names()),
            "pressure_outlet": list(bc.pressure_outlet.get_object_names()),
        }
        if "farfield" not in result["after"]["pressure_outlet"]:
            raise RuntimeError(f"farfield conversion failed: {result['after']}")
        outlet = bc.pressure_outlet["farfield"]
        result["pressure_outlet"] = describe(outlet)
        for child_name in ("momentum", "thermal", "turbulence", "species", "multiphase"):
            try:
                result["pressure_outlet"][child_name] = describe(getattr(outlet, child_name))
            except Exception as exc:
                result["pressure_outlet"][child_name + "_error"] = f"{type(exc).__name__}: {exc}"
        result["status"] = "ok"
        dump(work / "result.json", result)
        print("VORTEX_FARFIELD_PRESSURE_OUTLET_API_PROBE_OK", flush=True)
    except Exception:
        result["status"] = "failed"
        result["traceback"] = traceback.format_exc()
        dump(work / "result.json", result)
        raise
    finally:
        if session is not None:
            try:
                session.exit()
                print("VORTEX_PRESSURE_OUTLET_PROBE_SOLVER_EXIT_OK", flush=True)
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
