#!/usr/bin/env python3
"""Read-only Fluent 26.1 API probe for Vortex low-Mach initialization controls.

The current open ``pressure-outlet`` cold-air branch is numerically executable
but has not established mass closure.  This probe reads an explicitly supplied
case in an independent Solver session and records the available initialization,
density-based pseudo-time, solver-method, and control interfaces.  It never
assigns a setting, initializes, iterates, modifies the mesh, or writes case/data.
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

DEFAULT_CASE = Path(
    "/root/vortex_fluent/gmsh_fixed_v2/solver_cold_air_sst_pressure_outlet_config_r34/"
    "cold_air_sst_pressure_outlet_configured.cas.h5"
)


def stamp() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def dump(path: Path, document: dict[str, Any]) -> None:
    path.write_text(json.dumps(safe(document), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def describe(obj: Any) -> dict[str, Any]:
    """Collect introspection data only; do not invoke modifying commands."""
    result: dict[str, Any] = {"type": str(type(obj))}
    for name in (
        "get_active_child_names",
        "get_active_command_names",
        "get_active_query_names",
        "get_state",
        "allowed_values",
    ):
        try:
            result[name] = safe(getattr(obj, name)())
        except Exception as exc:
            result[name + "_error"] = f"{type(exc).__name__}: {exc}"
    return result


def resolve(root: Any, dotted_path: str) -> Any:
    current = root
    for component in dotted_path.split("."):
        current = getattr(current, component)
    return current


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, default=DEFAULT_CASE)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--processors", type=int, default=4)
    args = parser.parse_args()
    if args.processors < 1:
        raise ValueError("--processors must be positive")
    case, work = args.case.resolve(), args.work.resolve()
    if not case.is_file():
        raise FileNotFoundError(case)
    work.mkdir(parents=True, exist_ok=True)

    result: dict[str, Any] = {
        "schema": "vortex.fluent.low-mach-controls-api-probe/v1",
        "started_utc": stamp(),
        "status": "starting",
        "source_case": str(case),
        "source_is_read_only": True,
        "processors": args.processors,
        "prohibited_actions": [
            "setting assignment",
            "initialization",
            "iteration",
            "write_case",
            "write_data",
            "mesh modification",
        ],
    }
    dump(work / "result.json", result)
    session = None
    try:
        session = pyfluent.launch_fluent(
            precision="double", processor_count=args.processors, mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI, start_timeout=1200,
        )
        print("VORTEX_LOW_MACH_API_PROBE_SOLVER_LAUNCHED", flush=True)
        settings = session.settings
        settings.file.read_case(file_name=str(case))
        print("VORTEX_LOW_MACH_API_PROBE_CASE_READ_OK", flush=True)

        targets = (
            "setup.general.solver",
            "setup.general.solver.type",
            "setup.general.solver.time",
            "solution.initialization",
            "solution.initialization.initialization_type",
            "solution.initialization.defaults",
            "solution.run_calculation",
            "solution.methods",
            "solution.methods.pseudo_time_method",
            "solution.methods.pseudo_time_method.formulation",
            "solution.methods.pseudo_time_method.formulation.density_based_solver",
            "solution.methods.pseudo_time_method.formulation.pressure_based_solver",
            "solution.controls",
            "solution.controls.courant_number",
            "solution.controls.advanced",
            "mesh",
            "mesh.wall_distance_method",
        )
        result["interfaces"] = {}
        for path in targets:
            try:
                result["interfaces"][path] = describe(resolve(settings, path))
            except Exception as exc:
                result["interfaces"][path] = {
                    "error": f"{type(exc).__name__}: {exc}"
                }
        result["status"] = "ok"
        result["finished_utc"] = stamp()
        dump(work / "result.json", result)
        print("VORTEX_LOW_MACH_CONTROLS_API_PROBE_OK", flush=True)
    except Exception:
        result["status"] = "failed"
        result["finished_utc"] = stamp()
        result["traceback"] = traceback.format_exc()
        dump(work / "result.json", result)
        raise
    finally:
        if session is not None:
            try:
                session.exit()
                print("VORTEX_LOW_MACH_API_PROBE_SOLVER_EXIT_OK", flush=True)
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
