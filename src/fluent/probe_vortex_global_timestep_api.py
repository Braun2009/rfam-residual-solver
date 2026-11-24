#!/usr/bin/env python3
"""In-memory-only Fluent 26.1 probe for density-based global pseudo time step.

Reads an immutable Vortex pressure-outlet configuration in an independent
Solver session, changes only the in-memory pseudo-time formulation to
``global-time-step``, and records API states/active children.  It does not
initialize, iterate, modify the mesh, or write case/data, so the source case
remains unchanged.
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
        "schema": "vortex.fluent.global-time-step-api-probe/v1",
        "started_utc": stamp(),
        "status": "starting",
        "source_case": str(case),
        "source_is_read_only": True,
        "processors": args.processors,
        "in_memory_action": "set solution.methods.pseudo_time_method.formulation.density_based_solver=global-time-step",
        "prohibited_actions": [
            "initialization", "iteration", "write_case", "write_data", "mesh modification"
        ],
    }
    dump(work / "result.json", result)
    session = None
    try:
        session = pyfluent.launch_fluent(
            precision="double", processor_count=args.processors, mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI, start_timeout=1200,
        )
        print("VORTEX_GLOBAL_TIMESTEP_PROBE_SOLVER_LAUNCHED", flush=True)
        settings = session.settings
        settings.file.read_case(file_name=str(case))
        print("VORTEX_GLOBAL_TIMESTEP_PROBE_CASE_READ_OK", flush=True)
        method = settings.solution.methods.pseudo_time_method
        formulation = method.formulation
        density_based_leaf = formulation.density_based_solver
        result["before"] = {
            "pseudo_time_method": describe(method),
            "formulation": describe(formulation),
            "density_based_solver": describe(density_based_leaf),
            "solution_controls": describe(settings.solution.controls),
        }
        formulation.density_based_solver = "global-time-step"
        if formulation.density_based_solver() != "global-time-step":
            raise RuntimeError(
                f"global-time-step assignment lost: {formulation.density_based_solver()}"
            )
        result["after"] = {
            "pseudo_time_method": describe(method),
            "formulation": describe(formulation),
            "density_based_solver": describe(formulation.density_based_solver),
            "solution_controls": describe(settings.solution.controls),
            "solution_methods": describe(settings.solution.methods),
        }
        result["status"] = "ok"
        result["finished_utc"] = stamp()
        dump(work / "result.json", result)
        print("VORTEX_GLOBAL_TIMESTEP_API_PROBE_OK", flush=True)
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
                print("VORTEX_GLOBAL_TIMESTEP_PROBE_SOLVER_EXIT_OK", flush=True)
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
