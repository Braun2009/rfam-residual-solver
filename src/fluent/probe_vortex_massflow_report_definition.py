#!/usr/bin/env python3
"""只读验证 Vortex 冷流 case 的 Fluent flux report-definition 接口。

读取 case/data 后，在内存中创建质量流 report definition 并计算；绝不写回
源 case/data，也不初始化或迭代。其目的是让后续分段续算可以留下每段边界
质量通量证据，而不是仅依赖日志残差。
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

BOUNDARIES = ["inlet_spray", "inlet_torch", "farfield"]


def stamp() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def write(path: Path, document: dict[str, Any]) -> None:
    path.write_text(json.dumps(json_safe(document), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--processors", type=int, default=4)
    args = parser.parse_args()
    if args.processors < 1:
        raise ValueError("--processors must be positive")
    case, data, work = args.case.resolve(), args.data.resolve(), args.work.resolve()
    if not case.is_file() or not data.is_file():
        raise FileNotFoundError(f"missing source case/data: {case} / {data}")
    work.mkdir(parents=True, exist_ok=True)
    result: dict[str, Any] = {
        "schema": "vortex.fluent.massflow-report-definition-probe/v1",
        "started_utc": stamp(),
        "status": "starting",
        "source_case": str(case),
        "source_data": str(data),
        "source_is_read_only": True,
        "processors": args.processors,
        "boundaries": BOUNDARIES,
        "prohibited_actions": ["initialization", "iteration", "write_case", "write_data", "mesh modification"],
    }
    write(work / "result.json", result)
    session = None
    try:
        session = pyfluent.launch_fluent(
            precision="double", processor_count=args.processors, mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI, start_timeout=1200,
        )
        print("VORTEX_REPORT_PROBE_SOLVER_LAUNCHED", flush=True)
        settings = session.settings
        settings.file.read_case(file_name=str(case))
        settings.file.read_data(file_name=str(data))
        print("VORTEX_REPORT_PROBE_CASE_DATA_READ_OK", flush=True)
        reports = settings.solution.report_definitions
        name = "vortex_boundary_mass_flow_probe"
        flux = reports.flux
        before = list(flux.get_object_names())
        if name not in before:
            flux[name] = {}
        flux[name].boundaries = BOUNDARIES
        result["flux_report_names_before"] = before
        result["flux_report_state"] = flux[name].get_state()
        result["report_definition_compute"] = reports.compute(report_defs=[name])
        # Cross-check against the direct post-processing query.  Both values
        # should be available in Fluent 26.1 and use the same signed convention.
        result["direct_flux_query"] = settings.results.report.fluxes.get_mass_flow(zones=BOUNDARIES)
        result["status"] = "ok"
        result["finished_utc"] = stamp()
        write(work / "result.json", result)
        print("VORTEX_MASSFLOW_REPORT_DEFINITION_PROBE_OK", flush=True)
    except Exception:
        result["status"] = "failed"
        result["finished_utc"] = stamp()
        result["traceback"] = traceback.format_exc()
        write(work / "result.json", result)
        raise
    finally:
        if session is not None:
            try:
                session.exit()
                print("VORTEX_REPORT_PROBE_SOLVER_EXIT_OK", flush=True)
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
