#!/usr/bin/env python3
"""Read-only post-processing audit for a completed Vortex cold-air Fluent restart.

This script never reads the production mesh and never changes the supplied
checkpoint.  It opens an independent Solver session, reads that checkpoint,
executes mesh/check plus report queries, and writes JSON evidence in --work.

The scope remains explicitly limited to ideal-air, density-based SST
commissioning.  Its reported fluxes and extrema are not LOX/C3H8, GERG, EDC,
or reacting performance predictions.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import sys
import traceback
from typing import Any, Callable

os.environ.setdefault("AWP_ROOT261", "/ansys_inc/v261")
os.environ.setdefault("ANSYSLMD_LICENSE_FILE", "1055@localhost")
os.environ.setdefault("ANSYSCL_FNP_PATH", os.environ["ANSYSLMD_LICENSE_FILE"])

import ansys.fluent.core as pyfluent
from ansys.fluent.core import FluentMode, UIMode

BOUNDARIES = ["inlet_spray", "inlet_torch", "farfield"]
WALLS = ["wall_chamber", "wall_plate"]
FLUID_ZONES = [
    "vortex_full360_v2_zoned_ascii-prism-cells",
    "vortex_full360_v2_zoned_ascii-tet-cells",
]
BASE_EXPECTED = {
    "mass_flow_inlet": {"inlet_spray", "inlet_torch"},
    "wall": {"wall_chamber", "wall_plate"},
}
FARFIELD_GROUP_BY_TYPE = {
    "pressure-far-field": "pressure_far_field",
    "pressure-outlet": "pressure_outlet",
}


def utc() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def jsonify(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): jsonify(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonify(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def write(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(jsonify(payload), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def scalar_or_dict(value: Any) -> Any:
    """Preserve Fluent direct-query dicts and report-definition list payloads."""
    if isinstance(value, dict):
        return {str(k): scalar_or_dict(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [scalar_or_dict(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    try:
        return float(value)
    except (TypeError, ValueError):
        return repr(value)


def extract_net(value: Any) -> float | None:
    """Extract a scalar from direct flux dicts or Fluent 26.1 report-def results."""
    if isinstance(value, dict):
        for key in ("Net", "net"):
            if key in value:
                return extract_net(value[key])
        if len(value) == 1:
            return extract_net(next(iter(value.values())))
        return None
    if isinstance(value, (list, tuple)):
        # reports.compute() returns e.g. [{name: [0.009962, "kg/s"]}].
        if len(value) == 1:
            return extract_net(value[0])
        if value:
            try:
                return float(value[0])
            except (TypeError, ValueError):
                return None
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def query(label: str, func: Callable[[], Any], audit: dict[str, Any]) -> Any:
    try:
        result = func()
        audit["queries"][label] = {"ok": True, "value": scalar_or_dict(result)}
        return result
    except Exception as exc:  # audit should report version/API variation rather than discard all evidence
        audit["queries"][label] = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
        return None


def query_first_available_volume_field(
    label: str,
    report: Any,
    operation: str,
    candidates: tuple[str, ...],
    audit: dict[str, Any],
) -> Any:
    """Query a version-dependent Fluent cell-function name without hiding failures."""
    errors: list[str] = []
    method = getattr(report.volume_integrals, operation)
    for field in candidates:
        try:
            result = method(cell_function=field, cell_zones=FLUID_ZONES)
            audit["queries"][label] = {
                "ok": True,
                "cell_function": field,
                "value": scalar_or_dict(result),
            }
            return result
        except Exception as exc:
            errors.append(f"{field}: {type(exc).__name__}: {exc}")
    audit["queries"][label] = {"ok": False, "errors": errors}
    return None


def assert_topology(settings: Any, farfield_type: str) -> dict[str, list[str]]:
    bc = settings.setup.boundary_conditions
    farfield_group = FARFIELD_GROUP_BY_TYPE[farfield_type]
    expected = {
        **BASE_EXPECTED,
        "pressure_far_field": {"farfield"} if farfield_group == "pressure_far_field" else set(),
        "pressure_outlet": {"farfield"} if farfield_group == "pressure_outlet" else set(),
    }
    got: dict[str, list[str]] = {}
    for group, required_names in expected.items():
        names = set(getattr(bc, group).get_object_names())
        got[group] = sorted(names)
        if names != required_names:
            raise RuntimeError(f"{group}: expected {sorted(required_names)}, got {sorted(names)}")
    fluid = set(settings.setup.cell_zone_conditions.fluid.get_object_names())
    got["fluid"] = sorted(fluid)
    if fluid != set(FLUID_ZONES):
        raise RuntimeError(f"fluid zones: expected {FLUID_ZONES}, got {sorted(fluid)}")
    return got


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--processors", type=int, default=4)
    parser.add_argument(
        "--farfield-type", choices=tuple(FARFIELD_GROUP_BY_TYPE),
        default="pressure-far-field",
        help="actual Fluent type of the retained farfield physical zone",
    )
    args = parser.parse_args()
    if args.processors < 1:
        raise ValueError("--processors must be positive")
    case, data, work = args.case.resolve(), args.data.resolve(), args.work.resolve()
    if not case.is_file() or not data.is_file():
        raise FileNotFoundError(f"missing input: {case} / {data}")
    work.mkdir(parents=True, exist_ok=True)
    audit: dict[str, Any] = {
        "schema": "vortex.fluent.cold-air-postprocess-audit/v2",
        "started_utc": utc(),
        "status": "starting",
        "input": {"case": str(case), "data": str(data)},
        "processors": args.processors,
        "farfield_boundary_type": args.farfield_type,
        "scope": "ideal-air density-based implicit SST cold-flow commissioning only; no real-fluid/reaction/EDC claim",
        "queries": {},
    }
    write(work / "audit.json", audit)
    session = None
    try:
        session = pyfluent.launch_fluent(
            precision="double", processor_count=args.processors, mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI, start_timeout=1200,
        )
        print("VORTEX_AUDIT_SOLVER_LAUNCHED", flush=True)
        settings = session.settings
        settings.file.read_case(file_name=str(case))
        settings.file.read_data(file_name=str(data))
        print("VORTEX_AUDIT_CASE_DATA_READ_OK", flush=True)
        audit["topology"] = assert_topology(settings, args.farfield_type)
        setup = settings.setup
        farfield_group = getattr(
            setup.boundary_conditions, FARFIELD_GROUP_BY_TYPE[args.farfield_type]
        )
        audit["readback"] = {
            "solver": setup.general.solver.get_state(),
            "energy": setup.models.energy.get_state(),
            "viscous": setup.models.viscous.get_state(),
            "air": setup.materials.fluid["air"].get_state(),
            "inlets": {name: setup.boundary_conditions.mass_flow_inlet[name].get_state() for name in ("inlet_spray", "inlet_torch")},
            "farfield": farfield_group["farfield"].get_state(),
            "walls": {name: setup.boundary_conditions.wall[name].get_state() for name in WALLS},
            "initialization": settings.solution.initialization.get_state(),
            "solution_methods": settings.solution.methods.get_state(),
            "solution_controls": settings.solution.controls.get_state(),
            "wall_distance_method": settings.mesh.wall_distance_method(),
        }
        session.mesh.check()
        print("VORTEX_AUDIT_MESH_CHECK_OK", flush=True)

        report = settings.results.report
        # Create a collision-free flux report definition only in this independent
        # in-memory session, then cross-check it against the direct report query.
        # It is never written to the supplied source case/data.
        report_definitions = settings.solution.report_definitions
        flux_definitions = report_definitions.flux
        existing_flux_definitions = list(flux_definitions.get_object_names())
        report_definition_name = "vortex_audit_boundary_mass_flow"
        suffix = 1
        while report_definition_name in existing_flux_definitions:
            suffix += 1
            report_definition_name = f"vortex_audit_boundary_mass_flow_{suffix}"
        flux_definitions[report_definition_name] = {}
        flux_definitions[report_definition_name].boundaries = BOUNDARIES
        audit["mass_flow_report_definition"] = {
            "name": report_definition_name,
            "existing_flux_definition_names": existing_flux_definitions,
            "state": flux_definitions[report_definition_name].get_state(),
            "persistence": "created only in the independent in-memory audit session; source case/data are not written",
        }
        report_definition_flux = query(
            "boundary_mass_flow_report_definition",
            lambda: report_definitions.compute(report_defs=[report_definition_name]),
            audit,
        )
        # Signed values are intentionally retained.  In this Fluent case the
        # prescribed inlets are positive and a net-discharging farfield is negative.
        all_flux = query("boundary_mass_flow_all", lambda: report.fluxes.get_mass_flow(zones=BOUNDARIES), audit)
        spray_flux = query("boundary_mass_flow_inlet_spray", lambda: report.fluxes.get_mass_flow(zones=["inlet_spray"]), audit)
        torch_flux = query("boundary_mass_flow_inlet_torch", lambda: report.fluxes.get_mass_flow(zones=["inlet_torch"]), audit)
        far_flux = query("boundary_mass_flow_farfield", lambda: report.fluxes.get_mass_flow(zones=["farfield"]), audit)
        query("wall_heat_transfer", lambda: report.fluxes.get_heat_transfer(zones=WALLS), audit)
        for field, method in (("temperature", "temperature"), ("density", "density"), ("absolute-pressure", "absolute_pressure"), ("mach-number", "mach_number")):
            query(f"domain_{method}_minimum", lambda f=field: report.volume_integrals.get_minimum(cell_function=f, cell_zones=FLUID_ZONES), audit)
            query(f"domain_{method}_maximum", lambda f=field: report.volume_integrals.get_maximum(cell_function=f, cell_zones=FLUID_ZONES), audit)
        turbulent_viscosity_ratio_min = query_first_available_volume_field(
            "domain_turbulent_viscosity_ratio_minimum", report, "get_minimum",
            ("viscosity-ratio", "turb-viscosity-ratio", "turbulent-viscosity-ratio"), audit,
        )
        turbulent_viscosity_ratio_max = query_first_available_volume_field(
            "domain_turbulent_viscosity_ratio_maximum", report, "get_maximum",
            ("viscosity-ratio", "turb-viscosity-ratio", "turbulent-viscosity-ratio"), audit,
        )
        turbulent_viscosity_ratio_max_net = extract_net(turbulent_viscosity_ratio_max)
        audit["turbulence_limiter"] = {
            "fluent_viscosity_ratio_limit": 1.0e5,
            "minimum": extract_net(turbulent_viscosity_ratio_min),
            "maximum": turbulent_viscosity_ratio_max_net,
            "limit_is_reached_or_exceeded": (
                turbulent_viscosity_ratio_max_net is not None
                and turbulent_viscosity_ratio_max_net >= 1.0e5 * (1.0 - 1.0e-12)
            ),
        }
        query("farfield_temperature_mass_avg", lambda: report.surface_integrals.get_mass_weighted_avg(report_of="temperature", surface_names=["farfield"]), audit)
        query("wall_absolute_pressure_area_avg", lambda: report.surface_integrals.get_area_weighted_avg(report_of="absolute-pressure", surface_names=WALLS), audit)

        report_definition_net = extract_net(report_definition_flux)
        direct_all_flux_net = extract_net(all_flux)
        audit["mass_flow_report_definition"].update({
            "computed_net_kg_s": report_definition_net,
            "direct_all_boundaries_net_kg_s": direct_all_flux_net,
            "absolute_difference_to_direct_kg_s": (
                abs(report_definition_net - direct_all_flux_net)
                if report_definition_net is not None and direct_all_flux_net is not None
                else None
            ),
        })

        individual = [extract_net(item) for item in (spray_flux, torch_flux, far_flux)]
        if all(v is not None and math.isfinite(v) for v in individual):
            signed_sum = float(sum(v for v in individual if v is not None))
            inlet_magnitude = abs(float(individual[0])) + abs(float(individual[1]))
            # For these mass-flow inlets Fluent reports positive values while a
            # net discharging farfield reports negative.  A relative imbalance
            # near zero therefore indicates closure; a value near one means the
            # prescribed inlet mass is still accumulating in the domain.
            audit["mass_balance"] = {
                "sign_convention": "Fluent signed boundary fluxes retained; prescribed inlets are positive and net discharging farfield is negative in this case",
                "inlet_spray": individual[0],
                "inlet_torch": individual[1],
                "farfield": individual[2],
                "signed_sum_kg_s": signed_sum,
                "total_prescribed_inlet_magnitude_kg_s": inlet_magnitude,
                "relative_mass_imbalance_to_inlets": abs(signed_sum) / inlet_magnitude if inlet_magnitude else None,
                "farfield_discharge_fraction_of_inlets": (-float(individual[2]) / inlet_magnitude) if inlet_magnitude else None,
                "interpretation": "relative_mass_imbalance_to_inlets=0 is closed; =1 is full inlet accumulation",
            }
        else:
            audit["mass_balance"] = {"available": False, "reason": "at least one flux query failed or returned nonnumeric data"}
        if all_flux is not None:
            audit["mass_balance"]["combined_query_net_kg_s"] = extract_net(all_flux)

        audit["status"] = "ok"
        audit["finished_utc"] = utc()
        write(work / "audit.json", audit)
        print("VORTEX_3D_COLD_AIR_SST_AUDIT_OK", flush=True)
    except Exception:
        audit["status"] = "failed"
        audit["finished_utc"] = utc()
        audit["traceback"] = traceback.format_exc()
        write(work / "audit.json", audit)
        raise
    finally:
        if session is not None:
            try:
                session.exit()
                print("VORTEX_AUDIT_SOLVER_EXIT_OK", flush=True)
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
