#!/usr/bin/env python3
"""Create an immutable-source Vortex case with explicit hybrid-init options.

This utility changes only the hybrid-initialization settings stored in a new
case.  It is intended for isolated initialization diagnostics after the
pressure-outlet branch showed markedly better initial mass balance with Fluent
hybrid initialization than with standard initialization.

The source case, source data, and mesh are never modified.  This utility does
not initialize or iterate; a later independently tracked Solver run owns that
operation.  The scope remains ideal-air, density-based SST commissioning only.
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


def emit(message: str) -> None:
    print(message, flush=True)


def dump(path: Path, document: dict[str, Any]) -> None:
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def verify_topology_models(settings: Any, farfield_type: str) -> dict[str, Any]:
    setup = settings.setup
    bc = setup.boundary_conditions
    farfield_group = FARFIELD_GROUP_BY_TYPE[farfield_type]
    expected_boundaries = {
        **BASE_EXPECTED,
        "pressure_far_field": {"farfield"} if farfield_group == "pressure_far_field" else set(),
        "pressure_outlet": {"farfield"} if farfield_group == "pressure_outlet" else set(),
    }
    result: dict[str, Any] = {}
    for group, expected_names in expected_boundaries.items():
        actual_names = set(getattr(bc, group).get_object_names())
        result[group] = sorted(actual_names)
        if actual_names != expected_names:
            raise RuntimeError(f"{group}: expected {sorted(expected_names)}, got {sorted(actual_names)}")
    fluid = set(setup.cell_zone_conditions.fluid.get_object_names())
    result["fluid"] = sorted(fluid)
    if fluid != EXPECTED_FLUID:
        raise RuntimeError(f"fluid zones: expected {sorted(EXPECTED_FLUID)}, got {sorted(fluid)}")

    result["solver"] = setup.general.solver.get_state()
    result["energy"] = setup.models.energy.get_state()
    result["viscous"] = setup.models.viscous.get_state()
    result["air"] = setup.materials.fluid["air"].get_state()
    result["farfield"] = getattr(bc, farfield_group)["farfield"].get_state()
    result["wall_distance_method"] = settings.mesh.wall_distance_method()
    result["solution_methods"] = settings.solution.methods.get_state()
    result["solution_controls"] = settings.solution.controls.get_state()
    if result["solver"]["type"] != "density-based-implicit":
        raise RuntimeError(f"solver mismatch: {result['solver']}")
    if result["solver"]["time"] != "steady":
        raise RuntimeError(f"time formulation mismatch: {result['solver']}")
    if not result["energy"]["enabled"]:
        raise RuntimeError("energy is disabled")
    if result["air"]["density"]["option"] != "ideal-gas":
        raise RuntimeError(f"air density mismatch: {result['air']}")
    if result["viscous"]["model"] != "k-omega" or result["viscous"]["k_omega_model"] != "sst":
        raise RuntimeError(f"SST mismatch: {result['viscous']}")
    return result


def configure_hybrid_options(settings: Any, iterations: int) -> dict[str, Any]:
    if iterations < 1:
        raise ValueError("--hybrid-iterations must be positive")
    initialization = settings.solution.initialization
    before_type = initialization.initialization_type()
    before_state = initialization.get_state()
    initialization.initialization_type = "hybrid"
    if initialization.initialization_type() != "hybrid":
        raise RuntimeError("hybrid initialization type did not activate")
    hybrid = initialization.hybrid_init_options
    before_options = hybrid.get_state()
    hybrid.general_settings.iter_count = iterations
    after_state = initialization.get_state()
    active_options = after_state.get("hybrid_init_options")
    if active_options is None:
        raise RuntimeError(f"hybrid option state missing after activation: {after_state}")
    actual_iterations = int(active_options["general_settings"]["iter_count"])
    if actual_iterations != iterations:
        raise RuntimeError(f"hybrid iter_count readback mismatch: {actual_iterations} != {iterations}")
    return {
        "initialization_type_before": before_type,
        "initialization_state_before": before_state,
        "hybrid_options_before_iter_change": before_options,
        "initialization_state_after": after_state,
        "requested_hybrid_iterations": iterations,
        "readback_hybrid_iterations": actual_iterations,
        "unchanged_default_hybrid_options": {
            "const_velocity": active_options["general_settings"]["initialization_options"]["const_velocity"],
            "external_aero": active_options["general_settings"]["initialization_options"]["external_aero"],
            "initial_pressure": active_options["general_settings"]["initialization_options"]["initial_pressure"],
            "explicit_urf": active_options["general_settings"]["explicit_urf"],
            "averaged_turbulent_parameters": active_options["turbulent_setting"]["averaged_turbulent_parameters"],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--processors", type=int, default=4)
    parser.add_argument(
        "--farfield-type", choices=tuple(FARFIELD_GROUP_BY_TYPE), required=True,
        help="actual Fluent type of the retained farfield physical zone",
    )
    parser.add_argument("--hybrid-iterations", type=int, required=True)
    args = parser.parse_args()
    if args.processors < 1:
        raise ValueError("--processors must be positive")
    source = args.case.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    output_case = work / "cold_air_sst_hybrid_options_configured.cas.h5"
    document: dict[str, Any] = {
        "schema": "vortex.fluent.hybrid-initialization-options-configuration/v1",
        "started_utc": stamp(),
        "status": "starting",
        "source_case": str(source),
        "source_is_immutable": True,
        "processors": args.processors,
        "farfield_boundary_type": args.farfield_type,
        "requested_hybrid_iterations": args.hybrid_iterations,
        "changed_setting": "solution.initialization.hybrid_init_options.general_settings.iter_count",
        "explicitly_unchanged_options": [
            "const_velocity", "external_aero", "initial_pressure", "explicit_urf",
            "averaged_turbulent_parameters",
        ],
        "actions": [
            "read case", "verify topology/models", "activate hybrid option tree",
            "set hybrid iter_count", "write independent case",
        ],
        "prohibited_actions": [
            "initialization", "iteration", "mesh modification", "write_data", "source-case write"
        ],
        "scope": "ideal-air density-based implicit SST initialization-options diagnostic only",
    }
    dump(work / "configuration.json", document)
    session = None
    try:
        session = pyfluent.launch_fluent(
            precision="double", processor_count=args.processors, mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI, start_timeout=1200,
        )
        emit("HYBRID_OPTIONS_CONFIG_SOLVER_LAUNCHED")
        settings = session.settings
        settings.file.read_case(file_name=str(source))
        emit("HYBRID_OPTIONS_CONFIG_SOURCE_CASE_READ_OK")
        document["topology_and_models"] = verify_topology_models(settings, args.farfield_type)
        document["hybrid_options"] = configure_hybrid_options(settings, args.hybrid_iterations)
        document["output_case"] = str(output_case)
        document["status"] = "configured"
        dump(work / "configuration.json", document)
        settings.file.write_case(file_name=str(output_case))
        document["status"] = "ok"
        document["finished_utc"] = stamp()
        dump(work / "configuration.json", document)
        emit(f"HYBRID_OPTIONS_CONFIG_CASE_WRITTEN path={output_case}")
        emit("VORTEX_HYBRID_OPTIONS_CONFIG_OK")
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
                emit("HYBRID_OPTIONS_CONFIG_SOLVER_EXIT_OK")
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
