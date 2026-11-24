#!/usr/bin/env python3
"""Run the first 3D Fluent SST cold-air commissioning case on the final mesh.

This is deliberately *not* a LOX/C3H8 physical prediction.  It validates the
new full-360 mesh, its five physical boundaries, the automatic mm-to-m Solver
conversion, density-based SST, and solver restart artifacts before the GERG /
species / EDC path is enabled.

The default mass-flow scale (0.0025) turns the design mass flows into a
roughly velocity-matched ambient-air commissioning condition.  Do not label
its pressure, thrust, or heat transfer as a rocket operating-point result.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
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


MESH = Path("/root/vortex_fluent/gmsh_fixed_v2/strict_tui_prism_final.msh.h5")
MESH_SHA256 = "2640bb1637e2420b583920e7054906fab1545a53f0b1867287ef81cdbaba7a36"
DEFAULT_WORK = Path("/root/vortex_fluent/gmsh_fixed_v2/solver_cold_air_sst_v1")

# Simplified one-spray-cone geometry: main LOX+C3H8 total is intentionally a
# single zone; the separate centre torch remains a second inlet zone.
DESIGN_MDOT_KG_S = {
    "inlet_spray": 2.8 + 1.05,
    "inlet_torch": 0.1348,
}
REQUIRED_BOUNDARIES = {
    "mass_flow_inlet": {"inlet_spray", "inlet_torch"},
    "pressure_far_field": {"farfield"},
    "wall": {"wall_chamber", "wall_plate"},
}
REQUIRED_FLUID_ZONES = {
    "vortex_full360_v2_zoned_ascii-prism-cells",
    "vortex_full360_v2_zoned_ascii-tet-cells",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def log(message: str) -> None:
    print(message, flush=True)


def write_json(path: Path, document: dict[str, Any]) -> None:
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def named_objects(group: Any) -> set[str]:
    return set(group.get_object_names())


def assert_mesh_topology(settings: Any) -> dict[str, list[str]]:
    boundary_conditions = settings.setup.boundary_conditions
    inventory: dict[str, list[str]] = {}
    for group_name, required_names in REQUIRED_BOUNDARIES.items():
        actual = named_objects(getattr(boundary_conditions, group_name))
        inventory[group_name] = sorted(actual)
        if actual != required_names:
            raise RuntimeError(
                f"{group_name} mismatch: expected {sorted(required_names)}, got {sorted(actual)}"
            )
    cell_zones = settings.setup.cell_zone_conditions
    fluid_names = named_objects(cell_zones.fluid)
    inventory["fluid"] = sorted(fluid_names)
    if fluid_names != REQUIRED_FLUID_ZONES:
        raise RuntimeError(
            f"fluid-zone mismatch: expected {sorted(REQUIRED_FLUID_ZONES)}, got {sorted(fluid_names)}"
        )
    return inventory


def configure_models(settings: Any, courant: float) -> None:
    setup = settings.setup
    setup.materials.fluid["air"].density.option = "ideal-gas"
    setup.models.energy.enabled = True
    setup.models.viscous.model = "k-omega"
    setup.models.viscous.k_omega_model = "sst"
    setup.general.solver.time = "steady"
    # Fluent 26.1's exact density-based enum; not the legacy "density-based".
    setup.general.solver.type = "density-based-implicit"
    setup.general.operating_conditions.operating_pressure = 101325.0

    # Preserve the two conformal cell zones and make their common material
    # explicit, rather than relying on imported defaults.
    for zone_name in REQUIRED_FLUID_ZONES:
        setup.cell_zone_conditions.fluid[zone_name].general.material = "air"

    controls = settings.solution.controls
    controls.courant_number = courant
    controls.limits.min_pressure = 1.0
    controls.limits.min_temperature = 1.0
    controls.limits.max_temperature = 6000.0


def configure_boundaries(
    settings: Any,
    mass_scale: float,
    temperature_k: float,
    turbulence_intensity: float,
    turbulence_viscosity_ratio: float,
    farfield_pressure_pa_abs: float,
    farfield_mach: float,
) -> dict[str, Any]:
    if not 0.0 < mass_scale <= 1.0:
        raise ValueError(f"mass_scale must be in (0, 1], got {mass_scale}")
    if farfield_pressure_pa_abs <= 0.0:
        raise ValueError("farfield absolute pressure must be positive")
    if not 0.0 < farfield_mach < 1.0:
        raise ValueError(f"farfield_mach must be in (0, 1), got {farfield_mach}")

    setup = settings.setup
    bc = setup.boundary_conditions
    inlets = bc.mass_flow_inlet
    assigned_mdot: dict[str, float] = {}
    for zone_name, design_mdot in DESIGN_MDOT_KG_S.items():
        inlet = inlets[zone_name]
        assigned = design_mdot * mass_scale
        inlet.momentum.mass_flow_specification = "Mass Flow Rate"
        inlet.momentum.mass_flow_rate.value = assigned
        # A planar inlet at x=0: Fluent chooses the inward normal.  Cone/swirl
        # velocity vectors belong to the later species/real-fluid injector case.
        inlet.momentum.direction_specification = "Normal to Boundary"
        inlet.thermal.total_temperature.value = temperature_k
        inlet.turbulence.turbulence_specification = "Intensity and Viscosity Ratio"
        inlet.turbulence.turbulent_intensity = turbulence_intensity
        inlet.turbulence.turbulent_viscosity_ratio = turbulence_viscosity_ratio
        assigned_mdot[zone_name] = assigned

    farfield = bc.pressure_far_field["farfield"]
    farfield.momentum.gauge_pressure.value = (
        farfield_pressure_pa_abs - setup.general.operating_conditions.operating_pressure()
    )
    # Fluent pressure-far-field requires strictly positive Mach during flow initialization.
    farfield.momentum.mach_number.value = farfield_mach
    farfield.thermal.temperature.value = temperature_k
    farfield.turbulence.turbulence_specification = "Intensity and Viscosity Ratio"
    farfield.turbulence.turbulent_intensity = turbulence_intensity
    farfield.turbulence.turbulent_viscosity_ratio = turbulence_viscosity_ratio

    # Explicit adiabatic, stationary, no-slip walls for this non-reacting chain
    # validation.  The thermal CHT/wall-temperature model comes after cold flow.
    for zone_name in ("wall_chamber", "wall_plate"):
        wall = bc.wall[zone_name]
        wall.momentum.wall_motion = "Stationary Wall"
        wall.momentum.shear_condition = "No Slip"
        wall.thermal.thermal_condition = "Heat Flux"
        wall.thermal.heat_flux.value = 0.0

    return {
        "mass_scale_of_design": mass_scale,
        "inlet_mass_flow_kg_s": assigned_mdot,
        "inlet_total_temperature_k": temperature_k,
        "farfield_absolute_pressure_pa": farfield_pressure_pa_abs,
        "farfield_gauge_pressure_pa": farfield.momentum.gauge_pressure.value(),
        "farfield_mach": farfield.momentum.mach_number.value(),
        "farfield_mach_note": "small positive initialization value required by Fluent pressure-far-field", 
        "turbulence": {
            "specification": "Intensity and Viscosity Ratio",
            "intensity": turbulence_intensity,
            "viscosity_ratio": turbulence_viscosity_ratio,
        },
        "walls": "stationary_no_slip_adiabatic",
    }


def readback(settings: Any) -> dict[str, Any]:
    setup = settings.setup
    bc = setup.boundary_conditions
    result: dict[str, Any] = {
        "solver": setup.general.solver.get_state(),
        "energy": setup.models.energy.get_state(),
        "viscous": setup.models.viscous.get_state(),
        "air": setup.materials.fluid["air"].get_state(),
        "cell_materials": {
            zone: setup.cell_zone_conditions.fluid[zone].general.material()
            for zone in sorted(REQUIRED_FLUID_ZONES)
        },
        "inlets": {},
        "farfield": bc.pressure_far_field["farfield"].get_state(),
        "walls": {
            name: bc.wall[name].get_state()
            for name in ("wall_chamber", "wall_plate")
        },
        "solution_methods": settings.solution.methods.get_state(),
        "solution_controls": settings.solution.controls.get_state(),
    }
    for name in sorted(DESIGN_MDOT_KG_S):
        result["inlets"][name] = bc.mass_flow_inlet[name].get_state()
    if result["solver"]["type"] != "density-based-implicit":
        raise RuntimeError(f"density-based solver setting lost: {result['solver']}")
    if result["air"]["density"]["option"] != "ideal-gas":
        raise RuntimeError(f"ideal-gas setting lost: {result['air']}")
    if not result["energy"]["enabled"]:
        raise RuntimeError("energy equation is disabled")
    if result["viscous"]["model"] != "k-omega" or result["viscous"]["k_omega_model"] != "sst":
        raise RuntimeError(f"SST setting lost: {result['viscous']}")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--processors", type=int, default=8)
    parser.add_argument("--mass-scale", type=float, default=0.0025)
    parser.add_argument("--iterations", type=int, default=20)
    parser.add_argument("--courant", type=float, default=1.0)
    parser.add_argument(
        "--initialization",
        choices=("standard", "hybrid"),
        default="standard",
        help=(
            "initialization method; standard is the commissioning default because "
            "hybrid initialization entered an interrupted-compute state on this 7.27M-cell mesh"
        ),
    )
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--farfield-pressure-pa", type=float, default=101325.0)
    parser.add_argument(
        "--farfield-mach", type=float, default=0.01,
        help="strictly positive pressure-far-field Mach required by Fluent initialization",
    )
    parser.add_argument("--turbulence-intensity", type=float, default=0.05)
    parser.add_argument("--turbulence-viscosity-ratio", type=float, default=10.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.processors < 1:
        raise ValueError("--processors must be at least 1")
    if args.iterations < 1:
        raise ValueError("--iterations must be at least 1")
    if args.courant <= 0.0:
        raise ValueError("--courant must be positive")

    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    if not MESH.is_file():
        raise FileNotFoundError(MESH)
    actual_sha = sha256(MESH)
    if actual_sha != MESH_SHA256:
        raise RuntimeError(f"final mesh SHA256 mismatch: {actual_sha} != {MESH_SHA256}")

    manifest: dict[str, Any] = {
        "schema": "vortex.fluent.cold-air-sst-commissioning/v1",
        "started_utc": utc_now(),
        "status": "starting",
        "mesh": {"path": str(MESH), "sha256": actual_sha, "hdf5_units": "mm", "solver_units": "m"},
        "intent": (
            "Mesh/solver commissioning only; ideal air, no species, no reactions, "
            "no GERG UDRGM, no EDC, no physical LOX/C3H8 performance claim."
        ),
        "arguments": {
            **vars(args),
            "work": str(args.work),
        },
        "design_mass_flow_kg_s": DESIGN_MDOT_KG_S,
    }
    write_json(work / "manifest.json", manifest)

    session = None
    try:
        session = pyfluent.launch_fluent(
            precision="double",
            processor_count=args.processors,
            mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI,
            start_timeout=1200,
        )
        log("COLD_AIR_SOLVER_LAUNCHED")
        settings = session.settings
        settings.file.read_mesh(file_name=str(MESH))
        log("COLD_AIR_FINAL_MESH_READ_OK")
        manifest["topology"] = assert_mesh_topology(settings)

        configure_models(settings, args.courant)
        manifest["boundary_assignment"] = configure_boundaries(
            settings,
            args.mass_scale,
            args.temperature_k,
            args.turbulence_intensity,
            args.turbulence_viscosity_ratio,
            args.farfield_pressure_pa,
            args.farfield_mach,
        )
        manifest["readback_before_initialization"] = readback(settings)
        manifest["status"] = "configured"
        write_json(work / "manifest.json", manifest)
        settings.file.write_case(file_name=str(work / "cold_air_sst_configured.cas.h5"))
        log("COLD_AIR_CONFIGURED_CASE_WRITTEN")

        if args.initialization == "hybrid":
            settings.solution.initialization.initialization_type = "hybrid"
            settings.solution.initialization.hybrid_initialize()
        else:
            settings.solution.initialization.standard_initialize()
        log(f"COLD_AIR_{args.initialization.upper()}_INITIALIZATION_OK")
        settings.file.write_case(file_name=str(work / "cold_air_sst_initialized.cas.h5"))
        settings.file.write_data(file_name=str(work / "cold_air_sst_initialized.dat.h5"))
        log("COLD_AIR_INITIALIZED_CASE_DATA_WRITTEN")

        settings.solution.run_calculation.iterate(iter_count=args.iterations)
        log(f"COLD_AIR_ITERATION_OK count={args.iterations}")
        manifest["readback_after_iterations"] = readback(settings)
        manifest["status"] = "iterated"
        manifest["finished_utc"] = utc_now()
        write_json(work / "manifest.json", manifest)
        settings.file.write_case(file_name=str(work / "cold_air_sst_commissioning.cas.h5"))
        settings.file.write_data(file_name=str(work / "cold_air_sst_commissioning.dat.h5"))
        log("COLD_AIR_FINAL_CASE_DATA_WRITTEN")
        log("VORTEX_3D_COLD_AIR_SST_COMMISSIONING_OK")
    except Exception:
        manifest["status"] = "failed"
        manifest["finished_utc"] = utc_now()
        manifest["traceback"] = traceback.format_exc()
        write_json(work / "manifest.json", manifest)
        raise
    finally:
        if session is not None:
            try:
                session.exit()
                log("COLD_AIR_SOLVER_EXIT_OK")
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
