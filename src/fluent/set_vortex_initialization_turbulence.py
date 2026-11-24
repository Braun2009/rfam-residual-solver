#!/usr/bin/env python3
"""Create an immutable-source Vortex case with paired ambient SST init defaults.

The base ideal-air configurations use Fluent standard-initialization defaults
``k=1`` and ``omega=1``.  At 300 K / one atmosphere this implies an ambient
eddy-viscosity ratio of roughly 6.6e4, close to Fluent's 1e5 limiter.  This
utility changes only the *stored standard-initialization defaults* for k and
omega so their free-stream estimate has a requested finite viscosity ratio.

It neither initializes nor iterates.  It does not modify the mesh, source
case, or source data.  It is strictly an ideal-air low-Mach commissioning
configuration, not a real-fluid or reacting calculation.
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
from typing import Any

os.environ.setdefault("AWP_ROOT261", "/ansys_inc/v261")
os.environ.setdefault("ANSYSLMD_LICENSE_FILE", "1055@localhost")
os.environ.setdefault("ANSYSCL_FNP_PATH", os.environ["ANSYSLMD_LICENSE_FILE"])

import ansys.fluent.core as pyfluent
from ansys.fluent.core import FluentMode, UIMode

UNIVERSAL_GAS_CONSTANT_J_KMOL_K = 8314.462618
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


def get_value(state: dict[str, Any], key: str) -> float:
    entry = state[key]
    if isinstance(entry, dict):
        return float(entry["value"])
    return float(entry)


def verify_topology_and_models(settings: Any, farfield_type: str) -> dict[str, Any]:
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
    if result["solver"]["type"] != "density-based-implicit":
        raise RuntimeError(f"solver mismatch: {result['solver']}")
    if result["solver"]["time"] != "steady":
        raise RuntimeError(f"time formulation mismatch: {result['solver']}")
    if not result["energy"]["enabled"]:
        raise RuntimeError("energy is disabled")
    if result["viscous"]["model"] != "k-omega" or result["viscous"]["k_omega_model"] != "sst":
        raise RuntimeError(f"SST mismatch: {result['viscous']}")
    if result["air"]["density"]["option"] != "ideal-gas":
        raise RuntimeError(f"air density mismatch: {result['air']}")
    return result


def configure_initialization_turbulence(
    settings: Any, background_k: float, target_ratio: float
) -> dict[str, Any]:
    if not math.isfinite(background_k) or background_k <= 0.0:
        raise ValueError("--background-k must be finite and positive")
    if not math.isfinite(target_ratio) or target_ratio <= 0.0:
        raise ValueError("--target-turbulent-viscosity-ratio must be finite and positive")

    setup = settings.setup
    initialization = settings.solution.initialization
    state_before = initialization.get_state()
    defaults_before = state_before["defaults"]
    required_defaults = {"pressure", "temperature", "k", "omega"}
    missing = required_defaults.difference(defaults_before)
    if missing:
        raise RuntimeError(f"required standard-init defaults missing: {sorted(missing)}")
    if initialization.initialization_type() != "standard":
        initialization.initialization_type = "standard"
    operating_pressure_pa = float(setup.general.operating_conditions.operating_pressure())
    gauge_pressure_pa = float(defaults_before["pressure"])
    temperature_k = float(defaults_before["temperature"])
    molecular_weight_kg_kmol = get_value(setup.materials.fluid["air"].get_state(), "molecular_weight")
    viscosity_pa_s = get_value(setup.materials.fluid["air"].get_state(), "viscosity")
    absolute_pressure_pa = operating_pressure_pa + gauge_pressure_pa
    if absolute_pressure_pa <= 0.0 or temperature_k <= 0.0 or viscosity_pa_s <= 0.0:
        raise RuntimeError(
            "invalid source initialization thermodynamic state: "
            f"p_abs={absolute_pressure_pa}, T={temperature_k}, mu={viscosity_pa_s}"
        )
    density_kg_m3 = (
        absolute_pressure_pa * molecular_weight_kg_kmol
        / (UNIVERSAL_GAS_CONSTANT_J_KMOL_K * temperature_k)
    )
    omega_s_inverse = density_kg_m3 * background_k / (target_ratio * viscosity_pa_s)
    if not math.isfinite(omega_s_inverse) or omega_s_inverse <= 0.0:
        raise RuntimeError(f"invalid computed omega: {omega_s_inverse}")

    initialization.defaults["k"] = background_k
    initialization.defaults["omega"] = omega_s_inverse
    state_after = initialization.get_state()
    defaults_after = state_after["defaults"]
    actual_k = float(defaults_after["k"])
    actual_omega = float(defaults_after["omega"])
    actual_ratio = density_kg_m3 * actual_k / (actual_omega * viscosity_pa_s)
    tolerance = 1.0e-12
    if abs(actual_k - background_k) > max(tolerance, abs(background_k) * tolerance):
        raise RuntimeError(f"k-default readback mismatch: {actual_k} != {background_k}")
    if abs(actual_omega - omega_s_inverse) > max(tolerance, abs(omega_s_inverse) * tolerance):
        raise RuntimeError(f"omega-default readback mismatch: {actual_omega} != {omega_s_inverse}")
    if state_after["initialization_type"] != "standard":
        raise RuntimeError(f"standard initialization setting lost: {state_after}")
    return {
        "formula": "mu_t/mu ~= rho*k/(omega*mu)",
        "source_defaults_before": defaults_before,
        "defaults_after": defaults_after,
        "operating_pressure_pa": operating_pressure_pa,
        "standard_init_gauge_pressure_pa": gauge_pressure_pa,
        "standard_init_absolute_pressure_pa": absolute_pressure_pa,
        "standard_init_temperature_k": temperature_k,
        "air_molecular_weight_kg_kmol": molecular_weight_kg_kmol,
        "air_viscosity_pa_s": viscosity_pa_s,
        "ideal_gas_density_kg_m3": density_kg_m3,
        "requested_background_k_m2_s2": background_k,
        "requested_turbulent_viscosity_ratio": target_ratio,
        "computed_omega_s_inverse": omega_s_inverse,
        "readback_implied_turbulent_viscosity_ratio": actual_ratio,
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
    parser.add_argument(
        "--background-k", type=float, default=1.0e-4,
        help="standard-init ambient k in m^2/s^2 (default: 1e-4)",
    )
    parser.add_argument(
        "--target-turbulent-viscosity-ratio", type=float, default=10.0,
        help="paired free-stream mu_t/mu estimate (default: 10)",
    )
    args = parser.parse_args()
    if args.processors < 1:
        raise ValueError("--processors must be positive")
    source = args.case.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    output_case = work / "cold_air_sst_init_turbulence_configured.cas.h5"
    document: dict[str, Any] = {
        "schema": "vortex.fluent.initialization-turbulence-configuration/v1",
        "started_utc": stamp(),
        "status": "starting",
        "source_case": str(source),
        "source_is_immutable": True,
        "processors": args.processors,
        "farfield_boundary_type": args.farfield_type,
        "actions": [
            "read case",
            "verify topology/models",
            "set standard initialization defaults k and omega",
            "write independent case",
        ],
        "prohibited_actions": [
            "initialization", "iteration", "mesh modification", "write_data", "source-case write"
        ],
    }
    dump(work / "configuration.json", document)
    session = None
    try:
        session = pyfluent.launch_fluent(
            precision="double", processor_count=args.processors, mode=FluentMode.SOLVER,
            ui_mode=UIMode.NO_GUI, start_timeout=1200,
        )
        emit("INIT_TURB_CONFIG_SOLVER_LAUNCHED")
        settings = session.settings
        settings.file.read_case(file_name=str(source))
        emit("INIT_TURB_CONFIG_SOURCE_CASE_READ_OK")
        document["topology_and_models"] = verify_topology_and_models(
            settings, args.farfield_type
        )
        document["initialization_turbulence"] = configure_initialization_turbulence(
            settings, args.background_k, args.target_turbulent_viscosity_ratio
        )
        document["output_case"] = str(output_case)
        document["status"] = "configured"
        dump(work / "configuration.json", document)
        settings.file.write_case(file_name=str(output_case))
        document["status"] = "ok"
        document["finished_utc"] = stamp()
        dump(work / "configuration.json", document)
        emit(f"INIT_TURB_CONFIG_CASE_WRITTEN path={output_case}")
        emit("VORTEX_INIT_TURB_CONFIG_OK")
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
                emit("INIT_TURB_CONFIG_SOLVER_EXIT_OK")
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
