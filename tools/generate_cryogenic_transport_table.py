#!/usr/bin/env python3
"""Generate runtime cryogenic dense-fluid transport tables.

Produces a strict, provenance-bound (T, rho) grid table of mu and kappa for
pure O2 and C3H8 using the pinned CoolProp 7.2.0 HEOS backend.  This is the
runtime input for the C++ CryogenicTransportTable reader (Phase 1b).  It is
NOT a dense-mixture model, a two-phase effective closure, or a substitution
for the frozen reference matrix.

Table design:
  - coordinate: (T, rho) because the dense-fluid state is most naturally
    defined by density, and the solver owns rho from its conserved state.
  - O2  : T in [60, 200] K (Tc=154.59 K), rho in [50, 1250] kg/m3
          (covers liquid, gas, and supercritical branches).
  - C3H8: T in [90, 400] K (Tc=369.89 K), rho in [10, 650] kg/m3
          (covers liquid, gas, and supercritical branches).
  - mu and kappa evaluated with CoolProp PropsSI on (D,T) input mode.
  - Out-of-region or two-phase (T, rho) queries are recorded as NaN and must
    fail closed in the reader; no clamping is allowed.

Every generated row binds backend, version, git revision, CAS, correlation
references, and a canonical SHA-256 of the fluid's CoolProp TRANSPORT JSON.

Schema: rfam.cryogenic-transport-runtime-table/v1
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
from typing import Any, Iterable

import CoolProp
from CoolProp import CoolProp as CP


EXPECTED_COOLPROP_VERSION = "7.2.0"
EXPECTED_COOLPROP_GIT_REVISION = "98b3523d5daa98454618d381d2ae53f7471d216b"
SCHEMA_VERSION = "rfam.cryogenic-transport-runtime-table/v1"

DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "tests"
    / "reference"
    / "cryogenic_transport_runtime_coolprop72_v1.csv"
)

FLUIDS = {
    "O2": {
        "coolprop_name": "Oxygen",
        "cas": "7782-44-7",
        "viscosity_reference": "Lemmon-IJT-2004",
        "conductivity_reference": "Lemmon-IJT-2004",
        "transport_json_sha256": (
            "17d25c4a10f8ccd5ec69d14f0fb1d3327a687423b6bdeea77ee1f79bff449b8a"
        ),
        "temperature_nodes": list(range(60, 201, 5)),     # 60..200 K step 5
        "density_nodes": [
            50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 400.0, 500.0,
            600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1250.0,
        ],
    },
    "C3H8": {
        "coolprop_name": "Propane",
        "cas": "74-98-6",
        "viscosity_reference": "Vogel-JPCRD-1998",
        "conductivity_reference": "Marsh-JCED-2002",
        "transport_json_sha256": (
            "c3a43f7b18a507f67f04e08fc28f353c0ae5a8370503cb21485b2298782afbfb"
        ),
        "temperature_nodes": list(range(90, 451, 10)),    # 90..450 K step 10
        "density_nodes": [
            10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 150.0, 200.0, 250.0,
            300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0,
        ],
    },
}


def transport_json_sha256(fluid: str) -> str:
    """Canonical SHA-256 of the fluid's CoolProp TRANSPORT JSON block."""
    block = CP.get_fluid_param_string(fluid, "TRANSPORT")
    return hashlib.sha256(block.encode("utf-8")).hexdigest()


def props_dt(fluid: str, T: float, rho: float, prop: str) -> float:
    """CoolProp PropsSI on the (D,T) input pair. NaN on any failure."""
    try:
        value = CP.PropsSI(prop, "D", rho, "T", T, fluid)
        if not math.isfinite(value):
            return float("nan")
        return float(value)
    except Exception:
        return float("nan")


def phase_dt(fluid: str, T: float, rho: float) -> str:
    try:
        return CP.PhaseSI("D", rho, "T", T, fluid)
    except Exception:
        return "error"


def generate(fluid_key: str) -> list[dict[str, Any]]:
    spec = FLUIDS[fluid_key]
    coolprop_name = spec["coolprop_name"]

    # Stage 1: collect saturation endpoints for every temperature layer and
    # the global density node set (fixed nodes plus saturation-adjacent
    # nodes from every layer, so all layers share one density axis).
    layer_info: list[tuple[float, float, float]] = []   # (T, rho_vap, rho_liq)
    global_nodes: set[float] = set(spec["density_nodes"])
    for T in spec["temperature_nodes"]:
        rho_vapor_sat = float("nan")
        rho_liquid_sat = float("nan")
        try:
            rho_vapor_sat = float(CP.PropsSI("D", "T", T, "Q", 1.0, coolprop_name))
        except Exception:
            pass
        try:
            rho_liquid_sat = float(CP.PropsSI("D", "T", T, "Q", 0.0, coolprop_name))
        except Exception:
            pass
        if not math.isfinite(rho_vapor_sat) or not math.isfinite(rho_liquid_sat):
            rho_vapor_sat = float("nan")
            rho_liquid_sat = float("nan")
        else:
            # Saturation-adjacent nodes so a compressed-liquid or stable-
            # vapor query has a single-phase neighbour on both sides.  The
            # lower-side node must stay positive (low-T propane vapour
            # densities can be far below 0.5 kg/m3).
            for value in (rho_vapor_sat, rho_liquid_sat):
                if value - 0.5 > 0.0:
                    global_nodes.add(value - 0.5)
                global_nodes.add(value + 0.5)
        layer_info.append((T, rho_vapor_sat, rho_liquid_sat))
    density_nodes = sorted(global_nodes)

    rows: list[dict[str, Any]] = []
    for (T, rho_vapor_sat, rho_liquid_sat) in layer_info:
        for rho in density_nodes:
            phase = phase_dt(coolprop_name, T, rho)
            # (T, rho) may fall in the two-phase dome: CoolProp returns
            # nonphysical values there and the CFD topology is undefined.
            # Record the phase and force NaN so the C++ reader fails closed.
            in_two_phase = phase.lower().startswith("twophase")
            if in_two_phase:
                mu = float("nan")
                kappa = float("nan")
            else:
                mu = props_dt(coolprop_name, T, rho, "VISCOSITY")
                kappa = props_dt(coolprop_name, T, rho, "CONDUCTIVITY")
            rows.append({
                "schema_version": SCHEMA_VERSION,
                "species": fluid_key,
                "coolprop_fluid": coolprop_name,
                "cas": spec["cas"],
                "backend": "HEOS",
                "coolprop_version": EXPECTED_COOLPROP_VERSION,
                "coolprop_git_revision": EXPECTED_COOLPROP_GIT_REVISION,
                "input_mode": "DT",
                "temperature_k": T,
                "rho_mass_kg_m3": rho,
                "rho_vapor_sat_kg_m3": rho_vapor_sat,
                "rho_liquid_sat_kg_m3": rho_liquid_sat,
                "phase": phase,
                "viscosity_pa_s": mu,
                "conductivity_w_m_k": kappa,
                "viscosity_reference": spec["viscosity_reference"],
                "conductivity_reference": spec["conductivity_reference"],
                "transport_json_sha256": spec["transport_json_sha256"],
                "production_ready": 0,
            })
    return rows


def write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    if not rows:
        raise RuntimeError("no rows generated")
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def canonical_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(rows: list[dict[str, Any]]) -> list[str]:
    problems: list[str] = []
    for row in rows:
        if row["schema_version"] != SCHEMA_VERSION:
            problems.append(f"schema mismatch row {row}")
        if not math.isfinite(float(row["temperature_k"])):
            problems.append(f"nonfinite T {row}")
        if not math.isfinite(float(row["rho_mass_kg_m3"])):
            problems.append(f"nonfinite rho {row}")
        mu = float(row["viscosity_pa_s"])
        kappa = float(row["conductivity_w_m_k"])
        if math.isnan(mu) != math.isnan(kappa):
            problems.append(f"mu/kappa NaN mismatch {row}")
        if not math.isnan(mu) and not (mu > 0.0):
            problems.append(f"nonpositive mu {row}")
        if not math.isnan(kappa) and not (kappa > 0.0):
            problems.append(f"nonpositive kappa {row}")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    if args.check:
        if not args.output.is_file():
            print(f"missing table: {args.output}")
            return 2
        import csv as _csv
        with open(args.output, newline="", encoding="utf-8") as handle:
            rows = list(_csv.DictReader(handle))
        problems = verify(rows)
        if problems:
            for problem in problems[:20]:
                print("PROBLEM:", problem)
            return 2
        print(f"CRYOGENIC_TRANSPORT_TABLE_CHECK_PASS rows={len(rows)} "
              f"sha256={canonical_sha256(args.output)}")
        return 0

    all_rows: list[dict[str, Any]] = []
    for fluid_key in FLUIDS:
        all_rows.extend(generate(fluid_key))

    problems = verify(all_rows)
    if problems:
        for problem in problems[:20]:
            print("PROBLEM:", problem)
        return 2

    write_csv(args.output, all_rows)
    print(f"wrote {len(all_rows)} rows to {args.output}")
    print(f"sha256={canonical_sha256(args.output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
