#!/usr/bin/env python3
"""Generate the GERG-2008 saturation-line table (bubble pressure + roots).

Uses torch-flash binary_bubble_point (validated to 0.05% at the pure-O2
limit) plus CoolProp HEOS as the independent cross-check, over the regions
where the torch-flash root solver is healthy:

  - pure O2:  T = 60..154 K (bubble = saturation)
  - pure C3H8: T = 90..369 K
  - O2/C3H8 xO2 in {0.50, 0.75, 0.95}: T = 90..150 K (health-audited levels)

Every row: T, xO2, p_sat (torch-flash), rho_liquid, rho_vapor, y_O2,
p_sat_coolprop (HEOS cross-check) -> relative deviation.

The propane-rich region (xO2 < 0.5, deep cold) is EXCLUDED: both root
solvers diverge there (2026-08-03 dual-implementation finding); the 20 NIST
ThermoML points in that region stay pending teqp C++.

Schema: rfam.gerg-saturation-table/v2

v2 (2026-08-04): e_liquid/e_vapor converted from CoolProp absolute "U" to
the RFAM NASA7 sensible reference (h_sens(298.15 K)=0; e_sens(T,rho) =
U_cp(T,rho) - U_cp(298.15, rho_ig) - R*298.15/MW).  Gate 5 reference-state
unification: the HEM two-phase internal energies are now comparable with
the NASA7 sensible e carried by the rest of the solver.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path

import CoolProp.CoolProp as CP
import torch
import torch_flash as tf

SCHEMA_VERSION = "rfam.gerg-saturation-table/v2"


def tf_bubble(eos, T: float, xO2: float):
    """torch-flash bubble point; returns (p, y) or None."""
    try:
        out = tf.binary_bubble_point(
            eos, torch.tensor(float(T)),
            torch.tensor([xO2, 1.0 - xO2]))
        return float(out.pressure), out.vapor_composition
    except Exception:
        return None


def coolprop_bubble(T: float, xO2: float):
    """CoolProp HEOS bubble pressure; returns p or None."""
    try:
        mix = f"Oxygen[{xO2}]&Propane[{1 - xO2}]"
        return CP.PropsSI("P", "T", T, "Q", 0.0, mix)
    except Exception:
        return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(
        "tests/reference/gerg_saturation_table_v4.csv"))
    args = parser.parse_args()

    eos = tf.gerg2008(("oxygen", "propane"))

    # Reference-state unification (Gate 5, 2026-08-04): CoolProp "U" is an
    # absolute internal energy with an opaque reference; RFAM uses the
    # NASA7 sensible convention (h_sens(Tref)=0 at Tref=298.15 K, so
    # e_sens(Tref, ideal gas) = -p/rho = -R*Tref/MW).  Convert every stored
    # internal energy to the sensible convention via a per-mixture shift:
    #     e_sens(T,rho) = U_cp(T,rho) - U_cp(Tref, rho_ig) - R*Tref/MW
    TREF = 298.15
    P_REF = 101325.0
    R = 8.314462618

    def sensible_shift(x_o2: float) -> float:
        mw = 0.0319988 * x_o2 + 0.044097 * (1.0 - x_o2)
        rho_ig = P_REF * mw / (R * TREF)
        if x_o2 > 0.999:
            mix = "Oxygen"
        elif x_o2 < 0.001:
            mix = "Propane"
        else:
            mix = f"Oxygen[{x_o2}]&Propane[{1 - x_o2}]"
        # (T, p) input: TP flashes are supported for mixtures, (T, rho) is
        # not ("DHSU_T_flash does not support mixtures").  298.15 K /
        # 101325 Pa is deep in the single-phase gas region for all levels.
        u_ref = CP.PropsSI("U", "T", TREF, "P", P_REF, mix)
        e_ref = -P_REF / rho_ig  # = -R*TREF/MW (NASA7 h_sens(Tref) = 0)
        return e_ref - u_ref

    rows = []
    cases = []
    # pure O2: 60..154 K step 2
    for T in range(60, 155, 2):
        cases.append((T, 0.999999))
    # pure C3H8: 90..369 K step 3
    for T in range(90, 370, 3):
        cases.append((T, 1e-6))
    # O2/C3H8 mixtures (health-audited levels): 90..150 K step 5
    for xO2 in (0.50, 0.75, 0.95):
        for T in range(90, 151, 5):
            cases.append((T, xO2))

    for T, xO2 in cases:
        b = tf_bubble(eos, T, xO2)
        pc = coolprop_bubble(T, xO2)
        if b is None:
            # torch-flash failed; record CoolProp-only row if available
            if pc is not None:
                rows.append({
                    "schema_version": "rfam.gerg-saturation-table/v2",
                    "temperature_k": T,
                    "x_o2": xO2,
                    "p_sat_pa": pc,
                    "rho_liquid_kgm3": math.nan,
                    "rho_vapor_kgm3": math.nan,
                    "e_liquid_j_kg": math.nan,
                    "e_vapor_j_kg": math.nan,
                    "y_o2_vapor": math.nan,
                    "p_sat_coolprop_pa": pc,
                    "p_deviation": math.nan,
                    "source": "coolprop-only",
                })
            continue
        p, y = b
        dev = (p - pc) / pc if pc else math.nan
        # liquid/vapour densities AND internal energies via HEOS at (T, Q);
        # internal energies converted to the NASA7 sensible reference.
        rhoL = rhoV = eL = eV = math.nan
        try:
            mix = f"Oxygen[{xO2}]&Propane[{1 - xO2}]"
            shift = sensible_shift(xO2)
            rhoL = CP.PropsSI("D", "T", T, "Q", 0.0, mix)
            rhoV = CP.PropsSI("D", "T", T, "Q", 1.0, mix)
            eL = CP.PropsSI("U", "T", T, "Q", 0.0, mix) + shift
            eV = CP.PropsSI("U", "T", T, "Q", 1.0, mix) + shift
            # physical sanity: vapour root must be > 1e-3 kg/m3 (CoolProp
            # Q=1 can return absurdly small values in deep-cold mixes)
            if not (rhoL > rhoV > 1e-3 and eL < eV):
                rhoV = rhoL = eL = eV = math.nan
        except Exception:
            pass
        rows.append({
            "schema_version": "rfam.gerg-saturation-table/v2",
            "temperature_k": T,
            "x_o2": xO2,
            "p_sat_pa": p,
            "rho_liquid_kgm3": rhoL,
            "rho_vapor_kgm3": rhoV,
            "e_liquid_j_kg": eL,
            "e_vapor_j_kg": eV,
            "y_o2_vapor": float(y[0]) if y is not None else math.nan,
            "p_sat_coolprop_pa": pc if pc else math.nan,
            "p_deviation": dev,
            "source": "torch-flash+coolprop" if pc else "torch-flash-only",
        })
        pc_kpa = f"{pc/1e3:.2f}" if pc else "None"
        dev_s = f"{dev*100:.2f}%" if pc else "n/a"
        print(f"T={T}K xO2={xO2}: p={p/1e3:.2f} kPa "
              f"coolprop={pc_kpa} kPa dev={dev_s}")

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    data = Path(args.output).read_bytes()
    print(f"wrote {len(rows)} rows, sha256={hashlib.sha256(data).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
