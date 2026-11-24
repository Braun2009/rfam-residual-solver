#!/usr/bin/env python3
"""R2 (TwoStepCO) vs GRI-30 ignition-delay audit (2026-08-03).

Re-calibration status audit: quantifies the R2 diagnostic-mechanism IDT
gap against the Cantera GRI-30 detailed mechanism (rights-clean, checked
into tests/reference/) across a coarse (T, p, phi) grid.  Output feeds the
UNVALIDATED_KINETICS audit marker; it does NOT re-fit R2 parameters.

R2 rates (TwoStepCO defaults):
  C3H8 + 3.5 O2 -> 3 CO + 4 H2O    A1=1e10  Ea1=125.52 kJ/mol
  2 CO + O2    -> 2 CO2            A2=2e12  Ea2=167.36 kJ/mol
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

try:
    import cantera as ct
except ImportError:
    print("CANTERA_UNAVAILABLE")
    sys.exit(2)


R = 8.314462618


def r2_rates(T: float, Y: dict, rho: float):
    """Mass-basis R2 rates (mol/m3/s) from TwoStepCO defaults."""
    # 2026-08-07 GRI-30 calibration (tools/calibrate_twostepco.py):
    # A1=7.77e5, exponents a=1.303 (C3H8), b=0.263 (O2), A2 scaled 7.77e-5.
    A1, n1, Ea1 = 7.77e5, 0.0, 125520.0
    A2, n2, Ea2 = 1.55e8, 0.0, 167360.0
    MW = {"C3H8": 0.04409562, "O2": 0.0319988, "CO": 0.0280101,
          "CO2": 0.0440095, "H2O": 0.01801528}
    xC3 = Y.get("C3H8", 0.0) / MW["C3H8"]
    xO2 = Y.get("O2", 0.0) / MW["O2"]
    xCO = Y.get("CO", 0.0) / MW["CO"]
    r1 = A1 * math.exp(-Ea1 / (R * T)) * (rho * xC3) ** 1.303 * (rho * xO2) ** 0.263
    r2 = A2 * math.exp(-Ea2 / (R * T)) * (rho * xCO) ** 2.0 * (rho * xO2) ** 1.0
    return r1, r2


def r2_idt(T0: float, p0: float, phi: float):
    """R2 0D adiabatic ignition delay (C3H8 + 3.5 O2 -> 3 CO + 4 H2O)."""
    xC3 = 1.0
    xO2 = 3.5 / phi
    MW = {"C3H8": 0.04409562, "O2": 0.0319988, "CO": 0.0280101,
          "CO2": 0.0440095, "H2O": 0.01801528}
    n_tot = xC3 + xO2
    Y = {"C3H8": xC3 * MW["C3H8"] / (xC3 * MW["C3H8"] + xO2 * MW["O2"]),
         "O2": xO2 * MW["O2"] / (xC3 * MW["C3H8"] + xO2 * MW["O2"])}
    MW_mix = Y["C3H8"] / (xC3 / n_tot) if False else 0.0
    MW_mix = (xC3 * MW["C3H8"] + xO2 * MW["O2"]) / n_tot
    rho = p0 * MW_mix / (R * T0)
    dt = 1e-9
    T = T0
    t = 0.0
    Cp = 1750.0  # J/kg/K mixture estimate
    qchem = 46.0e6  # J/kg fuel (C3H8 LHV)
    while t < 1.0:
        r1, r2 = r2_rates(T, Y, rho)
        fuel_mol_rate = r1  # mol/m3/s C3H8
        fuel_mass_rate = fuel_mol_rate * MW["C3H8"]
        qdot = fuel_mass_rate * qchem
        dT = qdot / (rho * Cp)
        T += dT * dt
        if T - T0 > 400.0:
            return t, "r2"
        t += dt
    return None, "r2"


def gri30_idt(T0: float, p0: float, phi: float):
    gas = ct.Solution("tests/reference/cantera_gri30_v3_0.yaml")
    gas.TPX = T0, p0, {"C3H8": 1.0, "O2": 5.0 / phi, "N2": 18.8 / phi}
    r = ct.IdealGasReactor(gas)
    sim = ct.ReactorNet([r])
    sim.rtol = 1e-10
    sim.atol = 1e-16
    sim.max_time_step = 1e-5
    t = 0.0
    T0r = gas.T
    while sim.time < 1.0:
        sim.step()
        if gas.T > T0r + 400.0:
            return sim.time, "gri30"
    return None, "gri30"


def main() -> int:
    cases = []
    for T0 in (1100, 1300, 1600):
        for p0_atm in (1, 10, 20):
            for phi in (0.5, 1.0, 2.0):
                cases.append((T0, p0_atm * ct.one_atm, phi))
    rows = []
    for T0, p0, phi in cases:
        t_r2, _ = r2_idt(T0, p0, phi)
        t_g, _ = gri30_idt(T0, p0, phi)
        ratio = (t_r2 / t_g) if (t_r2 and t_g) else None
        rows.append({
            "T0": T0, "p0_atm": p0 / ct.one_atm, "phi": phi,
            "idt_r2_s": t_r2, "idt_gri30_s": t_g,
            "ratio_r2_over_gri30": ratio,
        })
        print(f"T={T0}K p={p0/ct.one_atm:.0f}atm phi={phi}: "
              f"R2={t_r2*1e3:.4f}ms GRI30={t_g*1e3:.4f}ms "
              f"ratio={ratio:.1f}x" if t_r2 and t_g else
              f"T={T0}K p={p0/ct.one_atm:.0f}atm phi={phi}: MISSING")
    out = Path("tmp/r2_gri30_idt_audit_20260803.json")
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps({"schema": "rfam.r2-kinetics-audit/v1",
                               "mechanism": "gri30-3.0",
                               "cases": rows}, indent=1))
    print(f"wrote {out} ({len(rows)} cases)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
