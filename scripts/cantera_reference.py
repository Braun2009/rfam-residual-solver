#!/usr/bin/env python3
"""
cantera_reference.py — Cantera reference values for RFAM 1D flame comparison.

Computes using GRI-Mech 3.0:
  1. Adiabatic flame temperature (Tad) for stoich C3H8/O2 at several p, T0
  2. Ignition delay (tau_ign) in a constant-pressure 0D reactor
  3. 1D laminar flame speed (Su) at 1 atm, 300K (reference point)

Output: JSON to stdout, consumed by test_cantera_compare.cpp driver.

Usage: python3 scripts/cantera_reference.py
"""

import sys
import json
import cantera as ct

def adiabatic_T(fuel, oxidizer, p_Pa, T0_K, phi=1.0):
    """Cantera equilibrium Tad at constant pressure."""
    g = ct.Solution("gri30.yaml")
    g.TP = T0_K, p_Pa
    g.set_equivalence_ratio(phi, fuel, oxidizer)
    g.equilibrate("HP")
    return g.T

def ignition_delay(fuel, oxidizer, p_Pa, T0_K, phi=1.0, t_end=1e-2):
    """0-D constant-pressure ignition delay (time to 400K rise)."""
    g = ct.Solution("gri30.yaml")
    g.TP = T0_K, p_Pa
    g.set_equivalence_ratio(phi, fuel, oxidizer)
    r = ct.IdealGasConstPressureReactor(g)
    sim = ct.ReactorNet([r])
    T_init = r.thermo.T
    t = 0.0
    tau = None
    while t < t_end:
        t = sim.step()
        if r.thermo.T - T_init > 400.0 and tau is None:
            tau = t
            break
    return tau  # None if no ignition within t_end

def laminar_flame_speed(fuel, oxidizer, p_Pa, T0_K, phi=1.0):
    """1D freely-propagating premixed flame speed [m/s]."""
    g = ct.Solution("gri30.yaml")
    g.TP = T0_K, p_Pa
    g.set_equivalence_ratio(phi, fuel, oxidizer)
    f = ct.FreeFlame(g, width=0.03)
    f.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)
    try:
        f.solve(loglevel=0, auto=True)
        return f.velocity[0]
    except Exception:
        return None

if __name__ == "__main__":
    fuel = "C3H8"
    oxidizer = "O2"  # pure O2 (no N2) — matches RFAM TwoStepCO

    results = {}

    # ---- 1. Tad sweep: several p, T0=300K, phi=1.0 ----
    tad = {}
    for p_MPa in [1.0, 4.0, 8.0]:
        T = adiabatic_T(fuel, oxidizer, p_MPa * 1e6, 300.0)
        tad[str(p_MPa)] = round(T, 1)
    results["tad_phi1_T300"] = tad

    # ---- 2. Ignition delay: T0=1200K, 4MPa, phi=1.0 (matches RFAM A-box IC) ----
    tau = ignition_delay(fuel, oxidizer, 4e6, 1200.0)
    results["tau_ign_1200K_4MPa"] = round(tau * 1e6, 2) if tau else None  # µs

    # ---- 3. 1D flame speed: 1atm, 300K, phi=1.0 ----
    Su = laminar_flame_speed(fuel, oxidizer, ct.one_atm, 300.0)
    results["Su_1atm_300K"] = round(Su, 4) if Su else None

    print(json.dumps(results, indent=2))
