#!/usr/bin/env python3
"""Calibrate TwoStepCO (R2) global-mechanism parameters to GRI-30 data.

Why (2026-08-07): the TwoStepCO defaults (A1=1e10, A2=2e12, exponents
[C3H8]^0.1 [O2]^1.65) were fitted to atmospheric air flames; in the LOX
rich-oxygen high-pressure chamber the O2^1.65 exponent explodes
(qdot ~1e14-1e16 W/m3 measured, vs ~1e10-1e12 physical) and throttles the
global explicit dt to ~1e-14 s.

Method (reusing the audit approach):
  - reference IDT grid from Cantera GRI-30 (CH4 approximates C3H8 for
    global-rate order-of-magnitude; rights-clean, checked in)
  - parameterised R2 0D ignition model with free exponents
      r1 = A1 * exp(-Ea1/RT) * [C3H8]^a * [O2]^b
  - fit log10(A1), a, b (Ea1 kept at the Westbrook-Dryer value) in
    log-space least squares over the (T, p, phi) grid
  - R2 (CO oxidation) scaled by the same A-factor ratio

Output: JSON with fitted parameters + residual report.  The fitted
parameters are DIAGNOSTIC calibration only; the production path remains
the licensed skeletal mechanism.

Schema: rfam.twostepco-calibration/v1
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

try:
    import cantera as ct
except ImportError:
    print("CANTERA_UNAVAILABLE")
    raise SystemExit(2)

R = 8.314462618
E1 = 125520.0   # R1 activation (Westbrook-Dryer C3H8, J/mol)
E2 = 167360.0   # R2 activation
MW = {"C3H8": 0.04409562, "O2": 0.0319988, "CO": 0.0280101,
      "CO2": 0.0440095, "H2O": 0.01801528, "CH4": 0.016043}
DH1 = 2.046e6   # J/mol C3H8
DH2 = 2.83e5    # J/mol CO


def gri30_idt(T0: float, p0: float, phi: float) -> float | None:
    """GRI-30 (CH4) 0D adiabatic ignition delay [s] (10% temperature rise)."""
    try:
        gas = ct.Solution("gri30.yaml")
        gas.TP = T0, p0
        gas.set_equivalence_ratio(phi, "CH4", "O2:1.0")
        gas.equilibrate("UV") if False else None
        r = ct.IdealGasReactor(gas)
        sim = ct.ReactorNet([r])
        sim.rtol, sim.atol = 1e-8, 1e-14
        T0g = gas.T
        t = 0.0
        while t < 10.0:
            sim.advance(t + 1e-4 * max(1.0, t))
            t = sim.time
            if r.thermo.T > T0g * 1.1:
                return t
        return None
    except Exception:
        return None


def r2_0d_idt(T0: float, p0: float, phi: float,
              A1: float, a: float, b: float) -> float:
    """TwoStepCO R1 0D ignition delay [s]: fuel-depletion time scale.

    r1 = A1 exp(-E1/RT) [C3H8]^a [O2]^b  (mol/m3/s)
    Fuel depletion: d[C3H8]/dt = -r1 -> tau = [C3H8]/(r1*a_fuel)
    with a_fuel = 1 for the fuel exponent scaling of the rate on C3H8.
    IDT ~ 1/(A1 exp(-E1/RT) [C3H8]^(a-1) [O2]^b)  (to 10% fuel consumed).
    """
    xC3 = 1.0
    xO2 = 3.5 / phi
    MWmix = (MW["C3H8"] * xC3 + MW["O2"] * xO2) / (xC3 + xO2)
    rho = p0 * MWmix / (R * T0)
    cC3 = rho * xC3 / MW["C3H8"]
    cO2 = rho * xO2 / MW["O2"]
    k = A1 * math.exp(-E1 / (R * T0))
    r1 = k * cC3 ** a * cO2 ** b
    if r1 <= 0.0:
        return float("inf")
    tau10 = 0.1 * cC3 / r1
    return tau10


def main() -> int:
    grid = []
    for T0 in (1000, 1200, 1500, 1800):
        for p0 in (101325.0, 4.0e6):
            for phi in (0.8, 1.0, 1.3):
                grid.append((T0, p0, phi))

    ref = []
    for (T0, p0, phi) in grid:
        idt = gri30_idt(T0, p0, phi)
        if idt is None or not math.isfinite(idt):
            continue
        ref.append((T0, p0, phi, idt))
    print(f"GRI-30 IDT grid: {len(ref)}/{len(grid)} valid")

    # fit log10(A1), a, b in log space
    Y = np.array([math.log(r[-1]) for r in ref])
    A = np.zeros((len(ref), 3))
    for i, (T0, p0, phi, idt) in enumerate(ref):
        xC3 = 1.0
        xO2 = 3.5 / phi
        MWmix = (MW["C3H8"] * xC3 + MW["O2"] * xO2) / (xC3 + xO2)
        rho = p0 * MWmix / (R * T0)
        cC3 = rho * xC3 / MW["C3H8"]
        cO2 = rho * xO2 / MW["O2"]
        A[i, 0] = 1.0                      # log10(A1) coefficient
        A[i, 1] = math.log(0.1 * cC3)      # (a-1) coefficient: tau = 0.1 cC3 / (A1 k cC3^a cO2^b)
        A[i, 2] = math.log(cO2)
        Y[i] = math.log(idt) - (-E1 / (R * T0)) - math.log(1.0)
        # tau = 0.1*cC3 / (A1*k*cC3^a*cO2^b) = (0.1*cC3)/A1/k * cC3^-a * cO2^-b
        # ln tau = ln(0.1*cC3) - ln A1 - ln k - a ln cC3 - b ln cO2
        # -> ln(0.1*cC3) - ln tau - ln k = ln A1 + a ln cC3 + b ln cO2
        A[i, 1] = math.log(cC3)
        A[i, 2] = math.log(cO2)
        Y[i] = math.log(0.1 * cC3) - math.log(idt) \
            + E1 / (R * T0) - math.log(1.0)

    coef, res, *_ = np.linalg.lstsq(A, Y, rcond=None)
    lnA1, a, b = coef
    A1 = math.exp(lnA1)
    print(f"fit: A1={A1:.3e}  a={a:.3f}  b={b:.3f}  (E1={E1/1e3:.1f} kJ/mol)")

    # residual report
    print(f"{'T0':>5} {'p':>8} {'phi':>4} {'IDT_GRI':>10} {'IDT_R2':>10} {'ratio':>8}")
    ratios = []
    for (T0, p0, phi, idt) in ref:
        idt_r2 = r2_0d_idt(T0, p0, phi, A1, a, b)
        ratio = idt_r2 / idt if idt > 0 else float("nan")
        ratios.append(ratio)
        print(f"{T0:5.0f} {p0:8.0f} {phi:4.1f} {idt:10.3e} {idt_r2:10.3e} {ratio:8.2f}")

    mean_ratio = float(np.mean(np.abs(np.log(np.clip(ratios, 1e-9, 1e9)))))
    out = {
        "schema": "rfam.twostepco-calibration/v1",
        "mechanism_reference": "gri30-3.0 (CH4 approximates C3H8)",
        "fitted": {"A1": A1, "a": a, "b": b, "E1": E1},
        "A2_scale": A1 / 1.0e10,
        "log_mean_abs_ratio": mean_ratio,
        "cases": [
            {"T0": T0, "p0": p0, "phi": phi, "idt_gri30": idt,
             "idt_r2": r2_0d_idt(T0, p0, phi, A1, a, b)}
            for (T0, p0, phi, idt) in ref],
    }
    outpath = Path("tools/tmp/twostepco_calibration.json")
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {outpath}; log-mean-abs-ratio={mean_ratio:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
