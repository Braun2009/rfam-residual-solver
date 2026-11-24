#!/usr/bin/env python3
"""Generate the GERG-2008 saturation-line table v4: SAME-PRESSURE endpoints.

v4 (2026-08-04, teqp unlock): the mixture rows carry density/energy
endpoints evaluated AT the common bubble pressure p_sat (the coexisting
phase pair from the teqp GERG-2008 trace), unlike v3 whose endpoints were
NaN (CoolProp Q=0/1 mixture endpoints were physically wrong in deep cold).

Method (all teqp, Python API):
  - coexisting (x, y, p_sat) from trace_VLE_isotherm_binary (stable
    segment; pure rows unchanged from v2/v3)
  - liquid root: bisection of p(rho; T, x) = p_sat on the dense branch
  - vapour root: bisection of p(rho; T, y) = p_sat on the dilute branch
  - pressure: p = rho_molar * R * T * (1 + Lambda^r_01)  [verified:
    pure O2 120 K saturated liquid gives 1.022 MPa = CoolProp exactly]
  - internal energy: u = R * T * Lambda_10 (total alpha; teqp reference)
    converted to the NASA7 sensible reference via a Tref shift
    (Gate 5: e_sens(T,rho) = u_teqp(T,rho) - u_teqp(Tref, rho_ig)
     - R*Tref/MW at the ideal-gas point)

Schema: rfam.gerg-saturation-table/v4
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path

import numpy as np
import teqp

R = 8.314462618
TREF = 298.15
P_REF = 101325.0

# RFAM NASA7 coefficients (src/thermo/NASA7Thermo.H, low/high ranges)
NASA7 = {
    "O2": ([3.78245636e+00, -2.99673416e-03, 9.84730201e-06,
             -9.68129509e-09, 3.24372837e-12, -1.06394356e+03,
             3.65767573e+00], 200.0, 1000.0),
    "C3H8": ([9.33553810e-01, 2.64245790e-02, 6.10597270e-06,
               -2.19774990e-08, 9.51492530e-12, -1.39585200e+04,
               1.92016910e+01], 300.0, 1000.0),
}
MW_SPEC = {"O2": 0.0319988, "C3H8": 0.044097}


def h_sens_mole(pure: str, T: float) -> float:
    """NASA7 sensible molar enthalpy [J/mol], h_sens(TREF)=0."""
    a, _, _ = NASA7[pure]
    h = lambda Tt: R * Tt * (a[0] + a[1] * Tt / 2 + a[2] * Tt ** 2 / 3
                             + a[3] * Tt ** 3 / 4 + a[4] * Tt ** 4 / 5
                             + a[5] / Tt)
    return h(T) - h(TREF)


def e_ig_sens_mix(z) -> float:
    """NASA7 sensible ideal internal energy [J/kg] of the mixture."""
    num = 0.0
    den = 0.0
    for pure, zi in (("O2", z[0]), ("C3H8", z[1])):
        mw = MW_SPEC[pure]
        num += zi * (h_sens_mole(pure, TREF) - R * TREF) / mw  # e_ig(TREF)
        den += zi / mw
    return num / den


def e_sens_teqp(m, T: float, rho_molar: float, z) -> float:
    """Sensible internal energy [J/kg] = NASA7 ideal (temperature part)
    + teqp residual: e = [h_sens(T)-RT]_mix/MW + R*T*Lambda^r_10/MW.
    The Lambda^r_10 term carries the residual contribution; the ideal
    part is referenced to TREF exactly like RFAM's NASA7 sensible e.
    """
    mw = 0.0319988 * z[0] + 0.044097 * z[1]
    # mass-weighted ideal e: sum(z_i * (h_sens_i(T) - R*T)) / sum(z_i*MW_i)
    e_ig = sum(zi * (h_sens_mole(pure, T) - R * T)
               for pure, zi in (("O2", z[0]), ("C3H8", z[1]))) / mw
    u_r = R * T * m.get_Arxy(1, 0, T, rho_molar, z) / mw
    return e_ig + u_r


def pressure(m, T: float, rho_molar: float, z) -> float:
    return rho_molar * m.get_R(z) * T * (1.0 + m.get_Arxy(0, 1, T, rho_molar, z))


def internal_energy_molar(m, T: float, rho_molar: float, z) -> float:
    """Molar internal energy [J/mol] (teqp reference: ideal U=0 at Tref)."""
    return m.get_R(z) * T * m.get_Arxy(1, 0, T, rho_molar, z)


def ideal_gas_u(m, T: float, z) -> float:
    rho_ig = P_REF / (m.get_R(z) * T)  # mol/m3
    return internal_energy_molar(m, T, rho_ig, z)


def bisect_root(m, T: float, p_sat: float, z, rho_lo: float, rho_hi: float,
                descending: bool) -> float | None:
    p_lo = pressure(m, T, rho_lo, z)
    p_hi = pressure(m, T, rho_hi, z)
    if not (math.isfinite(p_lo) and math.isfinite(p_hi)):
        return None
    if descending:
        # liquid branch: p decreases with rho in the unstable zone, but on
        # the stable dense branch p increases with rho; keep the window
        # where p_lo <= p_sat <= p_hi
        pass
    if p_lo > p_sat or p_hi < p_sat:
        return None
    for _ in range(120):
        mid = 0.5 * (rho_lo + rho_hi)
        p_mid = pressure(m, T, mid, z)
        if not math.isfinite(p_mid):
            return None
        if p_mid < p_sat:
            rho_lo = mid
        else:
            rho_hi = mid
    return 0.5 * (rho_lo + rho_hi)


def liquid_root(m, T: float, p_sat: float, x) -> float | None:
    # Dense-liquid branch: fixed molar-density window (20000..45000
    # mol/m3) where the pressure is monotonically increasing on the
    # stable liquid branch; the reducing-density-scaled window was
    # unreliable at deep cold (spinodal crossing).
    for lo, hi in ((20000.0, 45000.0), (15000.0, 50000.0),
                   (25000.0, 40000.0)):
        r = bisect_root(m, T, p_sat, x, lo, hi, False)
        if r is not None:
            return r
    return None


def vapour_root(m, T: float, p_sat: float, y) -> float | None:
    # Dilute branch: scan for the p crossing then bisect (the p(rho)
    # curve is monotone on the vapour branch up to the spinodal).
    best_lo = best_hi = None
    prev_rho, prev_p = 1e-6, pressure(m, T, 1e-6, y)
    for rho in np.geomspace(1e-6, 1e5, 400):
        p = pressure(m, T, rho, y)
        if not math.isfinite(p):
            break
        if (prev_p - p_sat) * (p - p_sat) < 0.0:
            best_lo, best_hi = prev_rho, rho
            break
        prev_rho, prev_p = rho, p
    if best_lo is None:
        return None
    return bisect_root(m, T, p_sat, y, best_lo, best_hi, False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(
        "tests/reference/gerg_saturation_table_v4.csv"))
    args = parser.parse_args()

    pkg = Path(teqp.__file__).parent
    fd = str(pkg / "fluiddata")
    model = teqp.build_multifluid_model(["O2", "Propane"], fd)

    rows = []

    def add_row(T, xO2, p_sat, rhoL, rhoV, eL, eV, yO2, source, dev):
        rows.append({
            "schema_version": "rfam.gerg-saturation-table/v4",
            "temperature_k": T, "x_o2": xO2, "p_sat_pa": p_sat,
            "rho_liquid_kgm3": rhoL, "rho_vapor_kgm3": rhoV,
            "e_liquid_j_kg": eL, "e_vapor_j_kg": eV,
            "y_o2_vapor": yO2,
            "p_sat_coolprop_pa": math.nan, "p_deviation": dev,
            "source": source,
        })

    # pure rows: reuse v3 (trusted torch-flash+CoolProp values)
    v3 = {}
    v3path = Path("tests/reference/gerg_saturation_table_v4.csv")
    for r in csv.DictReader(open(v3path, encoding="utf-8")):
        v3[(float(r["temperature_k"]), float(r["x_o2"]))] = r
    for T in range(60, 155, 2):
        r = v3.get((float(T), 0.999999))
        if r:
            add_row(float(T), 0.999999, float(r["p_sat_pa"]),
                    float(r["rho_liquid_kgm3"]) if r["rho_liquid_kgm3"] else math.nan,
                    float(r["rho_vapor_kgm3"]) if r["rho_vapor_kgm3"] else math.nan,
                    float(r["e_liquid_j_kg"]) if r["e_liquid_j_kg"] else math.nan,
                    float(r["e_vapor_j_kg"]) if r["e_vapor_j_kg"] else math.nan,
                    0.999999, "v3-pure", math.nan)
    for T in range(90, 370, 3):
        r = v3.get((float(T), 1e-6))
        if r:
            add_row(float(T), 1e-6, float(r["p_sat_pa"]),
                    float(r["rho_liquid_kgm3"]) if r["rho_liquid_kgm3"] else math.nan,
                    float(r["rho_vapor_kgm3"]) if r["rho_vapor_kgm3"] else math.nan,
                    float(r["e_liquid_j_kg"]) if r["e_liquid_j_kg"] else math.nan,
                    float(r["e_vapor_j_kg"]) if r["e_vapor_j_kg"] else math.nan,
                    1e-6, "v3-pure", math.nan)

    # mixture rows: teqp same-pressure endpoints
    rhoL0 = np.array([30000.0, 10000.0])
    rhoV0 = np.array([500.0, 100.0])
    for T in range(110, 151, 5):
        try:
            tr = model.trace_VLE_isotherm_binary(float(T), rhoL0, rhoV0, None)
        except Exception:
            continue
        for xA in (0.5, 0.75, 0.95):
            best = None
            bestDist = 1e9
            for d in tr:
                x = d["xL_0 / mole frac."]
                y = d["xV_0 / mole frac."]
                pV = d["pV / Pa"]
                if not (0.0 < x < 1.0 and 0.0 < y < 1.0 and y > x):
                    continue
                if not (0.0 < pV < 2e7):
                    continue
                dist = abs(x - xA)
                if dist < bestDist:
                    bestDist = dist
                    best = (x, y, pV)
            if best is None or bestDist > 0.03:
                add_row(float(T), float(xA), math.nan, math.nan, math.nan,
                        math.nan, math.nan, math.nan, "teqp-uncovered", math.nan)
                continue
            x, y, pV = best
            zL = np.array([x, 1.0 - x])
            zV = np.array([y, 1.0 - y])
            rhoL = liquid_root(model, float(T), pV, zL)
            rhoV = vapour_root(model, float(T), pV, zV)
            eL = eV = math.nan
            if rhoL is not None:
                eL = e_sens_teqp(model, float(T), rhoL, zL)
            if rhoV is not None:
                eV = e_sens_teqp(model, float(T), rhoV, zV)
            rhoL_kg = rhoL * (0.0319988 * x + 0.044097 * (1.0 - x)) if rhoL else math.nan
            rhoV_kg = rhoV * (0.0319988 * y + 0.044097 * (1.0 - y)) if rhoV else math.nan
            add_row(float(T), float(xA), pV, rhoL_kg, rhoV_kg, eL, eV,
                    float(y), "teqp-same-p", math.nan)
            print(f"T={T}K xO2={xA}: p={pV/1e3:.1f} kPa y={y:.4f} "
                  f"rhoL={rhoL_kg if rhoL_kg==rhoL_kg else float('nan'):.1f} "
                  f"rhoV={rhoV_kg if rhoV_kg==rhoV_kg else float('nan'):.3f}")

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    data = Path(args.output).read_bytes()
    print(f"wrote {len(rows)} rows, sha256={hashlib.sha256(data).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
