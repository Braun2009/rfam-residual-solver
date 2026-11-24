#!/usr/bin/env python3
"""PR vs GERG-2008 calibration audit (2026-08-04).

Quantifies the PR-EOS deviation against the GERG-2008/CoolProp reference
on the LOX/C3H8 propellant line, per the project decision:
  - production two-phase path   = GERG-based (v3 coexisting-phase table)
  - PR path (PRRhoEFlash)       = calibration/comparison tool

Audits:
  1. pure O2 saturated liquid density + bubble pressure, 60..154 K
  2. pure C3H8 same, 90..369 K
  3. binary xO2=0.75 bubble pressure, 110..150 K (PR vs teqp v3 table)

PR saturation uses the teqp canonical_PR pure-fluid VLE (liquid density
root) plus a classic-PR pressure evaluation at the vapour root.

Schema: rfam.pr-vs-gerg-audit/v1
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import CoolProp.CoolProp as CP
import numpy as np
import teqp

MW = {"O2": 0.0319988, "C3H8": 0.044097}
R = 8.314462618

CRIT = {
    "O2":   (154.581, 5.043e6, 0.0222),
    "C3H8": (369.89, 4.2512e6, 0.1521),
}


def classic_pr_p(rho_molar: float, T: float, pure: str) -> float:
    """Classic PR pressure [Pa] at molar density [mol/m3]."""
    Tc, Pc, omega = CRIT[pure]
    Tr = T / Tc
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega
    a = 0.45724 * R * R * Tc * Tc / Pc * (1.0 + kappa * (1.0 - math.sqrt(Tr))) ** 2
    b = 0.07780 * R * Tc / Pc
    v = 1.0 / rho_molar
    return R * T / (v - b) - a / (v * (v + b) + b * (v - b))


def pr_saturation(T: float, pure: str) -> tuple[float, float] | None:
    """PR saturated liquid density [kg/m3] and bubble p [Pa]."""
    Tc, Pc, ac = CRIT[pure]
    # multi-seed: deep-cold dense liquid roots need a liquid-side seed;
    # C3H8 liquid is much lighter than O2 liquid.
    seeds = [(35000.0, 500.0), (25000.0, 200.0), (40000.0, 1000.0),
             (30000.0, 300.0), (15000.0, 100.0), (20000.0, 150.0)]
    for rhoL_seed, rhoV_seed in seeds:
        try:
            model = teqp.canonical_PR([Tc], [Pc], [ac])
            r = model.pure_VLE_T(T, rhoL_seed, rhoV_seed, 200, None)
            rhoL_molar, rhoV_molar = float(r[0]), float(r[1])
            if rhoL_molar > 0 and rhoV_molar > 0 \
                    and rhoL_molar > rhoV_molar:
                p = classic_pr_p(rhoV_molar, T, pure)
                return rhoL_molar * MW[pure], p
        except Exception:
            continue
    return None


def main() -> int:
    report = {"schema": "rfam.pr-vs-gerg-audit/v1", "rows": []}

    print("=== PR vs GERG-2008 calibration audit ===")
    print(f"{'T[K]':>6} {'fluid':>6} {'PR rhoL':>10} {'HEOS rhoL':>10} "
          f"{'dRho%':>7} {'PR p_sat':>10} {'HEOS p_sat':>10} {'dP%':>7}")

    dRho_accum = []
    for pure, Tlo, Thi, step in (("O2", 60, 155, 5), ("C3H8", 90, 370, 10)):
        for T in range(Tlo, Thi, step):
            pr = pr_saturation(float(T), pure)
            if pr is None:
                continue
            rhoL_pr, p_pr = pr
            if rhoL_pr < 500.0:  # bad root (vapour-side) — discard
                continue
            try:
                rhoL_heos = CP.PropsSI("D", "T", T, "Q", 0.0, pure)
                p_heos = CP.PropsSI("P", "T", T, "Q", 0.0, pure)
            except Exception:
                continue
            dRho = 100.0 * (rhoL_pr - rhoL_heos) / rhoL_heos
            dP = 100.0 * (p_pr - p_heos) / p_heos
            dRho_accum.append((T, pure, dRho))
            print(f"{T:6.0f} {pure:>6} {rhoL_pr:10.1f} {rhoL_heos:10.1f} "
                  f"{dRho:7.1f} {p_pr:10.0f} {p_heos:10.0f} {dP:7.1f}")
            report["rows"].append({
                "T": T, "fluid": pure, "pr_rhoL_kgm3": rhoL_pr,
                "heos_rhoL_kgm3": rhoL_heos, "dRho_pct": dRho,
                "pr_p_sat_pa": p_pr, "heos_p_sat_pa": p_heos, "dP_pct": dP})

    if dRho_accum:
        mean = sum(x[2] for x in dRho_accum) / len(dRho_accum)
        print(f"\nmean PR rhoL deviation: {mean:+.1f}% "
              f"({len(dRho_accum)} points)")

    # 3: binary xO2=0.75 bubble pressure vs teqp v3 table
    print("\n=== xO2=0.75 bubble pressure: PR vs teqp v3 ===")
    v3 = {}
    for r in csv.DictReader(open(
            "tests/reference/gerg_saturation_table_v4.csv", encoding="utf-8")):
        if abs(float(r["x_o2"]) - 0.75) < 0.01:
            v3[float(r["temperature_k"])] = float(r["p_sat_pa"])

    for T in sorted(v3):
        try:
            model = teqp.canonical_PR(
                [154.581, 369.89], [5.043e6, 4.2512e6], [0.0222, 0.1521])
            z = np.array([0.75, 0.25])
            rhoL0 = np.array([20000.0, 8000.0]) * z
            rhoV0 = np.array([500.0, 150.0]) * z
            code, L, V = model.mix_VLE_Tx(
                T, rhoL0, rhoV0, z, 1e-8, 1e-8, 1e-8, 1e-8, 200)
            L = np.asarray(L); V = np.asarray(V)
            if not ((L > 0).all() and (V > 0).all()):
                print(f"{T:6.0f} PR binary flash non-physical")
                continue
            # binary classic-PR pressure at the vapour root (pseudo-pure
            # mixing: a and b from the vapour composition)
            y = V / V.sum()
            Tc_mix = 154.581 * y[0] + 369.89 * y[1]
            Pc_mix = 5.043e6 * y[0] + 4.2512e6 * y[1]
            omega_mix = 0.0222 * y[0] + 0.1521 * y[1]
            rhoV_tot = float(V.sum())
            Tr = T / Tc_mix
            kappa = 0.37464 + 1.54226 * omega_mix - 0.26992 * omega_mix ** 2
            a = 0.45724 * R * R * Tc_mix ** 2 / Pc_mix \
                * (1.0 + kappa * (1.0 - math.sqrt(Tr))) ** 2
            b = 0.07780 * R * Tc_mix / Pc_mix
            v = 1.0 / rhoV_tot
            p_pr = R * T / (v - b) - a / (v * (v + b) + b * (v - b))
            p_teqp = v3[T]
            dP = 100.0 * (p_pr - p_teqp) / p_teqp
            print(f"{T:6.0f} PR p~{p_pr:10.0f} teqp {p_teqp:10.0f} "
                  f"dP {dP:7.1f}%")
            report["rows"].append({
                "T": T, "fluid": "O2C3H8_x0.75",
                "pr_p_pa": p_pr, "teqp_p_sat_pa": p_teqp, "dP_pct": dP})
        except Exception as e:
            print(f"{T:6.0f} PR binary flash FAIL: {str(e)[:60]}")

    out = Path("tools/tmp/pr_vs_gerg_audit.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
