#!/usr/bin/env python3
"""GERG-2008 (torch-flash) vs PR vs NIST ThermoML O2/C3H8 bubble pressure.

Gate 2 of the GERG adoption plan.  Compares predicted bubble pressure at the
four experimental temperatures (110.22/120.13/130.58/139.95 K) against the
28-row NIST ThermoML P-x dataset (DOI 10.1021/je100574y).

Models:
  - GERG-2008 via torch-flash (independent Helmholtz multifluid backend)
  - Peng-Robinson via manual cubic solver with kij=0.055 (historical RFAM)
    and kij=0.0 (no-interaction control)

Usage: python tools/compare_gerg_o2_c3h8_vle.py \
           --csv F:/tmp/a408-o2-c3h8-thermoml-candidate.csv [--check]
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path

import torch
import torch_flash as tf


def gerg_bubble(T: float, xO2: float, p_hint: float | None = None) -> float:
    eos = tf.gerg2008(("oxygen", "propane"))
    try:
        out = tf.binary_bubble_point(
            eos, torch.tensor(float(T)), torch.tensor([xO2, 1.0 - xO2]))
        if out.converged:
            return float(out.pressure)
    except Exception:
        pass
    return float("nan")


def _pr_bubble_pressure_manual(T: float, xO2: float, kij: float) -> float:
    R = 8.314462618
    Tc = {"O2": 154.581, "C3H8": 369.89}
    Pc = {"O2": 5.043e6, "C3H8": 4.2512e6}
    w = {"O2": 0.0222, "C3H8": 0.152}
    names = ["O2", "C3H8"]
    x = [xO2, 1.0 - xO2]

    def alpha(name, t):
        kappa = 0.37464 + 1.54226*w[name] - 0.26992*w[name]**2
        tr = t / Tc[name]
        return (1.0 + kappa*(1.0 - math.sqrt(tr)))**2

    def a_ij(i, j):
        ai = 0.45724*R*R*Tc[names[i]]**2/Pc[names[i]]*alpha(names[i], T)
        aj = 0.45724*R*R*Tc[names[j]]**2/Pc[names[j]]*alpha(names[j], T)
        k = kij if i != j else 0.0
        return math.sqrt(ai*aj)*(1.0 - k)

    def fugacity_coeffs(p):
        am = 0.0
        for i in range(2):
            for j in range(2):
                am += x[i]*x[j]*a_ij(i, j)
        bm = sum(x[i]*0.07780*R*Tc[names[i]]/Pc[names[i]] for i in range(2))
        A = am*p/(R*R*T*T)
        B = bm*p/(R*T)
        c2 = -(1.0 - B)
        c1 = A - 3.0*B*B - 2.0*B
        c0 = -(A*B - B*B - B**3)
        roots = _cubic_roots(1.0, c2, c1, c0)
        if not roots:
            return None
        # Bubble: liquid-phase fugacity coefficients -> smallest Z root.
        Z = min(roots)
        if Z - B <= 0.0:
            return None
        phi = []
        for i in range(2):
            bi = 0.07780*R*Tc[names[i]]/Pc[names[i]]
            sumj = sum(x[j]*a_ij(i, j) for j in range(2))
            lnphi = bi/bm*(Z-1.0) - math.log(Z-B)
            term = (Z+(1.0+math.sqrt(2.0))*B)/(Z-(1.0-math.sqrt(2.0))*B)
            if term <= 0.0:
                return None
            lnphi -= A/(2.0*math.sqrt(2.0)*B) * (2.0*sumj/am - bi/bm) \
                * math.log(term)
            phi.append(math.exp(lnphi))
        return phi

    p_lo, p_hi = 1e2, 5e7
    for _ in range(80):
        p_mid = 0.5*(p_lo+p_hi)
        phi = fugacity_coeffs(p_mid)
        if phi is None:
            p_hi = p_mid
            continue
        s = sum(x[i]*phi[i] for i in range(2))
        if s > 1.0:
            p_hi = p_mid
        else:
            p_lo = p_mid
    return 0.5*(p_lo+p_hi)


def _cubic_roots(a, b, c, d):
    if abs(a) < 1e-30:
        return []
    p = (3.0*c/a - b*b/(a*a))/3.0
    q = (2.0*b**3/(a**3) - 9.0*b*c/(a*a) + 27.0*d/a)/27.0
    disc = q*q/4.0 + p**3/27.0
    roots = []
    if disc >= 0.0:
        u = math.copysign((abs(q/2.0)+math.sqrt(disc))**(1.0/3.0), q/2.0)
        v = -p/(3.0*u) if abs(u) > 1e-30 else 0.0
        roots.append(u+v - b/(3.0*a))
    else:
        r = math.sqrt(-p**3/27.0)
        theta = math.acos(-q/(2.0*r)) if r > 0 else 0.0
        for k in range(3):
            roots.append(2.0*math.sqrt(-p/3.0)*math.cos(
                theta/3.0 + 2.0*math.pi*k/3.0) - b/(3.0*a))
    return [r for r in roots if r > 0.0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    rows = []
    with open(args.csv, newline="", encoding="utf-8") as f:
        lines = [ln for ln in f if not ln.lstrip().startswith("#")]
        for row in csv.DictReader(lines):
            if row.get("liquid_x_o2", "").strip():
                rows.append(row)

    print(f"{'T[K]':>7} {'xO2':>7} {'Pexp[kPa]':>10} {'u(P)':>6} "
          f"{'P_GERG':>10} {'P_PR055':>10} {'eGERG%':>7} {'ePR%':>7}")
    err_gerg, err_pr055 = [], []
    for r in rows:
        T = float(r["temperature_k"])
        xO2 = float(r["liquid_x_o2"])
        Pexp = float(r["pressure_kpa"]) * 1e3
        Pg = gerg_bubble(T, xO2)
        Pp = _pr_bubble_pressure_manual(T, xO2, 0.055)
        eg = (Pg - Pexp)/Pexp*100.0
        ep = (Pp - Pexp)/Pexp*100.0
        if math.isfinite(eg):
            err_gerg.append(abs(eg))
        err_pr055.append(abs(ep))
        print(f"{T:7.2f} {xO2:7.4f} {Pexp/1e3:10.1f} "
              f"{float(r['expanded_pressure_uncertainty_kpa']):6.1f} "
              f"{Pg/1e3:10.1f} {Pp/1e3:10.1f} "
              f"{eg:7.2f} {ep:7.2f}")

    print("\nSummary (absolute relative error on bubble pressure):")
    if err_gerg:
        print(f"  GERG-2008 (torch-flash): mean={statistics.fmean(err_gerg):.2f}% "
              f"max={max(err_gerg):.2f}%")
    else:
        print("  GERG-2008: no converged points")
    print(f"  PR kij=0.055: mean={statistics.fmean(err_pr055):.2f}% "
          f"max={max(err_pr055):.2f}%")

    if args.check:
        ok = bool(err_gerg) and statistics.fmean(err_gerg) \
            < statistics.fmean(err_pr055)
        print("GERG_O2C3H8_VLE", "PASS" if ok else "FAIL")
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
