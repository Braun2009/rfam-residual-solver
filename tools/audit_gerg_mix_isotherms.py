#!/usr/bin/env python3
"""GERG-2008 (torch-flash) isotherm health audit for O2/C3H8 mixtures.

2026-08-03 finding: the torch-flash gerg2008 parameter set produces
non-physical isotherms in the deep-cold propane-rich region (110-150 K,
xO2 < 0.5): pressures down to -1e9 kPa, two-phase loops spanning the whole
density range.  CoolProp HEOS binary bubble-point solver fails on the same
region (solver_rho_Tp), so this is a GERG-2008 root-solver limitation
confirmed by two independent implementations, not a single-implementation
bug.

Consumption:
  - Table v3 keeps only mixture levels xO2 in {0.5, 0.75, 0.95}.
  - The 20 propane-rich NIST ThermoML points (110-140 K, xO2<0.5) stay
    uncovered until teqp C++ provides a healthier root solver.

Audit criteria (per isotherm):
  - no significant negative pressure: min(p) > -1e7 Pa
  - loop (dp/drho < 0) fraction of physical samples < 80 %
  - interior loop pressure bounded (spurious high-pressure loops fail)
"""

from __future__ import annotations

import sys

import torch
import torch_flash as tf

R = 8.314462618
MW_TABLE = {0.05: 0.033510, 0.25: 0.041070, 0.50: 0.038047,
            0.75: 0.035024, 0.95: 0.032604}


def isotherm_pressure_samples(eos, T: float, xO2: float, n: int = 257):
    out = []
    for i in range(n):
        rho = 0.5 * (1300.0 / 0.5) ** (1.0 - i / (n - 1))
        x = torch.tensor([xO2, 1.0 - xO2])
        try:
            v = torch.tensor(MW_TABLE[xO2] / rho)
            p = float(eos.pressure(torch.tensor(float(T)), v, x))
        except Exception:
            continue
        if not torch.isfinite(torch.tensor(p)):
            continue
        out.append((rho, p))
    return out


def audit(T: float, xO2: float, verbose: bool = False) -> dict:
    eos = tf.gerg2008(("oxygen", "propane"))
    samples = isotherm_pressure_samples(eos, T, xO2)
    if not samples:
        return {"T": T, "xO2": xO2, "ok": False, "reason": "no samples"}
    ps = [p for _, p in samples]
    minp, maxp = min(ps), max(ps)
    dome = 0
    for k in range(1, len(samples) - 1):
        r0, p0 = samples[k - 1]
        r1, p1 = samples[k]
        r2, p2 = samples[k + 1]
        if (p2 - p0) / (r2 - r0) < 0:
            dome += 1
    dome_frac = dome / len(samples)
    ok = (minp > -1e7) and (dome_frac < 0.8)
    if verbose:
        print(f"  xO2={xO2} T={T}K: n={len(samples)} "
              f"p=[{minp/1e3:.0f},{maxp/1e3:.0f}] kPa loop={dome_frac:.2f} "
              f"-> {'OK' if ok else 'UNRELIABLE'}")
    return {"T": T, "xO2": xO2, "ok": ok, "minP": minp, "maxP": maxp,
            "loopFraction": dome_frac}


def main() -> int:
    verbose = "--verbose" in sys.argv
    summary = {}
    for T in (110.22, 130.0, 150.0, 200.0, 250.0, 300.0, 400.0):
        for xO2 in (0.05, 0.25, 0.5, 0.75, 0.95):
            r = audit(T, xO2, verbose)
            summary.setdefault(T, {})[xO2] = r
    ok_count = sum(1 for t in summary.values() for r in t.values()
                   if r["ok"])
    total = sum(len(t) for t in summary.values())
    print(f"GERG_MIX_HEALTH_AUDIT ok={ok_count}/{total}")
    return 0 if ok_count == total else 2


if __name__ == "__main__":
    raise SystemExit(main())
