#!/usr/bin/env python3
"""Statistical-convergence audit for a long RFAM run.

Compares consecutive checkpoint-exported fields (via extract_hotstart_field)
and reports the field difference norm vs step count: a settled flow shows
the L2/Linf field deltas decaying and flattening at small values.

Usage:
  audit_run_convergence.py <chk_prefix...>  (each: <prefix>.chk at iter N)
  e.g. audit_run_convergence.py /root/long_run_500 /root/long_run_1000 ...
Requires the fields to be exported first:
  extract_hotstart_field <chk> <field.csv>
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path


def load_field(path: str):
    with open(path, encoding="utf-8") as f:
        lines = [l for l in f if not l.startswith("#")]
    rows = list(csv.DictReader(lines))
    out = {}
    for r in rows:
        # accept both hotstart-field (x,y,z,T,p,rho) and solver-export
        # (x_m,y_m,z_m,T_K,p_Pa,rho_kg_m3) column names
        xk = "x" if "x" in r else "x_m"
        yk = "y" if "y" in r else "y_m"
        zk = "z" if "z" in r else "z_m"
        tk = "T" if "T" in r else "T_K"
        pk = "p" if "p" in r else "p_Pa"
        rk = "rho" if "rho" in r else "rho_kg_m3"
        key = (float(r[xk]), float(r[yk]), float(r[zk]))
        out[key] = (float(r[tk]), float(r[pk]), float(r[rk]))
    return out


def compare(a: dict, b: dict):
    keys = set(a) & set(b)
    dT = 0.0
    dp = 0.0
    drho = 0.0
    n = len(keys)
    for k in keys:
        dT = max(dT, abs(b[k][0] - a[k][0]))
        dp = max(dp, abs(b[k][1] - a[k][1]) / max(abs(a[k][1]), 1.0))
        drho = max(drho, abs(b[k][2] - a[k][2]) / max(abs(a[k][2]), 1.0))
    return n, dT, dp, drho


def main() -> int:
    fields = [Path(p) for p in sys.argv[1:]]
    if len(fields) < 2:
        print("usage: audit_run_convergence.py field1.csv field2.csv ...")
        return 2
    prev = load_field(str(fields[0]))
    print(f"field 0: {len(prev)} cells")
    for i in range(1, len(fields)):
        cur = load_field(str(fields[i]))
        n, dT, dp, drho = compare(prev, cur)
        print(f"delta[{i-1}->{i}]: cells={n} dT_linf={dT:.4f} K "
              f"dp_rel_linf={dp:.3e} drho_rel_linf={drho:.3e}")
        prev = cur
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
