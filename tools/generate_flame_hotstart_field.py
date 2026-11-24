#!/usr/bin/env python3
"""Generate a physically-real hot-start field from a Cantera 1D flame.

A freely-propagating premixed C3H8/air flame (GRI-30, rights-clean) is
computed in 1D; the temperature/species profiles are mapped onto the RFAM
sector: the flame sheet is placed across the chamber, hot products fill
the downstream half, fresh mixture fills the injector-adjacent upstream,
and the near-wall layer carries the wall temperature profile.

Output: rfam.hotstart-field/v1 CSV (x,y,z,p,T,u,v,w,rho,e,Y...) for
solver_3d_sector --hotstart-field.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import cantera as ct

L = 0.4065      # chamber length [m]
RMAX = 0.065    # chamber radius [m]
SECTOR = math.radians(120.0)
P0 = 4.0e6      # chamber pressure [Pa]

NX, NR, NT = 60, 20, 12
SPECIES = ["O2", "N2", "H2", "CH4", "C2H6", "C3H8", "H2O", "CO", "CO2"]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=Path("tmp/hotstart_flame_field_v1.csv"))
    parser.add_argument("--phi", type=float, default=1.0)
    parser.add_argument("--tu", type=float, default=700.0)  # unburnt T [K]
    args = parser.parse_args()

    # 1D free flame at chamber pressure.
    gas = ct.Solution("tests/reference/cantera_gri30_v3_0.yaml")
    gas.TPX = args.tu, P0, {"C3H8": 1.0, "O2": 5.0 / args.phi,
                            "N2": 18.8 / args.phi}
    width = 0.05  # 5 cm domain
    f = ct.FreeFlame(gas, width=width)
    f.set_refine_criteria(ratio=3.0, slope=0.08, curve=0.12)
    f.solve(loglevel=0, auto=True)
    xf = f.grid * 1.0
    Tf = f.T
    # flame position: point of max temperature gradient
    dT = [Tf[i + 1] - Tf[i] for i in range(len(Tf) - 1)]
    x_flame = xf[dT.index(max(dT)) + 1]

    def interp(prof, xq):
        import numpy as np
        return float(np.interp(xq, xf, np.asarray(prof, dtype=float)))

    species_index = {sp: gas.species_index(sp) for sp in SPECIES}
    print(f"flame: phi={args.phi} T_unburnt={args.tu}K "
          f"x_flame={x_flame*1e3:.1f}mm "
          f"T_burnt={Tf[-1]:.0f}K")

    header = ["x", "y", "z", "p", "T", "u", "v", "w", "rho", "e"] + SPECIES
    rows = []
    for i in range(NX):
        x = L * (i + 0.5) / NX
        # upstream distance from flame sheet (negative = unburnt side)
        xrel = x - 0.20  # flame sheet at 200 mm
        for j in range(NR):
            r = RMAX * (j + 0.5) / NR
            for k in range(NT):
                theta = SECTOR * (k + 0.5) / NT
                y = r * math.cos(theta)
                z = r * math.sin(theta)
                # wall thermal layer: linear 700 K -> flame T over 8 mm
                wall_factor = max(0.0, (r - (RMAX - 0.008)) / 0.008)
                # sample the 1D flame at the upstream coordinate
                xs = max(xf[0], min(xf[-1], xrel))
                Tq = interp(Tf, xs)
                Tw = Tq + (700.0 - Tq) * wall_factor
                T = Tw
                # mass fractions from the flame profile (registry species)
                Y = {}
                for sp in SPECIES:
                    if sp in species_index:
                        Y[sp] = max(interp(f.Y[species_index[sp]], xs), 0.0)
                # velocity: burnt-side products move downstream ~ few m/s
                u = 1.0 if xrel > 0.0 else 0.0
                rows.append([x, y, z, P0, T, u, 0.0, 0.0, 0.0, 0.0]
                            + [Y.get(sp, 0.0) for sp in SPECIES])
    with open(args.output, "w", newline="", encoding="utf-8") as fout:
        fout.write("# rfam.hotstart-field/v1\n")
        w = csv.writer(fout)
        w.writerow(header)
        w.writerows(rows)
    print(f"wrote {args.output}: {len(rows)} cells")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
