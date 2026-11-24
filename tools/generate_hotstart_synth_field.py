#!/usr/bin/env python3
"""Generate a synthetic RFAM hot-start field CSV (rfam.hotstart-field/v1).

Creates a physically layered field for the smoke sector (60x20x12):
  - injector-adjacent cells: 130 K LOX (GERG liquid) and 300 K C3H8
  - chamber body: 2200 K combustion products (CO2/H2O)
  - walls: 1200 K products
Used to exercise --hotstart-field import in solver_3d_sector.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

L = 0.4065      # chamber length [m]
RMAX = 0.065    # chamber radius [m]
SECTOR = math.radians(120.0)

NX, NR, NT = 60, 20, 12


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(
        "tmp/hotstart_synth_field_v1.csv"))
    args = parser.parse_args()

    header = ["x", "y", "z", "p", "T", "u", "v", "w",
              "O2", "C3H8", "CO2", "H2O"]
    rows = []
    for i in range(NX):
        x = L * (i + 0.5) / NX
        for j in range(NR):
            r = RMAX * (j + 0.5) / NR
            for k in range(NT):
                theta = SECTOR * (k + 0.5) / NT
                y = r * math.cos(theta)
                z = r * math.sin(theta)
                p = 4.0e6
                u, v, w = 0.0, 0.0, 0.0
                if x < 0.020 and r < 0.006:
                    # injector core: 130 K LOX
                    T = 130.0
                    Y = {"O2": 1.0, "C3H8": 0.0, "CO2": 0.0, "H2O": 0.0}
                elif x < 0.020 and r < 0.012:
                    # injector annulus: 300 K C3H8
                    T = 300.0
                    Y = {"O2": 0.0, "C3H8": 1.0, "CO2": 0.0, "H2O": 0.0}
                elif r > 0.055:
                    # wall region: 1200 K products
                    T = 1200.0
                    Y = {"O2": 0.0, "C3H8": 0.0,
                         "CO2": 0.6468, "H2O": 0.3532}
                else:
                    # chamber body: 2200 K products
                    T = 2200.0
                    Y = {"O2": 0.0, "C3H8": 0.0,
                         "CO2": 0.6468, "H2O": 0.3532}
                rows.append([x, y, z, p, T, u, v, w,
                             Y["O2"], Y["C3H8"], Y["CO2"], Y["H2O"]])
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        f.write("# rfam.hotstart-field/v1\n")
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f"wrote {args.output}: {len(rows)} cells "
          f"({NX}x{NR}x{NT})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
