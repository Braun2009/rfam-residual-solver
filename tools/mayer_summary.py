#!/usr/bin/env python3
"""Mayer-case standard summary (rfam side) for the rfam-vs-Fluent comparison.

Reads solver_mayer_o2 outputs (mayer_centerline.dat, mayer_spreading.dat)
plus the run-time header and writes a frozen JSON summary with the standard
comparison quantities:
  - rho_jet / rho_ch (density ratio)
  - potential core length L_core/D
  - spreading rate d(y_05)/dx
  - centerline density decay table (x/D, rho/rho_jet)
The Fluent side is expected to export the same quantities in the same
schema for a direct diff (task cd8e1c9d).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def load_centerline(path: Path):
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 4:
                rows.append({
                    "x_m": float(parts[0]),
                    "x_D": float(parts[1]),
                    "rho_kgm3": float(parts[2]),
                    "rho_norm": float(parts[3]),
                })
    return rows


def load_spreading(path: Path):
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 4:
                rows.append({
                    "x_m": float(parts[0]),
                    "x_D": float(parts[1]),
                    "y05_m": float(parts[2]),
                    "y05_D": float(parts[3]),
                })
    return rows


def main() -> int:
    cwd = Path.cwd()
    center = cwd / "mayer_centerline.dat"
    spread = cwd / "mayer_spreading.dat"
    if not center.exists():
        print(f"MISSING {center}")
        return 2
    cl = load_centerline(center)
    sp = load_spreading(spread) if spread.exists() else []

    # Potential core: first x where rho_norm < 0.95
    L_core = None
    for r in cl:
        if r["rho_norm"] < 0.95:
            L_core = r["x_D"]
            break
    # Spreading rate: linear fit of y05_D vs x_D
    slope = None
    if len(sp) >= 3:
        xs = [r["x_D"] for r in sp]
        ys = [r["y05_D"] for r in sp]
        n = len(xs)
        sx = sum(xs)
        sy = sum(ys)
        sxx = sum(x * x for x in xs)
        sxy = sum(x * y for x, y in zip(xs, ys))
        denom = n * sxx - sx * sx
        if abs(denom) > 1e-30:
            slope = (n * sxy - sx * sy) / denom

    summary = {
        "schema": "rfam.mayer-summary/v1",
        "source": "solver_mayer_o2 (GERG-2008)",
        "potential_core_length_D": L_core,
        "spreading_rate_dy05_dx": slope,
        "centerline_points": len(cl),
        "spreading_points": len(sp),
        "centerline": cl,
        "spreading": sp,
    }
    out = cwd / "mayer_summary.json"
    out.write_text(json.dumps(summary, indent=1))
    print(f"wrote {out}")
    print(f"  L_core/D = {L_core}")
    print(f"  d(y05)/dx = {slope}")
    print(f"  reference (Mayer 1998): L_core/D ~ 8-12, d(y05)/dx ~ 0.08-0.11")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
