#!/usr/bin/env python3
"""
gen_sector_mesh.py - Generate 120° sector structured mesh for Vortex 3-pintle chamber.

Reads Vortex v5 chamber profile (config/chamber_profile_3pintle.json) and generates
a structured cylindrical (r, theta, x) mesh for RFAM 3D simulation.

Coordinate convention:
  x = axial (0 → L_total = 406.5mm)
  y = radial (0 → r_wall(x))
  z = azimuthal (0 → 2π/3 = 120°)

Output: mesh_info dict + optional VTK for ParaView inspection.

Usage:
  python3 tools/gen_sector_mesh.py --Nx 120 --Nr 30 --Ntheta 20 --level pathfinder
  python3 tools/gen_sector_mesh.py --Nx 400 --Nr 80 --Ntheta 60 --level prod
"""

import json, math, os, sys, argparse
import numpy as np

# ── locate chamber profile ────────────────────────────────────────────────────
def _find_profile():
    candidates = [
        "config/chamber_profile_3pintle.json",
        "/mnt/f/Vortex/config/chamber_profile_3pintle.json",
        os.path.join(os.path.dirname(__file__), "..", "..", "Vortex",
                     "config", "chamber_profile_3pintle.json"),
    ]
    for p in candidates:
        if os.path.exists(p):
            return p
    raise FileNotFoundError("chamber_profile_3pintle.json not found")

def load_profile():
    path = _find_profile()
    with open(path) as f:
        d = json.load(f)
    pts = d["profile"]
    x = np.array([p["x_mm"] for p in pts]) * 1e-3   # m
    r = np.array([p["r_mm"] for p in pts]) * 1e-3   # m
    return x, r, d["nodes"], d["x_breaks"]

# ── sector mesh generation ────────────────────────────────────────────────────
def gen_sector_mesh(Nx=120, Nr=30, Ntheta=20,
                   sector_deg=120.0, r_margin=1.0):
    """
    Generate structured cylindrical sector mesh.

    Returns dict with:
      x_cell[Nx]          axial cell centres [m]
      r_wall[Nx]          wall radius at each axial station [m]
      r_cell[Nr]          radial cell centres (0 → r_wall, normalised)
      theta_cell[Ntheta]  azimuthal cell centres [rad]
      dx, dr_min, dtheta  cell sizes
      r_throat, i_throat  throat geometry
      total_cells         Nx * Nr * Ntheta
    """
    x_prof, r_prof, nodes, breaks = load_profile()

    L = x_prof[-1]
    sector_rad = math.radians(sector_deg)

    # Axial: uniform
    xf = np.linspace(0.0, L, Nx + 1)
    xc = 0.5 * (xf[:-1] + xf[1:])

    # Wall radius at each axial station
    r_wall = np.interp(xc, x_prof, r_prof)

    # Radial: uniform 0 → r_wall(x) per column (normalised 0→1, scaled per column)
    # For the mesh we store normalised radial coordinate eta ∈ [0,1]
    eta_f = np.linspace(0.0, 1.0, Nr + 1)
    eta_c = 0.5 * (eta_f[:-1] + eta_f[1:])

    # Azimuthal: uniform 0 → sector_rad
    theta_f = np.linspace(0.0, sector_rad, Ntheta + 1)
    theta_c = 0.5 * (theta_f[:-1] + theta_f[1:])

    # Key dimensions
    r_throat = r_wall.min()
    i_throat = int(np.argmin(r_wall))
    r_chamber = r_wall.max()

    dx = L / Nx
    dr_min = r_throat / Nr   # smallest radial cell (at throat)
    dtheta = sector_rad / Ntheta

    info = {
        "Nx": Nx, "Nr": Nr, "Ntheta": Ntheta,
        "total_cells": Nx * Nr * Ntheta,
        "L": L, "sector_deg": sector_deg,
        "x_cell": xc,
        "r_wall": r_wall,
        "eta_cell": eta_c,   # normalised radial [0,1]
        "theta_cell": theta_c,
        "dx": dx, "dr_min": dr_min, "dtheta": dtheta,
        "r_throat": r_throat, "i_throat": i_throat,
        "r_chamber": r_chamber,
    }
    return info

def print_summary(info):
    print(f"  Sector mesh: {info['Nx']}×{info['Nr']}×{info['Ntheta']} "
          f"= {info['total_cells']:,} cells")
    print(f"  Domain: L={info['L']*1e3:.1f}mm, sector={info['sector_deg']}°")
    print(f"  dx={info['dx']*1e3:.2f}mm  dr_min={info['dr_min']*1e3:.3f}mm  "
          f"dtheta={math.degrees(info['dtheta']):.2f}°")
    print(f"  r_throat={info['r_throat']*1e3:.1f}mm @ i={info['i_throat']}  "
          f"r_chamber={info['r_chamber']*1e3:.1f}mm")
    print(f"  Cells/throat-radius (radial): {info['r_throat']/info['dr_min']:.0f}")

def write_vtk(info, path="/tmp/sector_mesh.vtk"):
    """Write VTK structured grid for ParaView inspection."""
    Nx, Nr, Nt = info["Nx"], info["Nr"], info["Ntheta"]
    xc, r_wall, eta_c, theta_c = (info["x_cell"], info["r_wall"],
                                   info["eta_cell"], info["theta_cell"])
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\nSector mesh\nASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        n_pts = Nx * Nr * Nt
        f.write(f"POINTS {n_pts} float\n")
        for i in range(Nx):
            for j in range(Nr):
                r = eta_c[j] * r_wall[i]
                for k in range(Nt):
                    theta = theta_c[k]
                    # cylindrical → Cartesian: y=r*cos(theta), z=r*sin(theta)
                    f.write(f"{xc[i]*1e3:.3f} {r*math.cos(theta)*1e3:.3f} "
                            f"{r*math.sin(theta)*1e3:.3f}\n")
    print(f"  VTK written: {path}")

def write_rfam_init(info, path="/tmp/sector_ic.csv"):
    """Write initial condition CSV for RFAM Domain initialisation."""
    with open(path, "w") as f:
        f.write("# i,j,k,x[m],r[m],theta[rad],r_wall[m]\n")
        for i in range(info["Nx"]):
            for j in range(info["Nr"]):
                r = info["eta_cell"][j] * info["r_wall"][i]
                for k in range(info["Ntheta"]):
                    f.write(f"{i},{j},{k},"
                            f"{info['x_cell'][i]:.6f},"
                            f"{r:.6f},"
                            f"{info['theta_cell'][k]:.6f},"
                            f"{info['r_wall'][i]:.6f}\n")
    print(f"  IC CSV written: {path}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--Nx",     type=int,   default=120)
    ap.add_argument("--Nr",     type=int,   default=30)
    ap.add_argument("--Ntheta", type=int,   default=20)
    ap.add_argument("--sector", type=float, default=120.0)
    ap.add_argument("--level",  type=str,   default="pathfinder",
                    choices=["pathfinder","prod"])
    ap.add_argument("--vtk",    action="store_true")
    ap.add_argument("--csv",    action="store_true")
    args = ap.parse_args()

    if args.level == "prod":
        args.Nx, args.Nr, args.Ntheta = 400, 80, 60

    print(f"=== Vortex 3-pintle sector mesh ({args.level}) ===\n")
    info = gen_sector_mesh(args.Nx, args.Nr, args.Ntheta, args.sector)
    print_summary(info)

    if args.vtk:
        write_vtk(info)
    if args.csv:
        write_rfam_init(info)

    print("\nMESH_GEN_PASS")
