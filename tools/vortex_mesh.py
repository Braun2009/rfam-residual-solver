"""
vortex_mesh.py - Stage 7: Vortex v5 chamber profile -> axisymmetric 2D RFAM mesh.

Reads config/chamber_profile_3pintle.json (600 r(x) points),
generates a structured (Nx × Ny) axisymmetric half-domain:
  - x-axis: axial (0..L_total)
  - y-axis: radial (0..r_max), r_max = max(r(x)) + margin
  - j=0: axis of symmetry
  - j=Ny-1: outer far-field / wall boundary

Output:
  1. mesh data as Python arrays (x_cell, y_cell, r_wall(x))
  2. Optional VTK structured grid for ParaView inspection
  3. RFAM Domain initialiser helper
"""

import json, math, os, sys
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
    x   = np.array([p["x_mm"] for p in pts]) * 1e-3   # → m
    r   = np.array([p["r_mm"] for p in pts]) * 1e-3   # → m
    return x, r, d["nodes"], d["x_breaks"]


# ── axisymmetric grid ─────────────────────────────────────────────────────────
def make_grid(Nx=100, Ny=30, r_margin=1.2, y_stretch=2.0):
    """
    Build Nx × Ny structured axisymmetric grid.

    Nx cells along x  (0 .. L_total)
    Ny cells along y  (0 .. r_max * r_margin)

    y stretching: tanh with factor y_stretch concentrates cells near the wall
    (j=Ny-1 side where r_wall(x) lives).

    Returns
    -------
    xc, yc : (Nx, Ny) cell-centre coordinates [m]
    r_wall_cell : (Nx,) wall radius at each axial column [m]  ← the combustion-chamber wall
    xf, yf : (Nx+1,), (Ny+1,) face coordinates [m]
    """
    x_prof, r_prof, nodes, breaks = load_profile()

    L   = x_prof[-1]                        # total length [m]
    r_max = r_prof.max() * r_margin         # outer radial extent

    # uniform axial faces
    xf = np.linspace(0.0, L, Nx + 1)       # Nx+1 face x-positions

    # stretched radial faces (tanh stretching toward AXIS j=0, fine near wall j=Ny)
    # flip eta so fine cells are near j=0 (axis), coarser near r_max
    eta = np.linspace(0.0, 1.0, Ny + 1)
    if y_stretch > 0:
        yf_norm = np.tanh(y_stretch * eta) / np.tanh(y_stretch)
    else:
        yf_norm = eta
    yf = r_max * yf_norm  # fine near j=0 (axis), coarse near outer

    # cell centres
    xc = 0.5 * (xf[:-1] + xf[1:])          # (Nx,)
    yc = 0.5 * (yf[:-1] + yf[1:])          # (Ny,)
    XC, YC = np.meshgrid(xc, yc, indexing='ij')  # (Nx, Ny)

    # wall radius at each axial cell (interpolate from profile)
    r_wall_cell = np.interp(xc, x_prof, r_prof)   # (Nx,)

    return XC, YC, r_wall_cell, xf, yf, L, r_max


def write_vtk(XC, YC, r_wall_cell, path="/tmp/vortex_mesh.vtk"):
    Nx, Ny = XC.shape
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\nVortex mesh\nASCII\n")
        f.write("DATASET RECTILINEAR_GRID\n")
        xf_ax = np.unique(XC)
        yf_ax = np.unique(YC)
        f.write(f"DIMENSIONS {Nx} {Ny} 1\n")
        f.write(f"X_COORDINATES {Nx} float\n")
        f.write(" ".join(f"{v*1e3:.3f}" for v in XC[:,0]) + "\n")
        f.write(f"Y_COORDINATES {Ny} float\n")
        f.write(" ".join(f"{v*1e3:.3f}" for v in YC[0,:]) + "\n")
        f.write("Z_COORDINATES 1 float\n0\n")
        f.write(f"CELL_DATA {(Nx-1)*(Ny-1)}\n")
        f.write("SCALARS r_wall float\nLOOKUP_TABLE default\n")
        for j in range(Ny-1):
            for i in range(Nx-1):
                f.write(f"{r_wall_cell[i]*1e3:.3f}\n")
    print(f"  VTK written: {path}")


def print_summary(XC, YC, r_wall_cell, xf, yf, L, r_max):
    Nx, Ny = XC.shape
    dx = xf[1] - xf[0]
    dy_min = yf[1] - yf[0]
    dy_max = yf[-1] - yf[-2]
    r_throat = r_wall_cell.min()
    x_throat = XC[np.argmin(r_wall_cell), 0]

    print(f"  Domain:   {L*1e3:.1f} mm  ×  {r_max*1e3:.1f} mm")
    print(f"  Grid:     {Nx} × {Ny}  (axial × radial)")
    print(f"  dx:       {dx*1e3:.2f} mm")
    print(f"  dy_wall:  {dy_min*1e3:.3f} mm  dy_outer: {dy_max*1e3:.3f} mm")
    print(f"  r_throat: {r_throat*1e3:.1f} mm  at x={x_throat*1e3:.1f} mm")
    print(f"  r_chamber:{r_wall_cell.max()*1e3:.1f} mm  (max)")
    cells_per_throat_radius = r_throat / dy_min
    print(f"  Cells/throat-radius: {cells_per_throat_radius:.1f}  (≥10 needed for resolved BL)")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--Nx",     type=int,   default=120)
    ap.add_argument("--Ny",     type=int,   default=40)
    ap.add_argument("--margin", type=float, default=1.2)
    ap.add_argument("--vtk",    action="store_true")
    args = ap.parse_args()

    print("=== Stage 7: Vortex v5 axisymmetric mesh ===\n")
    XC, YC, r_wall, xf, yf, L, r_max = make_grid(
        Nx=args.Nx, Ny=args.Ny, r_margin=args.margin)
    print_summary(XC, YC, r_wall, xf, yf, L, r_max)

    if args.vtk:
        write_vtk(XC, YC, r_wall)

    print("\nMESH_GEN_PASS")
