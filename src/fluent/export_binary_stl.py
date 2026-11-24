#!/usr/bin/env python3
"""Export all triangles of a Gmsh surface mesh into ONE binary STL (view only)."""
import argparse
import struct
from pathlib import Path

import gmsh


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(str(args.input))
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    pts = {int(t): (float(coords[3 * i]), float(coords[3 * i + 1]), float(coords[3 * i + 2]))
           for i, t in enumerate(node_tags)}
    triangles = []
    for dim, tag in gmsh.model.getEntities(2):
        types, _, conns = gmsh.model.mesh.getElements(dim, tag)
        for typ, conn in zip(types, conns):
            props = gmsh.model.mesh.getElementProperties(typ)
            nper = props[3]
            if not props[0].startswith("Triangle"):
                continue
            for i in range(0, len(conn), nper):
                triangles.append([pts[int(x)] for x in conn[i:i + 3]])
    with args.output.open("wb") as f:
        f.write(b"vortex split-segments surface mesh" + b"\0" * (80 - 34))
        f.write(struct.pack("<I", len(triangles)))
        for a, b, c in triangles:
            ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
            vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
            nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
            mag = (nx * nx + ny * ny + nz * nz) ** 0.5
            f.write(struct.pack("<3f", nx / mag, ny / mag, nz / mag))
            for p in (a, b, c):
                f.write(struct.pack("<3f", *p))
            f.write(struct.pack("<H", 0))
    print("TRIANGLES", len(triangles))
    gmsh.finalize()


if __name__ == "__main__":
    main()
