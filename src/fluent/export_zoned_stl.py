#!/usr/bin/env python3
"""Export a Gmsh surface mesh as a multi-solid ASCII STL keyed by physical group."""
import argparse
import json
import math
import sys
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
    counts = {}
    with args.output.open("w", encoding="ascii", newline="\n") as f:
        for dim, pg in sorted(gmsh.model.getPhysicalGroups(2),
                              key=lambda x: gmsh.model.getPhysicalName(*x)):
            name = gmsh.model.getPhysicalName(dim, pg)
            count = 0
            f.write(f"solid {name}\n")
            for ent in gmsh.model.getEntitiesForPhysicalGroup(dim, pg):
                types, _, conns = gmsh.model.mesh.getElements(dim, int(ent))
                for typ, conn in zip(types, conns):
                    props = gmsh.model.mesh.getElementProperties(typ)
                    nper = props[3]
                    if not props[0].startswith("Triangle"):
                        continue
                    for i in range(0, len(conn), nper):
                        a, b, c = (pts[int(x)] for x in conn[i:i + 3])
                        ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
                        vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
                        nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
                        mag = math.sqrt(nx * nx + ny * ny + nz * nz)
                        nx, ny, nz = nx / mag, ny / mag, nz / mag
                        f.write(f"  facet normal {nx:.9g} {ny:.9g} {nz:.9g}\n    outer loop\n")
                        for p in (a, b, c):
                            f.write(f"      vertex {p[0]:.12g} {p[1]:.12g} {p[2]:.12g}\n")
                        f.write("    endloop\n  endfacet\n")
                        count += 1
            f.write(f"endsolid {name}\n")
            counts[name] = count
    print(json.dumps(counts, sort_keys=True))
    print("TOTAL", sum(counts.values()))
    gmsh.finalize()


if __name__ == "__main__":
    main()
