#!/usr/bin/env python3
"""Plume-domain Vortex geometry: open farfield downstream of the nozzle exit.

This is the pressure-based outer segment of the decoupled two-segment
workflow.  The domain starts at the nozzle exit plane (x = 406.5 mm) and
extends downstream to x = 1406.5 mm inside a cylindrical outer boundary of
radius 300 mm.  The nozzle exit disk (r <= 60 mm) is the profile inlet
`inlet_nozzle_exit`; the annular disk region r in [60, 300] mm at the same
plane stands in for the nozzle aft-body surface (`wall_nozzle_aft`) until
the remodelled outer wall geometry replaces it.

Fluid: one continuous meridian face revolved around the x axis.  Units: mm.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import gmsh


CHAMBER_END_MM = 406.5
PLUME_END_MM = 1406.5
EXIT_RADIUS_MM = 60.0
FARFIELD_RADIUS_MM = 300.0

PHYSICAL_IDS = {
    "fluid": 1,
    "inlet_nozzle_exit": 21,
    "wall_nozzle_aft": 22,
    "farfield": 25,
}

PRESETS = {
    "coarse": {
        "global_min": 2.0,
        "global_max": 60.0,
        "inlet": 3.0,
        "aft_wall": 8.0,
        "core": 10.0,
        "farfield": 50.0,
    },
    "production": {
        "global_min": 0.5,
        "global_max": 45.0,
        "inlet": 1.0,
        "aft_wall": 4.0,
        "core": 2.0,
        "farfield": 40.0,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    parser.add_argument("--mode", choices=("geometry", "surface", "volume"), default="surface")
    parser.add_argument("--preset", choices=tuple(PRESETS), default="production")
    parser.add_argument("--terminal", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def build_geometry() -> int:
    occ = gmsh.model.occ
    p_axis_in = occ.addPoint(CHAMBER_END_MM, 0.0, 0.0)
    p_exit = occ.addPoint(CHAMBER_END_MM, 0.0, EXIT_RADIUS_MM)
    p_far_up = occ.addPoint(CHAMBER_END_MM, 0.0, FARFIELD_RADIUS_MM)
    p_far_end = occ.addPoint(PLUME_END_MM, 0.0, FARFIELD_RADIUS_MM)
    p_axis_end = occ.addPoint(PLUME_END_MM, 0.0, 0.0)
    curves = [
        occ.addLine(p_axis_in, p_exit),      # nozzle exit disk radius
        occ.addLine(p_exit, p_far_up),       # aft-body stand-in annulus
        occ.addLine(p_far_up, p_far_end),    # cylindrical farfield
        occ.addLine(p_far_end, p_axis_end),  # downstream cap
        occ.addLine(p_axis_end, p_axis_in),  # axis
    ]
    meridian = occ.addPlaneSurface([occ.addCurveLoop(curves)])
    revolved = occ.revolve([(2, meridian)], 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0 * math.pi)
    volumes = [tag for dim, tag in revolved if dim == 3]
    occ.synchronize()
    if len(volumes) != 1 or len(gmsh.model.getEntities(3)) != 1:
        raise RuntimeError("expected exactly one revolved volume")
    volume = volumes[0]
    boundary = {tag for dim, tag in gmsh.model.getBoundary([(3, volume)], oriented=False, recursive=False) if dim == 2}
    orphans = [(dim, tag) for dim, tag in gmsh.model.getEntities(2) if tag not in boundary]
    if orphans != [(2, meridian)]:
        raise RuntimeError(f"unexpected orphan surfaces: {orphans}")
    occ.remove(orphans, recursive=False)
    occ.synchronize()
    return volume


def classify_boundaries(volume: int) -> dict[str, list[int]]:
    boundaries = gmsh.model.getBoundary([(3, volume)], oriented=False, recursive=False)
    surface_tags = [tag for dim, tag in boundaries if dim == 2]
    groups = {name: [] for name in PHYSICAL_IDS if name != "fluid"}
    tol = 1.0e-4
    for tag in surface_tags:
        box = gmsh.model.getBoundingBox(2, tag)
        xmin, ymin, zmin, xmax, ymax, zmax = box
        max_radius = max(abs(ymin), abs(ymax), abs(zmin), abs(zmax))
        if abs(xmin - CHAMBER_END_MM) < tol and abs(xmax - CHAMBER_END_MM) < tol:
            name = "inlet_nozzle_exit" if abs(max_radius - EXIT_RADIUS_MM) < tol else "wall_nozzle_aft"
        else:
            name = "farfield"
        groups[name].append(tag)
    assigned = [t for tags in groups.values() for t in tags]
    if sorted(assigned) != sorted(surface_tags):
        raise RuntimeError("boundary classification incomplete/overlapping")
    expected = {"inlet_nozzle_exit": 1, "wall_nozzle_aft": 1, "farfield": 2}
    counts = {name: len(tags) for name, tags in groups.items()}
    if counts != expected:
        raise RuntimeError(f"unexpected boundary topology: {counts}, expected {expected}")
    return groups


def add_physical_groups(volume: int, groups: dict[str, list[int]]) -> None:
    gmsh.model.addPhysicalGroup(3, [volume], PHYSICAL_IDS["fluid"])
    gmsh.model.setPhysicalName(3, PHYSICAL_IDS["fluid"], "fluid")
    for name, tags in groups.items():
        gmsh.model.addPhysicalGroup(2, tags, PHYSICAL_IDS[name])
        gmsh.model.setPhysicalName(2, PHYSICAL_IDS[name], name)


def add_threshold_field(field_tag, faces, size_min, size_max, distance_max):
    distance = field_tag
    threshold = field_tag + 1
    gmsh.model.mesh.field.add("Distance", distance)
    gmsh.model.mesh.field.setNumbers(distance, "FacesList", faces)
    gmsh.model.mesh.field.setNumber(distance, "Sampling", 200)
    gmsh.model.mesh.field.add("Threshold", threshold)
    gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMin", size_min)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMax", size_max)
    gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(threshold, "DistMax", distance_max)
    return threshold


def configure_mesh(groups: dict[str, list[int]], preset_name: str) -> None:
    preset = PRESETS[preset_name]
    gmsh.option.setNumber("Mesh.MeshSizeMin", preset["global_min"])
    gmsh.option.setNumber("Mesh.MeshSizeMax", preset["global_max"])
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)

    fields = [
        add_threshold_field(1, groups["inlet_nozzle_exit"], preset["inlet"], preset["farfield"], 60.0),
        add_threshold_field(3, groups["wall_nozzle_aft"], preset["aft_wall"], preset["farfield"], 60.0),
    ]
    # Plume core box: high velocity gradient region downstream of the exit.
    gmsh.model.mesh.field.add("Box", 5)
    for key, value in {"VIn": preset["core"], "VOut": preset["farfield"],
                       "XMin": 400.0, "XMax": 620.0, "YMin": -85.0, "YMax": 85.0,
                       "ZMin": -85.0, "ZMax": 85.0, "Thickness": 30.0}.items():
        gmsh.model.mesh.field.setNumber(5, key, value)
    fields.append(5)
    gmsh.model.mesh.field.add("Min", 7)
    gmsh.model.mesh.field.setNumbers(7, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(7)


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", int(args.terminal))
        gmsh.model.add("vortex_plume_406p5")
        volume = build_geometry()
        groups = classify_boundaries(volume)
        add_physical_groups(volume, groups)
        topology = {"volume_mm3": gmsh.model.occ.getMass(3, volume),
                    "boundary_surfaces": groups}
        print("TOPOLOGY", json.dumps(topology, sort_keys=True))
        gmsh.write(str(output_dir / "vortex_plume_406p5.step"))
        gmsh.write(str(output_dir / "vortex_plume_406p5.brep"))
        if args.mode == "geometry":
            print("VORTEX_PLUME_GEOMETRY_OK")
            return
        configure_mesh(groups, args.preset)
        gmsh.model.mesh.generate(2)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.write(str(output_dir / "vortex_plume_406p5_surface.msh"))
        print("VORTEX_PLUME_SURFACE_MESH_OK")
        if args.mode == "volume":
            gmsh.model.mesh.generate(3)
            gmsh.write(str(output_dir / "vortex_plume_406p5_volume.msh"))
            print("VORTEX_PLUME_VOLUME_MESH_OK")
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    main()
