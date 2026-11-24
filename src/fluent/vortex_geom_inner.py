#!/usr/bin/env python3
"""Inner-domain Vortex geometry: combustion chamber + nozzle, x in [0, 406.5] mm.

This is the density-based inner segment of the decoupled two-segment
workflow (chamber/nozzle -> profile -> plume).  The former single-domain
mesh (chamber + nozzle + one-metre farfield) is NOT reused: the farfield
is removed and the nozzle exit plane (x = 406.5, r_e = 60 mm) becomes an
explicit physical boundary `outlet_nozzle`.

Fluid: one continuous meridian face revolved around the x axis.  No pintle
solids, no internal interfaces, no orphan sheets.  Units: mm.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import gmsh


CHAMBER_END_MM = 406.5
TORCH_RADIUS_MM = 2.0
SPRAY_RADIUS_MM = 7.5
PLATE_RADIUS_MM = 24.0
EXIT_RADIUS_MM = 60.0

PHYSICAL_IDS = {
    "fluid": 1,
    "inlet_torch": 11,
    "inlet_spray": 12,
    "wall_plate": 13,
    "wall_chamber": 14,
    "outlet_nozzle": 16,
}

PRESETS = {
    "coarse": {
        "global_min": 1.5,
        "global_max": 50.0,
        "inlet": 1.5,
        "wall": 4.0,
        "wall_far": 30.0,
        "wall_transition": 45.0,
        "throat": 3.0,
        "exit": 8.0,
    },
    "production": {
        "global_min": 0.20,
        "global_max": 20.0,
        "inlet": 0.30,
        "wall": 0.50,
        "wall_far": 8.0,
        "wall_transition": 8.0,
        "throat": 0.50,
        "exit": 1.0,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path,
                        help="chamber_profile_3pintle.json (auto-detected when omitted)")
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    parser.add_argument("--mode", choices=("geometry", "surface", "volume"), default="surface")
    parser.add_argument("--preset", choices=tuple(PRESETS), default="production")
    parser.add_argument("--terminal", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def resolve_profile(requested: Path | None) -> Path:
    candidates = [
        requested,
        Path(os.environ["VORTEX_PROFILE"]) if "VORTEX_PROFILE" in os.environ else None,
        Path("/root/vortex_fluent/chamber_profile_3pintle.json"),
        Path("F:/Vortex/config/chamber_profile_3pintle.json"),
        Path.cwd() / "chamber_profile_3pintle.json",
    ]
    for candidate in candidates:
        if candidate is not None and candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError("chamber profile not found")


def profile_segments(profile_path: Path) -> list[list[tuple[float, float]]]:
    with profile_path.open(encoding="utf-8") as stream:
        document = json.load(stream)
    raw = [(float(p["x_mm"]), float(p["r_mm"])) for p in document["profile"]]
    breaks = [0.0, 123.0, 273.0, 276.5, CHAMBER_END_MM]
    radii = [PLATE_RADIUS_MM, 65.0, 23.6, 23.6, EXIT_RADIUS_MM]
    if len(raw) < 5 or any(raw[i][0] >= raw[i + 1][0] for i in range(len(raw) - 1)):
        raise ValueError("profile x coordinates must be strictly increasing")
    if abs(raw[0][0]) > 1.0e-9 or abs(raw[-1][0] - CHAMBER_END_MM) > 1.0e-6:
        raise ValueError("profile must span x=0..406.5 mm")
    segments: list[list[tuple[float, float]]] = []
    for index, (x0, x1) in enumerate(zip(breaks[:-1], breaks[1:])):
        inside = [(x, r) for x, r in raw if x0 < x < x1]
        stride = max(1, math.ceil(len(inside) / 40))
        sampled = inside[::stride]
        segments.append([(x0, radii[index]), *sampled, (x1, radii[index + 1])])
    return segments


def add_spline(occ, samples, first_point=None):
    tags = []
    for index, (x, radius) in enumerate(samples):
        if index == 0 and first_point is not None:
            tags.append(first_point)
        else:
            tags.append(occ.addPoint(x, 0.0, radius))
    return occ.addSpline(tags), tags[-1]


def build_geometry(profile_path: Path) -> int:
    occ = gmsh.model.occ
    segments = profile_segments(profile_path)
    p0 = occ.addPoint(0.0, 0.0, 0.0)
    p_torch = occ.addPoint(0.0, 0.0, TORCH_RADIUS_MM)
    p_spray = occ.addPoint(0.0, 0.0, SPRAY_RADIUS_MM)
    p_plate = occ.addPoint(0.0, 0.0, PLATE_RADIUS_MM)
    curves = [
        occ.addLine(p0, p_torch),
        occ.addLine(p_torch, p_spray),
        occ.addLine(p_spray, p_plate),
    ]
    previous = p_plate
    for segment in segments:
        spline, previous = add_spline(occ, segment, previous)
        curves.append(spline)
    # Nozzle exit plane: radial line at x = CHAMBER_END_MM from the wall
    # endpoint (which is exactly (406.5, 0, r_e)) down to the axis.
    p_exit_axis = occ.addPoint(CHAMBER_END_MM, 0.0, 0.0)
    curves.extend([
        occ.addLine(previous, p_exit_axis),
        occ.addLine(p_exit_axis, p0),
    ])
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


def radial_extent(box) -> float:
    return max(abs(box[1]), abs(box[2]), abs(box[4]), abs(box[5]))


def classify_boundaries(volume: int) -> dict[str, list[int]]:
    boundaries = gmsh.model.getBoundary([(3, volume)], oriented=False, recursive=False)
    surface_tags = [tag for dim, tag in boundaries if dim == 2]
    groups = {name: [] for name in PHYSICAL_IDS if name != "fluid"}
    tol = 1.0e-4
    for tag in surface_tags:
        box = gmsh.model.getBoundingBox(2, tag)
        xmin, _, _, xmax, _, _ = box
        radius = radial_extent(box)
        if abs(xmin) < tol and abs(xmax) < tol:
            radii = {"inlet_torch": TORCH_RADIUS_MM,
                     "inlet_spray": SPRAY_RADIUS_MM,
                     "wall_plate": PLATE_RADIUS_MM}
            name = min(radii, key=lambda key: abs(radius - radii[key]))
            if abs(radius - radii[name]) > tol:
                raise RuntimeError(f"unrecognized x=0 surface {tag}: bbox={box}")
        elif abs(xmin - CHAMBER_END_MM) < tol and abs(xmax - CHAMBER_END_MM) < tol:
            name = "outlet_nozzle"
        else:
            name = "wall_chamber"
        groups[name].append(tag)
    assigned = [t for tags in groups.values() for t in tags]
    if sorted(assigned) != sorted(surface_tags):
        raise RuntimeError("boundary classification incomplete/overlapping")
    expected = {"inlet_torch": 1, "inlet_spray": 1, "wall_plate": 1,
                "wall_chamber": 4, "outlet_nozzle": 1}
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
        add_threshold_field(1, groups["wall_chamber"], preset["wall"], preset["wall_far"], preset["wall_transition"]),
        add_threshold_field(3, groups["inlet_torch"] + groups["inlet_spray"] + groups["wall_plate"],
                            preset["inlet"], preset["wall_far"], preset["wall_transition"]),
    ]
    gmsh.model.mesh.field.add("Box", 5)
    for key, value in {"VIn": preset["throat"], "VOut": preset["global_max"],
                       "XMin": 258.0, "XMax": 292.0, "YMin": -32.0, "YMax": 32.0,
                       "ZMin": -32.0, "ZMax": 32.0, "Thickness": 8.0}.items():
        gmsh.model.mesh.field.setNumber(5, key, value)
    fields.append(5)
    gmsh.model.mesh.field.add("Box", 6)
    for key, value in {"VIn": preset["exit"], "VOut": preset["global_max"],
                       "XMin": 360.0, "XMax": 406.5, "YMin": -70.0, "YMax": 70.0,
                       "ZMin": -70.0, "ZMax": 70.0, "Thickness": 20.0}.items():
        gmsh.model.mesh.field.setNumber(6, key, value)
    fields.append(6)
    gmsh.model.mesh.field.add("Min", 7)
    gmsh.model.mesh.field.setNumbers(7, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(7)


def main() -> None:
    args = parse_args()
    profile_path = resolve_profile(args.profile)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", int(args.terminal))
        gmsh.model.add("vortex_inner_406p5")
        volume = build_geometry(profile_path)
        groups = classify_boundaries(volume)
        add_physical_groups(volume, groups)
        topology = {"profile": str(profile_path),
                    "volume_mm3": gmsh.model.occ.getMass(3, volume),
                    "boundary_surfaces": groups}
        print("TOPOLOGY", json.dumps(topology, sort_keys=True))
        gmsh.write(str(output_dir / "vortex_inner_406p5.step"))
        gmsh.write(str(output_dir / "vortex_inner_406p5.brep"))
        if args.mode == "geometry":
            print("VORTEX_INNER_GEOMETRY_OK")
            return
        configure_mesh(groups, args.preset)
        gmsh.model.mesh.generate(2)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.write(str(output_dir / "vortex_inner_406p5_surface.msh"))
        print("VORTEX_INNER_SURFACE_MESH_OK")
        if args.mode == "volume":
            gmsh.model.mesh.generate(3)
            gmsh.write(str(output_dir / "vortex_inner_406p5_volume.msh"))
            print("VORTEX_INNER_VOLUME_MESH_OK")
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    main()
