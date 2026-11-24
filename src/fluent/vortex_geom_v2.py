#!/usr/bin/env python3
"""Generate the corrected full-360 Vortex chamber and plume mesh.

The fluid is built from one continuous meridian face and one full revolution.
There are no pintle solids, chamber/plume fusions, or internal CAD interfaces.
This removes the segment/facet intersections that made the previous PLC fail.

Coordinates and mesh sizes are in mm. The axial direction is +x.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
from typing import Iterable

import gmsh


CHAMBER_END_MM = 406.5
PLUME_END_MM = 1406.5
FARFIELD_RADIUS_MM = 300.0
TORCH_RADIUS_MM = 2.0
SPRAY_RADIUS_MM = 7.5
PLATE_RADIUS_MM = 24.0

PHYSICAL_IDS = {
    "fluid": 1,
    "inlet_torch": 11,
    "inlet_spray": 12,
    "wall_plate": 13,
    "wall_chamber": 14,
    "farfield": 15,
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
        "shock": 8.0,
        "farfield": 40.0,
    },
    "production": {
        "global_min": 0.20,
        "global_max": 45.0,
        "inlet": 0.30,
        "wall": 0.50,
        "wall_far": 10.0,
        "wall_transition": 10.0,
        "throat": 0.50,
        "shock": 1.0,
        "farfield": 45.0,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile",
        type=Path,
        help="chamber_profile_3pintle.json (auto-detected when omitted)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path.cwd(),
        help="output directory (default: current directory)",
    )
    parser.add_argument(
        "--mode",
        choices=("geometry", "surface", "volume"),
        default="volume",
        help="highest generation stage (default: volume)",
    )
    parser.add_argument(
        "--preset",
        choices=tuple(PRESETS),
        default="coarse",
        help="mesh sizing preset (default: coarse)",
    )
    parser.add_argument(
        "--terminal",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="show Gmsh terminal output",
    )
    return parser.parse_args()


def resolve_profile(requested: Path | None) -> Path:
    script = Path(__file__).resolve()
    repo = script.parents[2]
    candidates = [
        requested,
        Path(os.environ["VORTEX_PROFILE"]) if "VORTEX_PROFILE" in os.environ else None,
        repo.parent / "Vortex" / "config" / "chamber_profile_3pintle.json",
        Path("F:/Vortex/config/chamber_profile_3pintle.json"),
        Path("/root/vortex_fluent/chamber_profile_3pintle.json"),
        Path.cwd() / "chamber_profile_3pintle.json",
    ]
    for candidate in candidates:
        if candidate is not None and candidate.is_file():
            return candidate.resolve()
    checked = "\n  ".join(str(path) for path in candidates if path is not None)
    raise FileNotFoundError(f"chamber profile not found; checked:\n  {checked}")


def profile_segments(profile_path: Path) -> list[list[tuple[float, float]]]:
    with profile_path.open(encoding="utf-8") as stream:
        document = json.load(stream)
    raw = [(float(point["x_mm"]), float(point["r_mm"])) for point in document["profile"]]
    breaks_doc = document.get("x_breaks", {})
    breaks = [
        float(breaks_doc.get("x0", 0.0)),
        float(breaks_doc.get("x1", 123.0)),
        float(breaks_doc.get("x2", 273.0)),
        float(breaks_doc.get("x3", 276.5)),
        float(breaks_doc.get("x4", CHAMBER_END_MM)),
    ]
    radii = [PLATE_RADIUS_MM, 65.0, 23.6, 23.6, 60.0]

    if len(raw) < 5 or any(raw[index][0] >= raw[index + 1][0] for index in range(len(raw) - 1)):
        raise ValueError("profile x coordinates must be strictly increasing")
    if abs(raw[0][0]) > 1.0e-9 or abs(raw[-1][0] - CHAMBER_END_MM) > 1.0e-6:
        raise ValueError("profile must span x=0..406.5 mm")

    segments: list[list[tuple[float, float]]] = []
    for segment_index, (x0, x1) in enumerate(zip(breaks[:-1], breaks[1:])):
        inside = [(x, radius) for x, radius in raw if x0 < x < x1]
        # About 40 interpolation points per smooth design segment preserve the
        # source contour without creating hundreds of CAD faces.
        stride = max(1, math.ceil(len(inside) / 40))
        sampled = inside[::stride]
        segment = [(x0, radii[segment_index]), *sampled, (x1, radii[segment_index + 1])]
        segments.append(segment)
    return segments


def add_spline(
    occ,
    samples: Iterable[tuple[float, float]],
    first_point: int | None = None,
) -> tuple[int, int]:
    point_tags = []
    for index, (x, radius) in enumerate(samples):
        if index == 0 and first_point is not None:
            point_tags.append(first_point)
        else:
            point_tags.append(occ.addPoint(x, 0.0, radius))
    return occ.addSpline(point_tags), point_tags[-1]


def build_geometry(profile_path: Path) -> int:
    occ = gmsh.model.occ
    segments = profile_segments(profile_path)

    p_axis_start = occ.addPoint(0.0, 0.0, 0.0)
    p_torch = occ.addPoint(0.0, 0.0, TORCH_RADIUS_MM)
    p_spray = occ.addPoint(0.0, 0.0, SPRAY_RADIUS_MM)
    p_plate = occ.addPoint(0.0, 0.0, PLATE_RADIUS_MM)

    curves = [
        occ.addLine(p_axis_start, p_torch),
        occ.addLine(p_torch, p_spray),
        occ.addLine(p_spray, p_plate),
    ]

    previous = p_plate
    for segment in segments:
        spline, previous = add_spline(occ, segment, previous)
        curves.append(spline)

    p_far_up = occ.addPoint(CHAMBER_END_MM, 0.0, FARFIELD_RADIUS_MM)
    p_far_end = occ.addPoint(PLUME_END_MM, 0.0, FARFIELD_RADIUS_MM)
    p_axis_end = occ.addPoint(PLUME_END_MM, 0.0, 0.0)
    curves.extend(
        [
            occ.addLine(previous, p_far_up),
            occ.addLine(p_far_up, p_far_end),
            occ.addLine(p_far_end, p_axis_end),
            occ.addLine(p_axis_end, p_axis_start),
        ]
    )

    loop = occ.addCurveLoop(curves)
    meridian = occ.addPlaneSurface([loop])
    revolved = occ.revolve(
        [(2, meridian)],
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        2.0 * math.pi,
    )
    volumes = [tag for dimension, tag in revolved if dimension == 3]
    occ.synchronize()
    all_volumes = gmsh.model.getEntities(3)
    if len(volumes) != 1 or len(all_volumes) != 1:
        raise RuntimeError(f"expected one revolved volume, got {all_volumes}")

    volume = volumes[0]
    boundary_surfaces = {
        tag
        for dimension, tag in gmsh.model.getBoundary(
            [(3, volume)], oriented=False, recursive=False
        )
        if dimension == 2
    }
    orphan_surfaces = [
        (dimension, tag)
        for dimension, tag in gmsh.model.getEntities(2)
        if tag not in boundary_surfaces
    ]
    if orphan_surfaces != [(2, meridian)]:
        raise RuntimeError(f"unexpected orphan surfaces after revolution: {orphan_surfaces}")
    occ.remove(orphan_surfaces, recursive=False)
    occ.synchronize()
    remaining_surfaces = {tag for _, tag in gmsh.model.getEntities(2)}
    if remaining_surfaces != boundary_surfaces:
        raise RuntimeError("export geometry contains non-boundary sheet bodies")
    return volume


def radial_extent(box: tuple[float, float, float, float, float, float]) -> float:
    return max(abs(box[1]), abs(box[2]), abs(box[4]), abs(box[5]))


def classify_boundaries(volume: int) -> dict[str, list[int]]:
    boundaries = gmsh.model.getBoundary([(3, volume)], oriented=False, recursive=False)
    surface_tags = [tag for dimension, tag in boundaries if dimension == 2]
    groups = {name: [] for name in PHYSICAL_IDS if name != "fluid"}
    tolerance = 1.0e-4

    for tag in surface_tags:
        box = gmsh.model.getBoundingBox(2, tag)
        xmin, _, _, xmax, _, _ = box
        radius = radial_extent(box)
        if abs(xmin) < tolerance and abs(xmax) < tolerance:
            inlet_radii = {
                "inlet_torch": TORCH_RADIUS_MM,
                "inlet_spray": SPRAY_RADIUS_MM,
                "wall_plate": PLATE_RADIUS_MM,
            }
            name = min(inlet_radii, key=lambda key: abs(radius - inlet_radii[key]))
            if abs(radius - inlet_radii[name]) > tolerance:
                raise RuntimeError(f"unrecognized x=0 surface {tag}: bbox={box}")
        elif (
            abs(xmin - CHAMBER_END_MM) < tolerance
            and abs(xmax - CHAMBER_END_MM) < tolerance
        ):
            name = "farfield"
        elif xmax <= CHAMBER_END_MM + tolerance:
            name = "wall_chamber"
        else:
            name = "farfield"
        groups[name].append(tag)

    assigned = [tag for tags in groups.values() for tag in tags]
    if sorted(assigned) != sorted(surface_tags) or len(assigned) != len(set(assigned)):
        raise RuntimeError("boundary classification is incomplete or overlapping")
    expected_counts = {
        "inlet_torch": 1,
        "inlet_spray": 1,
        "wall_plate": 1,
        "wall_chamber": 4,
        "farfield": 3,
    }
    counts = {name: len(tags) for name, tags in groups.items()}
    if counts != expected_counts:
        raise RuntimeError(f"unexpected boundary topology: {counts}, expected {expected_counts}")
    return groups


def add_physical_groups(volume: int, groups: dict[str, list[int]]) -> None:
    gmsh.model.addPhysicalGroup(3, [volume], PHYSICAL_IDS["fluid"])
    gmsh.model.setPhysicalName(3, PHYSICAL_IDS["fluid"], "fluid")
    for name, tags in groups.items():
        gmsh.model.addPhysicalGroup(2, tags, PHYSICAL_IDS[name])
        gmsh.model.setPhysicalName(2, PHYSICAL_IDS[name], name)


def add_threshold_field(
    field_tag: int,
    faces: list[int],
    size_min: float,
    size_max: float,
    distance_max: float,
) -> int:
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

    fields = []
    fields.append(
        add_threshold_field(
            1,
            groups["wall_chamber"],
            preset["wall"],
            preset["wall_far"],
            preset["wall_transition"],
        )
    )
    fields.append(
        add_threshold_field(
            3,
            groups["inlet_torch"] + groups["inlet_spray"] + groups["wall_plate"],
            preset["inlet"],
            preset["wall_far"],
            preset["wall_transition"],
        )
    )

    # Throat and near-plume boxes retain the requested local resolution while
    # allowing the one-metre farfield to coarsen without a point explosion.
    gmsh.model.mesh.field.add("Box", 5)
    for key, value in {
        "VIn": preset["throat"],
        "VOut": preset["global_max"],
        "XMin": 258.0,
        "XMax": 292.0,
        "YMin": -32.0,
        "YMax": 32.0,
        "ZMin": -32.0,
        "ZMax": 32.0,
        "Thickness": 8.0,
    }.items():
        gmsh.model.mesh.field.setNumber(5, key, value)
    fields.append(5)

    gmsh.model.mesh.field.add("Box", 6)
    for key, value in {
        "VIn": preset["shock"],
        "VOut": preset["farfield"],
        "XMin": 360.0,
        "XMax": 520.0,
        "YMin": -70.0,
        "YMax": 70.0,
        "ZMin": -70.0,
        "ZMax": 70.0,
        "Thickness": 20.0,
    }.items():
        gmsh.model.mesh.field.setNumber(6, key, value)
    fields.append(6)

    gmsh.model.mesh.field.add("Min", 7)
    gmsh.model.mesh.field.setNumbers(7, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(7)


def mesh_quality() -> dict[str, float | int]:
    element_types, element_tags_by_type, _ = gmsh.model.mesh.getElements(3)
    tetra_tags = []
    for element_type, tags in zip(element_types, element_tags_by_type):
        if gmsh.model.mesh.getElementProperties(element_type)[0].startswith("Tetrahedron"):
            tetra_tags.extend(tags)
    if not tetra_tags:
        raise RuntimeError("volume mesh contains no tetrahedra")
    qualities = gmsh.model.mesh.getElementQualities(tetra_tags, "minSICN")
    ordered = sorted(float(value) for value in qualities)
    p01 = ordered[max(0, int(0.01 * (len(ordered) - 1)))]
    result = {
        "tetrahedra": len(tetra_tags),
        "min_sicn": ordered[0],
        "p01_sicn": p01,
        "mean_sicn": sum(ordered) / len(ordered),
    }
    if result["min_sicn"] <= 0.0:
        raise RuntimeError(f"inverted tetrahedron detected: {result}")
    return result


def main() -> None:
    args = parse_args()
    profile_path = resolve_profile(args.profile)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", int(args.terminal))
        gmsh.model.add("vortex_full360_v2")
        volume = build_geometry(profile_path)
        groups = classify_boundaries(volume)
        add_physical_groups(volume, groups)

        volume_mm3 = gmsh.model.occ.getMass(3, volume)
        topology = {
            "profile": str(profile_path),
            "volume_mm3": volume_mm3,
            "boundary_surfaces": {name: tags for name, tags in groups.items()},
        }
        print("TOPOLOGY", json.dumps(topology, sort_keys=True))

        gmsh.write(str(output_dir / "vortex_full360_v2.step"))
        gmsh.write(str(output_dir / "vortex_full360_v2.brep"))
        if args.mode == "geometry":
            print("VORTEX_GEOMETRY_OK")
            return

        configure_mesh(groups, args.preset)
        gmsh.model.mesh.generate(2)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.write(str(output_dir / "vortex_full360_v2_surface.msh"))
        if args.mode == "surface":
            print("VORTEX_SURFACE_MESH_OK")
            return

        gmsh.model.mesh.generate(3)
        quality = mesh_quality()
        print("QUALITY", json.dumps(quality, sort_keys=True))
        gmsh.write(str(output_dir / "vortex_full360_v2_volume.msh"))
        with (output_dir / "vortex_full360_v2_report.json").open("w", encoding="utf-8") as stream:
            json.dump({**topology, "preset": args.preset, "quality": quality}, stream, indent=2)
            stream.write("\n")
        print("VORTEX_VOLUME_MESH_OK")
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    main()
