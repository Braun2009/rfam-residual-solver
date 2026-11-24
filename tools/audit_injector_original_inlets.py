#!/usr/bin/env python3
"""Audit BCAD-selected solid faces against original injector inlet apertures.

The solid STL has no fluid inlet caps at its open holes. Selected lip/wall
triangles anchor semantics; only a reconstructed fluid domain can expose BCs.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
import trimesh


AXIS_XY = np.array([49.21534822760748, 48.35004867381059])
BODY_SHA256 = '89a9579104f1ef532c7814fa6c612205e5a69c4ce1fa037163e5c6ba7a67f7d6'
SELECTOR_SHA256 = 'd3186d9b9df60b15ec7ddf029b7414f50764c9ec2af0c982f39de872a4b23ed1'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as src:
        for block in iter(lambda: src.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def boundary_loops(triangles):
    """Return vertex loops of a selected triangular solid-surface patch."""
    edges = defaultdict(int)
    for tri in triangles:
        verts = [tuple(np.round(v, 7)) for v in tri]
        for a, b in ((0, 1), (1, 2), (2, 0)):
            edges[tuple(sorted((verts[a], verts[b])))] += 1
    if any(count > 2 for count in edges.values()):
        raise ValueError('nonmanifold selected-face edge')
    neighbours = defaultdict(set)
    for (a, b), count in edges.items():
        if count == 1:
            neighbours[a].add(b)
            neighbours[b].add(a)
    if not neighbours or any(len(other) != 2 for other in neighbours.values()):
        raise ValueError('selected face patch has open/branched rim')
    components = []
    visited = set()
    for first in neighbours:
        if first in visited:
            continue
        stack = [first]
        visited.add(first)
        component = []
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for other in neighbours[vertex]:
                if other not in visited:
                    visited.add(other)
                    stack.append(other)
        components.append(np.asarray(component))
    return sorted(components, key=len, reverse=True)


def section_holes(mesh, z_mm):
    section = mesh.section(plane_origin=[0, 0, z_mm], plane_normal=[0, 0, 1])
    if section is None:
        raise ValueError(f'no STL cross-section at z={z_mm}')
    holes = []
    for path in section.discrete:
        if len(path) < 4 or np.linalg.norm(path[0] - path[-1]) > 1e-5:
            raise ValueError(f'cross-section has an open contour at z={z_mm}')
        area = abs(float(np.sum(path[:-1, 0] * path[1:, 1]
                                - path[1:, 0] * path[:-1, 1]) / 2))
        xy = np.mean(path[:-1, :2], axis=0)
        radius = float(np.linalg.norm(xy - AXIS_XY))
        if 4 < area < 7 and 21 < radius < 26:
            angle = float(np.degrees(np.arctan2(xy[1] - AXIS_XY[1],
                                                xy[0] - AXIS_XY[0])))
            holes.append((area, radius, angle))
    if len(holes) != 48:
        raise ValueError(f'expected 48 propane hole contours at z={z_mm}, got {len(holes)}')
    holes.sort(key=lambda item: item[2])
    angles = np.asarray([item[2] for item in holes])
    gaps = np.diff(np.r_[angles, angles[0] + 360])
    if np.max(np.abs(gaps - 7.5)) > 0.3:
        raise ValueError(f'propane hole pitch is not 7.5 degrees at z={z_mm}')
    return {'z_mm': z_mm, 'count': 48,
            'contour_area_mm2_min_median_max': np.percentile([x[0] for x in holes], [0, 50, 100]).tolist(),
            'contour_centre_radius_mm_min_median_max': np.percentile([x[1] for x in holes], [0, 50, 100]).tolist(),
            'angular_pitch_deg_min_median_max': np.percentile(gaps, [0, 50, 100]).tolist(),
            'note': 'Internal STL cross-section; neither fluid supply boundary area nor injector exit area.'}


def audit(body_path, selector_path, *, expected_body=BODY_SHA256, expected_selector=SELECTOR_SHA256):
    if digest(body_path) != expected_body or digest(selector_path) != expected_selector:
        raise ValueError('source body or BCAD selector hash changed')
    selector = json.loads(Path(selector_path).read_text(encoding='utf-8'))
    if (selector.get('schema_version') != 'bcad.stable-selector/v0'
            or selector.get('mesh_id') != 'injector_v3' or selector.get('kind') != 'faces'):
        raise ValueError('unexpected BCAD selector identity')
    raw = trimesh.load_mesh(str(body_path), process=False)
    indices = np.asarray(selector['faces'], dtype=np.int64)
    if len(indices) != 1478 or len(np.unique(indices)) != 1478 or indices.min() < 0 or indices.max() >= len(raw.faces):
        raise ValueError('selected face count/indices changed')
    selected = raw.triangles[indices]
    coordinates = np.asarray(selector['face_coordinates'], dtype=float)
    if coordinates.shape != selected.shape or not np.allclose(coordinates, selected, atol=1e-5, rtol=0):
        raise ValueError('BCAD selector coordinates do not match STL face order')
    z = selected.mean(axis=1)[:, 2]
    cap, wall = np.abs(z) < 1e-6, (z > 250) & (z < 270)
    if cap.sum() != 76 or wall.sum() != 1402 or not np.all(cap | wall):
        raise ValueError('LOX ring and propane-hole selection changed')
    cap_triangles = selected[cap]
    if not np.all(np.abs(cap_triangles[:, :, 2]) < 1e-6):
        raise ValueError('LOX lip is not planar at z=0')
    cap_loops = boundary_loops(cap_triangles)
    if len(cap_loops) != 2:
        raise ValueError('LOX ring must have two boundary loops')
    sizes = sorted([len(loop) for loop in cap_loops])
    if sizes != [38, 38]:
        raise ValueError('LOX lip loop discretization changed')
    radii = sorted((np.linalg.norm(loop[:, :2] - AXIS_XY, axis=1) for loop in cap_loops),
                   key=np.mean)
    inner, outer = radii
    if not (np.max(np.abs(inner - 4.55)) < 2e-4
            and np.max(np.abs(outer - 6.35)) < 2e-4):
        raise ValueError('LOX ring radii changed')
    cap_area = np.linalg.norm(np.cross(cap_triangles[:, 1] - cap_triangles[:, 0],
                                       cap_triangles[:, 2] - cap_triangles[:, 0]), axis=1).sum() / 2
    if not np.isclose(cap_area, 61.348230199721, atol=1e-6):
        raise ValueError('LOX solid ring area changed')
    propane_loops = boundary_loops(selected[wall])
    if sorted(map(len, propane_loops)) != [36, 80]:
        raise ValueError('propane lip/wall patch rims changed')
    verts, inverse = np.unique(raw.vertices, axis=0, return_inverse=True)
    solid = trimesh.Trimesh(vertices=verts, faces=inverse[raw.faces], process=False)
    if not solid.is_watertight or not solid.is_winding_consistent:
        raise ValueError('original injector solid is not closed/oriented')
    sections = [section_holes(solid, z_mm) for z_mm in (261.0, 261.5)]
    return {
        'schema': 'rfam.injector-original-inlet-geometry/v1',
        'status': 'SELECTED_SOLID_LIPS_CONFIRMED_FLUID_BOUNDARIES_NOT_BUILT',
        'body': {'path': str(body_path), 'sha256': expected_body},
        'selector': {'path': str(selector_path), 'sha256': expected_selector,
                     'faces': len(indices)},
        'lox': {'selected_solid_faces': int(cap.sum()), 'solid_annular_lip_area_mm2': float(cap_area),
                'lip_z_mm': 0, 'inner_radius_mm_min_max': [float(inner.min()), float(inner.max())],
                'outer_radius_mm_min_max': [float(outer.min()), float(outer.max())],
                'inner_void_cross_section_area_mm2': float(np.pi * 4.55 ** 2),
                'inner_void_polygon_area_mm2': float(38 * .5 * 4.55 ** 2 * np.sin(2 * np.pi / 38)),
                'semantic': 'user-selected metal annulus around central hole; CFD inlet is a NEW cap inside inner rim, not these STL faces',
                'fluid_cap_built': False, 'passage_connectivity_validated': False},
        'propane': {'selected_solid_faces': int(wall.sum()), 'selected_feature': 'one rounded opening and its hole wall',
                    'selected_boundary_loop_vertices': sorted(map(len, propane_loops)),
                    'sections': sections, 'fluid_caps_built': False,
                    'passage_connectivity_validated': False},
        'delivered_mesh_inlet_z_mm': 269.2,
        'delivered_mesh_contains_original_supply_openings': False,
        'torch_two_feed_geometry_resolved': False,
        'production_validated': False, 'ignition_authorized': False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--body', required=True)
    parser.add_argument('--selector', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    result = audit(args.body, args.selector)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print('ORIGINAL_INLET_GEOMETRY_PASS', out)


if __name__ == '__main__':
    main()
