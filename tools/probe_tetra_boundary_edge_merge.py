#!/usr/bin/env python3
"""Probe one boundary-node merge on a tetra NPZ without changing the source mesh.

The operation is intentionally explicit: it merges one supplied node pair,
removes only degenerate tetrahedra/boundary triangles caused by that merge,
and emits a temporary NPZ for the normal RFAM mesh/export audits.
"""
import argparse
import json
from pathlib import Path

import numpy as np


LOCAL_FACES = np.array([[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]])


def verify_boundary_faces(cells, boundary, owners):
    # Match the entire exterior face table, including its owner, before export.
    keys = np.sort(cells[:, LOCAL_FACES].reshape(-1, 3), axis=1)
    unique, first, counts = np.unique(keys, axis=0, return_index=True,
                                      return_counts=True)
    if np.any(counts > 2):
        raise ValueError('nonmanifold face after merge')
    exterior = counts == 1
    boundary_keys = np.sort(boundary, axis=1)
    order = np.lexsort(boundary_keys.T[::-1])
    if (not np.array_equal(unique[exterior], boundary_keys[order])
            or not np.array_equal(first[exterior] // 4, owners[order])):
        raise ValueError('candidate boundary faces or owners differ from reconstructed exterior')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--keep', type=int, required=True)
    ap.add_argument('--remove', type=int, required=True)
    ap.add_argument('--report', required=True)
    args = ap.parse_args()
    data = np.load(args.input)
    points = np.array(data['points'], dtype=np.float64, copy=True)
    cells = np.array(data['cells'], dtype=np.int64, copy=True)
    boundary = np.array(data['boundary'], dtype=np.int64, copy=True)
    boundary_zone = np.array(data['boundary_zone'], dtype=np.int32, copy=True)
    boundary_owner = np.array(data['boundary_owner'], dtype=np.int64, copy=True)
    inlet_patch = np.array(data['inlet_patch'], dtype=np.int32, copy=True)
    def area_by_zone(pts, faces, zone_ids):
        tri = pts[faces]
        area = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0],
                                       tri[:, 2] - tri[:, 0]), axis=1) / 2.
        return {int(zone): float(area[zone_ids == zone].sum())
                for zone in np.unique(zone_ids)}

    original_areas = area_by_zone(points, boundary, boundary_zone)
    keep, remove = args.keep, args.remove
    if keep == remove or not (0 <= keep < len(points)) or not (0 <= remove < len(points)):
        raise ValueError('invalid node pair')
    if not (np.any(boundary == keep) and np.any(boundary == remove)):
        raise ValueError('both nodes must be on the boundary')
    touching = np.any(np.isin(boundary, [keep, remove]), axis=1)
    if np.any(boundary_zone[touching] != 0):
        raise ValueError('merge may only affect wall zone 0')
    separation = float(np.linalg.norm(points[keep] - points[remove]))
    if not 0 < separation < .001:
        raise ValueError('boundary edge must be shorter than 1 micron')
    points[keep] = .5 * (points[keep] + points[remove])
    cells[cells == remove] = keep
    cell_key = np.sort(cells, axis=1)
    degenerate_cells = np.any(np.diff(cell_key, axis=1) == 0, axis=1)
    if np.any(degenerate_cells):
        old_to_new = np.full(len(cells), -1, dtype=np.int64)
        old_to_new[~degenerate_cells] = np.arange(np.count_nonzero(~degenerate_cells))
        cells = cells[~degenerate_cells]
    else:
        old_to_new = np.arange(len(cells), dtype=np.int64)
    boundary = boundary.copy()
    boundary[boundary == remove] = keep
    boundary_key = np.sort(boundary, axis=1)
    degenerate_faces = np.any(np.diff(boundary_key, axis=1) == 0, axis=1)
    duplicate_faces = np.zeros(len(boundary), dtype=bool)
    if np.count_nonzero(~degenerate_faces):
        keys = boundary_key[~degenerate_faces]
        _, first = np.unique(keys, axis=0, return_index=True)
        duplicate_faces[~degenerate_faces] = True
        valid_indices = np.flatnonzero(~degenerate_faces)[np.sort(first)]
        duplicate_faces[:] = True
        duplicate_faces[valid_indices] = False
    face_keep = ~(degenerate_faces | duplicate_faces)
    boundary = boundary[face_keep]
    boundary_zone = boundary_zone[face_keep]
    inlet_patch = inlet_patch[face_keep]
    boundary_owner = old_to_new[boundary_owner[face_keep]]
    if np.any(boundary_owner < 0):
        raise ValueError('boundary owner was removed by the merge')
    used = np.unique(cells)
    remap = np.full(len(points), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    points = points[used]
    cells = remap[cells]
    boundary = remap[boundary]
    if np.any(cells < 0) or np.any(boundary < 0):
        raise ValueError('dangling node after remap')
    cell_keys = np.sort(cells, axis=1)
    if len(np.unique(cell_keys, axis=0)) != len(cells):
        raise ValueError('merge created duplicate tetrahedra')
    face_keys = np.sort(boundary, axis=1)
    if len(np.unique(face_keys, axis=0)) != len(boundary):
        raise ValueError('merge left duplicate boundary faces')
    verify_boundary_faces(cells, boundary, boundary_owner)
    new_areas = area_by_zone(points, boundary, boundary_zone)
    if original_areas.keys() != new_areas.keys():
        raise ValueError('boundary zones changed')
    for zone, old_area in original_areas.items():
        delta = abs(new_areas[zone] - old_area)
        if (zone != 0 and delta > 1e-10 * max(old_area, 1.)) or (zone == 0 and delta > 1e-3):
            raise ValueError(f'boundary zone {zone} area changed by {delta} mm2')
    x = points[cells]
    signed = np.einsum('ij,ij->i', x[:, 1] - x[:, 0],
                       np.cross(x[:, 2] - x[:, 0], x[:, 3] - x[:, 0])) / 6.
    if np.any(signed <= 0):
        raise ValueError('merge produced nonpositive tetra volume')
    # The exporter reconstructs shared faces from cells and maps boundary faces
    # by sorted node key, so this output keeps the normal audited NPZ contract.
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, points=points, cells=cells, boundary=boundary,
                        boundary_zone=boundary_zone, boundary_owner=boundary_owner,
                        inlet_patch=inlet_patch)
    report = {
        'input': args.input,
        'output': str(out),
        'keep_node': keep,
        'remove_node': remove,
        'original_separation_mm': separation,
        'removed_cells': int(np.count_nonzero(degenerate_cells)),
        'removed_boundary_faces': int(np.count_nonzero(degenerate_faces)),
        'removed_duplicate_boundary_faces': int(np.count_nonzero(duplicate_faces)),
        'output_points': int(len(points)),
        'output_cells': int(len(cells)),
        'output_boundary_faces': int(len(boundary)),
        'positive_volumes': True,
        'exterior_faces_and_owners_match': True,
        'area_change_mm2_by_zone': {str(zone): new_areas[zone] - old_area
                                    for zone, old_area in original_areas.items()},
        'scope': 'wall-node merge trial only; not an authorized production mesh',
    }
    Path(args.report).write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('TETRA_BOUNDARY_EDGE_MERGE_PROBE', json.dumps(report, sort_keys=True))


if __name__ == '__main__':
    main()
