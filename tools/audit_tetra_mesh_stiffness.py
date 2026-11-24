#!/usr/bin/env python3
"""Audit the delivered A tetra mesh's small-cell stiffness and inlet hotspots.

The RFAMTET1 mesh is checked against the paired NPZ before reporting metrics.
This is a read-only geometry audit; acoustic dt assumes uniform stationary gas.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

EDGES = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
INLET_ZONES = [10, 11, 12, 20, 21, 22]


def percentiles(values):
    return dict(zip(['min', 'p001', 'p01', 'p1', 'p50', 'max'],
                    [float(v) for v in np.percentile(values, [0, .001, .01, 1, 50, 100])]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--npz', required=True)
    parser.add_argument('--rfamtet', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--sound-speed', type=float, default=352.922)
    parser.add_argument('--cfl', type=float, default=.5)
    parser.add_argument('--hotspot-cells', type=int, nargs='*', default=[])
    args = parser.parse_args()
    if args.sound_speed <= 0 or args.cfl <= 0:
        parser.error('sound speed and CFL must be positive')
    with np.load(args.npz) as data:
        points = data['points']
        cells = data['cells']
        boundary = data['boundary']
        owners = data['boundary_owner']
        zones = data['boundary_zone']
    npoints, ncells = len(points), len(cells)
    raw = np.memmap(args.rfamtet, dtype=np.uint8, mode='r')
    if bytes(raw[:8]) != b'RFAMTET1':
        raise ValueError('not a RFAMTET1 mesh')
    counts = np.frombuffer(raw, dtype='<u8', count=3, offset=8)
    if int(counts[0]) != npoints or int(counts[1]) != ncells:
        raise ValueError('NPZ and RFAMTET1 mesh counts differ')
    offset = 32
    binary_points = np.ndarray((npoints, 3), '<f8', raw, offset)
    offset += npoints * 24
    binary_cells = np.ndarray((ncells, 4), '<i4', raw, offset)
    offset += ncells * 16
    binary_volumes = np.ndarray(ncells, '<f8', raw, offset)
    offset += ncells * 8
    if raw.size != offset + int(counts[2]) * 24:
        raise ValueError('RFAMTET1 face table size mismatch')
    if not np.array_equal(binary_cells, cells):
        raise ValueError('NPZ and RFAMTET1 cell ordering differ')
    if not np.allclose(binary_points * 1000., points, rtol=0., atol=1e-9):
        raise ValueError('NPZ and RFAMTET1 coordinates differ')

    volumes = np.empty(ncells)
    quality = np.empty(ncells)
    min_heights = np.empty(ncells)
    geometric_dt = np.empty(ncells)
    max_edge_ratio = np.empty(ncells)
    for start in range(0, ncells, 100000):
        end = min(start + 100000, ncells)
        x = points[cells[start:end]]
        v = np.abs(np.einsum('ij,ij->i', x[:, 1] - x[:, 0],
                             np.cross(x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]))) / 6.
        edges = np.linalg.norm(x[:, EDGES[:, 0]] - x[:, EDGES[:, 1]], axis=2)
        area = np.stack([np.linalg.norm(np.cross(x[:, j] - x[:, i],
                                               x[:, k] - x[:, i]), axis=1) / 2.
                         for i, j, k in FACES], axis=1)
        if np.any(v <= 0) or np.any(area <= 0):
            raise ValueError('nonpositive cell volume or face area')
        volumes[start:end] = v
        quality[start:end] = 12. * np.power(3. * v, 2. / 3.) / np.sum(edges**2, axis=1)
        min_heights[start:end] = np.min(3. * v[:, None] / area, axis=1)
        geometric_dt[start:end] = args.cfl * v * 1e-9 / (area.sum(axis=1) * 1e-6 * args.sound_speed)
        max_edge_ratio[start:end] = np.max(edges, axis=1) / np.min(edges, axis=1)
    if not np.allclose(binary_volumes, volumes * 1e-9, rtol=1e-6, atol=1e-24):
        raise ValueError('NPZ and RFAMTET1 cell volumes differ')
    binary_faces = np.ndarray((int(counts[2]), 6), '<i4', raw, offset)
    exterior = binary_faces[binary_faces[:, 4] == -1]
    if len(exterior) != len(boundary):
        raise ValueError('NPZ and RFAMTET1 boundary face counts differ')
    binary_key = np.sort(exterior[:, :3], axis=1)
    npz_key = np.sort(boundary, axis=1)
    binary_order = np.lexsort(binary_key.T[::-1])
    npz_order = np.lexsort(npz_key.T[::-1])
    if (not np.array_equal(binary_key[binary_order], npz_key[npz_order])
            or not np.array_equal(exterior[binary_order, 3], owners[npz_order])
            or not np.array_equal(exterior[binary_order, 5], zones[npz_order])):
        raise ValueError('NPZ and RFAMTET1 boundary faces, owners or zones differ')
    inlet_mask = np.isin(zones, INLET_ZONES)
    inlet_owners = np.unique(owners[inlet_mask])

    def cell_row(index):
        index = int(index)
        if not 0 <= index < ncells:
            raise ValueError(f'invalid hotspot cell {index}')
        own_zones = np.unique(zones[owners == index]).astype(int).tolist()
        return {'cell': index,
                'centroid_mm': points[cells[index]].mean(axis=0).tolist(),
                'volume_mm3': float(volumes[index]),
                'mean_ratio': float(quality[index]),
                'min_height_mm': float(min_heights[index]),
                'geometric_cfl_s': float(geometric_dt[index]),
                'boundary_zones': own_zones,
                'inlet_owner': bool(np.isin(index, inlet_owners))}

    with open(args.rfamtet, 'rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    report = {
        'rfamtet_sha256': digest,
        'mesh_identity_verified': True,
        'boundary_owner_and_zone_verified': True,
        'cells': ncells,
        'cfl': args.cfl,
        'sound_speed_m_s': args.sound_speed,
        'geometric_cfl_assumption': 'uniform stationary gas, constant sound speed; excludes inlet trace and viscosity',
        'quality_mean_ratio': percentiles(quality),
        'min_height_mm': percentiles(min_heights),
        'geometric_cfl_s': percentiles(geometric_dt),
        'max_edge_ratio': percentiles(max_edge_ratio),
        'quality_below': {str(v): int(np.count_nonzero(quality < v)) for v in [.01, .02, .05, .1]},
        'height_below_mm': {str(v): int(np.count_nonzero(min_heights < v)) for v in [.0005, .001, .005, .01]},
        'cfl_below_s': {str(v): int(np.count_nonzero(geometric_dt < v)) for v in [1e-10, 1e-9, 1e-8]},
        'inlet_owner_count': len(inlet_owners),
        'min_inlet_owner_geometric_cfl_s': float(geometric_dt[inlet_owners].min()),
        'lowest_geometric_cfl_cells': [cell_row(i) for i in np.argsort(geometric_dt)[:10]],
        'hotspots': [cell_row(i) for i in args.hotspot_cells],
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('TETRA_MESH_STIFFNESS_AUDIT', digest,
          'geometric_cfl_min_s', geometric_dt.min(),
          'height_below_1um', report['height_below_mm']['0.001'])
    for row in report['hotspots']:
        print('HOTSPOT', row)


if __name__ == '__main__':
    main()
