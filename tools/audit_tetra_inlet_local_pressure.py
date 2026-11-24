#!/usr/bin/env python3
"""Diagnostic owner-pressure sensitivity of inlet-face supply capacity.

The solver imposes a single area-weighted pressure per zone, not each owner's
pressure at its boundary face. Local capacity is a model sensitivity, NOT an
admitted face boundary condition or proof that the current zone BC failed.
The ideal capacity is an upper bound, not a measured discharge coefficient.
"""
import argparse
import json
import struct
import sys
from pathlib import Path

import numpy as np
from CoolProp.CoolProp import PropsSI

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tetra_field_to_vtu import read_field, read_mesh

SUPPLIES = {
    10: ('Oxygen', 6.0e6, 133.5, 2.8 / 3),
    11: ('Oxygen', 6.0e6, 130., 2.8 / 3),
    12: ('Oxygen', 6.0e6, 130., 2.8 / 3),
    20: ('Propane', 4.5e6, 453.5, 1.05 / 3),
    21: ('Propane', 4.5e6, 453.5, 1.05 / 3),
    22: ('Propane', 4.5e6, 453.5, 1.05 / 3),
}


def audit(mesh_path, field_path, supplies=SUPPLIES):
    points, cells, _ = read_mesh(mesh_path)
    centroids, fields, _, _ = read_field(field_path)
    if len(cells) != len(fields):
        raise ValueError('field/mesh cell counts differ')
    raw = Path(mesh_path).read_bytes()
    npoints, ncells, nfaces = struct.unpack_from('<QQQ', raw, 8)
    if npoints != len(points) or ncells != len(cells):
        raise ValueError('mesh header inconsistent')
    faces = np.frombuffer(raw, dtype='<i4', count=nfaces * 6,
                          offset=32 + npoints * 24 + ncells * 24).reshape(-1, 6)
    result = {'mesh': str(mesh_path), 'field': str(field_path),
              'converged': False,
              'scope': 'diagnostic owner-pressure capacity, not the zone-uniform solver inlet',
              'zones': {}}
    for zone, (fluid, supply_p, supply_t, design_flow) in supplies.items():
        selected = faces[(faces[:, 4] == -1) & (faces[:, 5] == zone)]
        if not len(selected):
            raise ValueError(f'missing inlet zone {zone}')
        x = points[selected[:, :3]]
        areas = np.linalg.norm(np.cross(x[:, 1] - x[:, 0], x[:, 2] - x[:, 0]),
                               axis=1) * .5
        p = fields[selected[:, 3], 0].astype(np.float64)
        blocked = p >= supply_p
        h0 = PropsSI('H', 'P', supply_p, 'T', supply_t, fluid)
        s0 = PropsSI('S', 'P', supply_p, 'T', supply_t, fluid)
        capacity = 0.
        face_capacities = np.zeros(len(p), dtype=float)
        for index, (pressure, area) in enumerate(zip(p, areas)):
            if blocked[index]:
                continue
            density = PropsSI('D', 'P', float(pressure), 'S', s0, fluid)
            h = PropsSI('H', 'P', float(pressure), 'S', s0, fluid)
            face_capacities[index] = area * density * np.sqrt(max(0., 2 * (h0 - h)))
            capacity += face_capacities[index]
        fixed_face_flow = design_flow * areas / areas.sum()
        insufficient = face_capacities < fixed_face_flow
        result['zones'][str(zone)] = {
            'scope': 'owner-pressure alternative; not the imposed zone-uniform face pressure',
            'supply_T_K': supply_t,
            'faces': len(selected), 'area_m2': float(areas.sum()),
            'area_weighted_p_Pa': float(np.average(p, weights=areas)),
            'p_min_Pa': float(p.min()), 'p_p01_Pa': float(np.percentile(p, 1)),
            'p_p50_Pa': float(np.percentile(p, 50)),
            'p_p99_Pa': float(np.percentile(p, 99)), 'max_local_p_Pa': float(p.max()),
            'capacity_min_face_kg_s': float(face_capacities[~blocked].min()) if np.any(~blocked) else 0.,
            'capacity_p01_face_kg_s': float(np.percentile(face_capacities[~blocked], 1)) if np.any(~blocked) else 0.,
            'capacity_p50_face_kg_s': float(np.percentile(face_capacities[~blocked], 50)) if np.any(~blocked) else 0.,
            'capacity_p99_face_kg_s': float(np.percentile(face_capacities[~blocked], 99)) if np.any(~blocked) else 0.,
            'blocked_area_fraction': float(areas[blocked].sum() / areas.sum()),
            'fixed_face_local_capacity_shortfall_count': int(np.count_nonzero(insufficient)),
            'fixed_face_local_capacity_shortfall_area_fraction': float(areas[insufficient].sum() / areas.sum()),
            'fixed_face_local_capacity_shortfall_kg_s': float(np.maximum(fixed_face_flow - face_capacities, 0).sum()),
            'ideal_available_capacity_kg_s': float(capacity),
            'design_target_kg_s': design_flow,
            'target_below_ideal_capacity': bool(design_flow <= capacity),
        }
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh', required=True)
    ap.add_argument('--field', required=True)
    ap.add_argument('--out')
    ap.add_argument('--lox-p', type=float, default=6e6)
    ap.add_argument('--lox-T', type=float, default=133.5)
    ap.add_argument('--lox-mdot', type=float, default=2.8)
    ap.add_argument('--fuel-p', type=float, default=4.5e6)
    ap.add_argument('--fuel-T', type=float, default=453.5)
    ap.add_argument('--fuel-mdot', type=float, default=1.05)
    args = ap.parse_args()
    inputs = [args.lox_p, args.lox_T, args.lox_mdot,
              args.fuel_p, args.fuel_T, args.fuel_mdot]
    if not all(np.isfinite(v) and v > 0 for v in inputs):
        ap.error('all supply values must be finite and positive')
    supplies = {z: ('Oxygen', args.lox_p, args.lox_T, args.lox_mdot/3)
                for z in (10, 11, 12)}
    supplies.update({z: ('Propane', args.fuel_p, args.fuel_T, args.fuel_mdot/3)
                     for z in (20, 21, 22)})
    result = audit(args.mesh, args.field, supplies)
    text = json.dumps(result, indent=2)
    if args.out:
        Path(args.out).write_text(text + '\n', encoding='utf-8')
    print(text)


if __name__ == '__main__':
    main()
