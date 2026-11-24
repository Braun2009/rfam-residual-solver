#!/usr/bin/env python3
"""Gate whether a tetra delivery has semantic CAD injector exit surfaces.

This tool deliberately does not invent surfaces from z cuts or connectivity
mergers. It reports the currently encoded boundary zones and declares the
exit-control-surface gate open only when explicit named exit labels exist.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh', required=True, help='fluid_volume_mm.npz')
    ap.add_argument('--geometry-manifest', required=True)
    ap.add_argument('--connectivity', required=False)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    mesh_path = Path(args.mesh)
    with np.load(mesh_path) as d:
        zones = np.asarray(d['boundary_zone'], dtype=np.int32)
        points = np.asarray(d['points'], dtype=float)
        boundary = np.asarray(d['boundary'], dtype=np.int64)
    manifest = json.loads(Path(args.geometry_manifest).read_text(encoding='utf-8'))
    zone_names = {0: 'wall', 1: 'external_outlet', 10: 'lox_0_inlet',
                  11: 'lox_1_inlet', 12: 'lox_2_inlet', 20: 'fuel_0_inlet',
                  21: 'fuel_1_inlet', 22: 'fuel_2_inlet', 30: 'torch_inlet'}
    counts = {str(int(z)): int(np.count_nonzero(zones == z))
              for z in np.unique(zones)}
    unknown = sorted(set(int(z) for z in np.unique(zones)) - set(zone_names))
    connectivity = None
    if args.connectivity:
        connectivity = json.loads(Path(args.connectivity).read_text(encoding='utf-8'))
    report = {
        'mesh_sha256': sha256(mesh_path),
        'mesh_points': int(len(points)),
        'mesh_boundary_triangles': int(len(boundary)),
        'encoded_boundary_zones': counts,
        'encoded_boundary_zone_names': {str(k): zone_names.get(k, 'unknown') for k in counts},
        'unknown_boundary_zones': unknown,
        'true_exit_surfaces_available': False,
        'exit_surface_labels': [],
        'scope': 'semantic geometry gate; no axial cut or connectivity merger is promoted to an exit',
        'tip_z_mm': [float(p.get('tip_cad_z_mm')) for p in manifest.get('pintles', [])],
        'connectivity_evidence': connectivity,
        'required_exit_labels': ['lox_slot_exit_0', 'lox_slot_exit_1',
                                 'lox_slot_exit_2', 'fuel_hole_exit_0..59',
                                 'torch_exit'],
        'action': 'obtain named CAD/Fluent exit faces or an independently validated curved-surface semantic map',
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('TETRA_EXIT_SEMANTICS_MISSING',
          'zones=' + ','.join(sorted(counts)),
          'mesh_sha256=' + report['mesh_sha256'])


if __name__ == '__main__':
    main()
