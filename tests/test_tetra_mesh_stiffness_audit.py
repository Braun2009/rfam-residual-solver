#!/usr/bin/env python3
"""The RFAMTET1/NPZ audit must reject altered boundary ownership and zones."""
import gc
import struct
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from audit_tetra_mesh_stiffness import main as audit_main


def run_audit(npz, binary, report):
    saved = sys.argv
    sys.argv = ['audit_tetra_mesh_stiffness.py', '--npz', str(npz),
                '--rfamtet', str(binary), '--out', str(report)]
    try:
        audit_main()
    finally:
        sys.argv = saved
        gc.collect()


def main():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.],
                           [0., 0., 1.], [0., 0., -1.]])
        cells = np.array([[0, 1, 2, 3], [0, 2, 1, 4]], dtype='<i4')
        faces = np.array([[0, 2, 1, 0, 1, 0],
                          [0, 1, 3, 0, -1, 10], [1, 2, 3, 0, -1, 0],
                          [2, 0, 3, 0, -1, 0], [0, 2, 4, 1, -1, 0],
                          [2, 1, 4, 1, -1, 0], [1, 0, 4, 1, -1, 0]], dtype='<i4')
        boundary = faces[1:, :3]
        zones = faces[1:, 5].copy()
        owners = faces[1:, 3].copy()
        npz = root / 'mesh.npz'
        binary = root / 'mesh.rfamtet'
        np.savez(npz, points=points * 1000., cells=cells,
                 boundary=boundary, boundary_owner=owners, boundary_zone=zones)
        with binary.open('wb') as out:
            out.write(b'RFAMTET1')
            out.write(struct.pack('<QQQ', len(points), len(cells), len(faces)))
            out.write(points.astype('<f8').tobytes())
            out.write(cells.tobytes())
            out.write(np.full(len(cells), 1. / 6., dtype='<f8').tobytes())
            out.write(faces.tobytes())
        run_audit(npz, binary, root / 'ok.json')
        changed = zones.copy()
        changed[0] = 11
        np.savez(npz, points=points * 1000., cells=cells,
                 boundary=boundary, boundary_owner=owners, boundary_zone=changed)
        try:
            run_audit(npz, binary, root / 'bad.json')
        except ValueError as error:
            if 'boundary faces, owners or zones differ' not in str(error):
                raise
        else:
            raise AssertionError('tampered boundary zone was admitted')
        print('TETRA_MESH_STIFFNESS_AUDIT_PASS')


if __name__ == '__main__':
    main()
