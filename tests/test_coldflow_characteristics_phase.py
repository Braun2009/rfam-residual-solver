#!/usr/bin/env python3
"""Check EOS phase labels against a small RFAMTETF1 field fixture."""
from __future__ import annotations
import json
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        mesh, field = root / 'mesh.rfamtet', root / 'field.rfamtetf'
        points = np.array([[0, 0, .27090], [1e-5, 0, .27090],
                           [0, 1e-5, .27090], [0, 0, .27100]], dtype='<f8')
        cells = np.array([[0, 1, 2, 3], [0, 1, 2, 3], [0, 1, 2, 3]], dtype='<i4')
        vol = np.array([1e-15, 2e-15, 3e-15], dtype='<f8')
        mesh.write_bytes(b'RFAMTET1' + struct.pack('<QQQ', 4, 3, 0)
                         + points.tobytes() + cells.tobytes() + vol.tobytes())
        centroids = np.array([[0, 0, .27090], [0, 0, .27100], [0, 0, .27101]], dtype='<f4')
        values = np.array([[4e6, 130, 0, 0, 2, 930, 500, -100000],
                           [4e6, 130, 0, 0, 3, 930, 500, -100000],
                           [4e6, 400, 0, 0, 4, 40, 350, 10000]], dtype='<f4')
        y = np.array([[1, 0], [1, 0], [.784, .216]], dtype='<f4')
        field.write_bytes(b'RFAMTETF1' + struct.pack('<QQ', 3, 2)
                          + centroids.tobytes() + values.tobytes() + y.tobytes())
        csv = Path(str(field) + '.phase.csv')
        csv.write_text('cell,region,vapor_mass_fraction,vapor_volume_fraction,acoustic_evidence,projected\n'
                       '0,5,0.3,0.6,2,0\n1,2,nan,nan,1,0\n2,1,nan,nan,1,0\n')
        out = root / 'summary.json'
        command = [sys.executable, str(ROOT / 'tools/coldflow_characteristics.py'),
                   '--mesh', str(mesh), '--field', str(field), '--json', str(out)]
        run = subprocess.run(command, capture_output=True, text=True)
        if run.returncode:
            raise AssertionError(run.stderr)
        result = json.loads(out.read_text())
        assert result['two_phase']['cells'] == 1, result['two_phase']
        assert result['two_phase']['board_downstream_cells'] == 0, result['two_phase']
        assert result['two_phase']['post_tip_cells'] == 0, result['two_phase']
        assert result['eos_liquid']['cells'] == 1, result['eos_liquid']
        assert result['cryogenic_or_dense']['cells'] == 2, result['cryogenic_or_dense']
        assert result['streams']['dense_cells'] == 2
        assert result['streams']['pure_o2_post_tip_cells'] == 0
        assert 'liquid_cells' not in result['streams']
        csv.unlink()
        subprocess.run(command, check=True, capture_output=True, text=True)
        result = json.loads(out.read_text())
        assert result['two_phase']['cells'] is None, result['two_phase']
        assert result['cryogenic_or_dense']['cells'] == 2
        csv.write_text('cell,region\n1,5\n0,2\n2,1\n')
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode != 0 and 'phase sidecar/mesh cell order differs' in result.stderr
        periodic=Path(str(field)+'.s90')
        periodic.write_bytes(field.read_bytes())
        phase_periodic=Path(str(periodic)+'.phase.csv')
        phase_periodic.write_text('cell,region,vapor_mass_fraction,vapor_volume_fraction,acoustic_evidence,projected\n'
                                  '0,5,0.3,0.6,2,0\n1,2,nan,nan,1,0\n2,1,nan,nan,1,0\n')
        periodic_cmd=[sys.executable,str(ROOT/'tools/coldflow_characteristics.py'),
                      '--mesh',str(mesh),'--field',str(periodic),'--json',str(out)]
        result=subprocess.run(periodic_cmd,capture_output=True,text=True)
        assert result.returncode!=0 and 'matching step metadata' in result.stderr
        status=Path(str(periodic)+'.status.json')
        status.write_text(json.dumps({'steps':91,'phase_sidecar_available':True}))
        result=subprocess.run(periodic_cmd,capture_output=True,text=True)
        assert result.returncode!=0 and 'phase/status step mismatch' in result.stderr
        status.write_text(json.dumps({'steps':90,'phase_sidecar_available':True}))
        subprocess.run(periodic_cmd,check=True,capture_output=True,text=True)
        assert json.loads(out.read_text())['two_phase']['cells']==1
    print('COLD_FLOW_PHASE_PROVENANCE_PASS')


if __name__ == '__main__':
    main()
