#!/usr/bin/env python3
"""Convert an RFAMTETF1 cell field dump plus its mesh into a VTU file.

The dump is cell-centred (no duplicated point data), so ParaView shows one
uniform value per tetrahedron. For a smooth surface display use
Filters -> Cell Data to Point Data in ParaView.
"""
from __future__ import annotations
import argparse
import struct
import json
import re
from pathlib import Path
import numpy as np
import meshio

MAGIC = b'RFAMTETF1'
FIELD_NAMES = ['p', 'T', 'u', 'v', 'w', 'rho', 'sound', 'h']


def verify_phase_pair(field_path, phase_path):
    """A periodic phase CSV is usable only with same-step dump metadata."""
    if not Path(phase_path).exists():
        return False
    match = re.search(r'\.s(\d+)$', str(field_path))
    if match:
        status_path=Path(str(field_path)+'.status.json')
        if not status_path.exists():
            raise ValueError('periodic EOS phase requires matching step metadata')
        status=json.loads(status_path.read_text(encoding='utf-8'))
        if status.get('steps')!=int(match.group(1)) or status.get('phase_sidecar_available') is not True:
            raise ValueError('periodic EOS phase/status step mismatch')
    return True


def read_mesh(path):
    raw = Path(path).read_bytes()
    if raw[:8] != b'RFAMTET1':
        raise ValueError('not a RFAMTET1 mesh: ' + str(path))
    offset = 8
    npoints, ncells, nfaces = struct.unpack_from('<QQQ', raw, offset)
    offset += 24
    points = np.frombuffer(raw, '<f8', int(npoints) * 3, offset).reshape(-1, 3)
    offset += int(npoints) * 24
    cells = np.frombuffer(raw, '<i4', int(ncells) * 4, offset).reshape(-1, 4)
    offset += int(ncells) * 16
    volumes = np.frombuffer(raw, '<f8', int(ncells), offset)
    return points, cells, volumes


def read_field(path):
    raw = Path(path).read_bytes()
    if raw[:9] != MAGIC:
        raise ValueError('not a RFAMTETF1 dump: ' + str(path))
    offset = 9
    ncells, nspecies = struct.unpack_from('<QQ', raw, offset)
    offset += 16
    ncells, nspecies = int(ncells), int(nspecies)
    centroids = np.frombuffer(raw, '<f4', ncells * 3, offset).reshape(-1, 3)
    offset += ncells * 12
    fields = np.frombuffer(raw, '<f4', ncells * 8, offset).reshape(-1, 8)
    offset += ncells * 32
    composition = np.frombuffer(raw, '<f4', ncells * nspecies, offset).reshape(ncells, nspecies)
    return centroids, fields, composition, nspecies


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', required=True)
    parser.add_argument('--field', required=True, help='RFAMTETF1 dump')
    parser.add_argument('--out', required=True, help='output .vtu')
    parser.add_argument('--species', default='O2,C3H8')
    args = parser.parse_args()
    points, cells, volumes = read_mesh(args.mesh)
    centroids, fields, composition, nspecies = read_field(args.field)
    if len(cells) != len(fields):
        raise ValueError('mesh and field cell counts differ')
    names = [s.strip() for s in args.species.split(',')]
    if len(names) != nspecies:
        raise ValueError(f'--species expects {nspecies} names')
    mach = np.sqrt(fields[:, 2]**2 + fields[:, 3]**2 + fields[:, 4]**2) / np.maximum(fields[:, 6], 1e-30)
    cell_data = {
        'p': [fields[:, 0]],
        'T': [fields[:, 1]],
        'velocity': [fields[:, 2:5]],
        'rho': [fields[:, 5]],
        'sound': [fields[:, 6]],
        'enthalpy': [fields[:, 7]],
        'Mach': [mach],
        'volume_m3': [volumes.astype(np.float32)],
    }
    for index, name in enumerate(names):
        cell_data[name] = [composition[:, index]]
    status_path=Path(args.field+'.status.json')
    status=json.loads(status_path.read_text()) if status_path.exists() else {'converged':False,'production_validated':False}
    cell_data['converged']=[np.full(len(cells),int(status.get('converged',False)),dtype=np.uint8)]
    cell_data['production_validated']=[np.full(len(cells),int(status.get('production_validated',False)),dtype=np.uint8)]
    # -1 means not audited; never infer boundary compatibility from a completed run.
    boundary_valid = int(status['boundary_state_valid']) if 'boundary_state_valid' in status else -1
    cell_data['boundary_state_valid']=[np.full(len(cells),boundary_valid,dtype=np.int8)]
    phase_path=Path(args.field+'.phase.csv')
    if phase_path.exists():
        verify_phase_pair(args.field,phase_path)
        phase=np.genfromtxt(phase_path,delimiter=',',names=True)
        if len(phase)!=len(cells) or not np.array_equal(phase['cell'],np.arange(len(cells))):
            raise ValueError('phase sidecar and field cell order differ')
        for name in ['region','vapor_mass_fraction','vapor_volume_fraction','acoustic_evidence','projected']:
            cell_data[name]=[phase[name]]
    meshio.write(args.out, meshio.Mesh(points, [('tetra', cells)], cell_data=cell_data))
    speed = np.linalg.norm(fields[:, 2:5], axis=1)
    print(f'VTU_WRITTEN {args.out} cells={len(cells)}')
    print(f'  p  [{fields[:,0].min():.4g}, {fields[:,0].max():.4g}] Pa')
    print(f'  T  [{fields[:,1].min():.4g}, {fields[:,1].max():.4g}] K')
    print(f'  |u|[{speed.min():.4g}, {speed.max():.4g}] m/s   Mach_max={mach.max():.4g}')
    print(f'  rho[{fields[:,5].min():.4g}, {fields[:,5].max():.4g}] kg/m3')
    for index, name in enumerate(names):
        print(f'  Y_{name} [{composition[:,index].min():.4g}, {composition[:,index].max():.4g}]')


if __name__ == '__main__':
    main()
