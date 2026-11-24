#!/usr/bin/env python3
"""Prime a cold-flow field with the equivalent-injector mixture near the inlet.

Initial-condition tool only: it does not change the boundary conditions or the
solver. The B-route equivalent injector delivers the propellant as a mixed
stream (LOX 2.8 kg/s at 130 K + C3H8 1.05 kg/s at 450 K); at 4.112 MPa that
mixture is T=235.53 K, rho=93.04 kg/m3, Q=0.824, Y_O2=0.72727 (mass weighted).
Cells within a small distance of the inlet faces are set to that state with the
axial injection velocity so the jets start developed instead of crawling
through the sub-0.3 mm plate refinement for thousands of pseudo-time steps.
"""
import argparse
import struct
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ZP = 0.2709527893066406 + 0.014          # inlet plane (board + 14 mm)
AX = (0.04921534822760748, 0.04835004867381059)


def read_mesh(path):
    raw = Path(path).read_bytes()
    np_, nc, nf = struct.unpack_from('<QQQ', raw, 8)
    points = np.frombuffer(raw, '<f8', np_*3, 32).reshape(-1, 3)
    face_offset = 32 + np_*24 + nc*24
    faces = np.frombuffer(raw, '<i4', nf*6, face_offset).reshape(-1, 6)
    return points, faces


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh', required=True)
    ap.add_argument('--field', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--radius', type=float, default=2.5e-3,
                    help='distance to the nearest inlet face centroid [m]')
    ap.add_argument('--depth', type=float, default=6e-3,
                    help='also prime cells within this distance of the inlet plane [m]')
    ap.add_argument('--footprint', type=float, default=20e-3,
                    help='radial footprint of the depth criterion [m]')
    ap.add_argument('--T', type=float, default=235.5262)
    ap.add_argument('--rho', type=float, default=93.0440)
    ap.add_argument('--yO2', type=float, default=0.7272727)
    ap.add_argument('--uaxial', type=float, default=63.45)
    ap.add_argument('--table', default='',
                    help='RFAMTAB2 table used to store consistent rho/h at each cell pressure')
    args = ap.parse_args()

    points, faces = read_mesh(args.mesh)
    raw = Path(args.field).read_bytes()
    assert raw[:9] == b'RFAMTETF1'
    ncells, nsp = struct.unpack_from('<QQ', raw, 9)
    off = 25
    centroids = np.frombuffer(raw, '<f4', ncells*3, off).reshape(-1, 3).astype(np.float64)
    off += ncells*12
    fields = np.frombuffer(raw, '<f4', ncells*8, off).reshape(-1, 8).copy()
    off += ncells*32
    composition = np.frombuffer(raw, '<f4', ncells*nsp, off).reshape(-1, nsp).copy()
    assert nsp == 2

    inlets = np.vstack([points[faces[faces[:, 5] == z][:, :3]].mean(axis=1)
                        for z in (10, 11, 12, 20, 21, 22)])
    tree = cKDTree(inlets)
    distance, _ = tree.query(centroids)
    zmm = (centroids[:, 2] - ZP) * 1e3
    radial = np.hypot(centroids[:, 0]-AX[0], centroids[:, 1]-AX[1])
    selected = (distance < args.radius) | ((zmm >= -1.) & (zmm < args.depth*1e3) & (radial < args.footprint))
    print('PRIME selected', int(selected.sum()), 'of', ncells)
    fields[selected, 1] = args.T
    fields[selected, 2] = 0.0
    fields[selected, 3] = 0.0
    fields[selected, 4] = args.uaxial
    composition[selected, 0] = args.yO2
    composition[selected, 1] = 1.0 - args.yO2
    # Direct CoolProp PT evaluation at each selected cell's own pressure keeps
    # the stored rho/h consistent with the solver's seed check (2e-5 relative).
    try:
        from CoolProp.CoolProp import AbstractState
        from CoolProp import CoolProp as CP
        href = {}
        for name in ('Oxygen', 'Propane'):
            st = AbstractState('HEOS', name)
            st.update(CP.DmassT_INPUTS, 1e-10, 298.15)
            href[name] = st.hmass_idealgas()
        offset = args.yO2*href['Oxygen'] + (1.0-args.yO2)*href['Propane']
        mix = AbstractState('HEOS', 'Oxygen&Propane')
        mix.set_mass_fractions([args.yO2, 1.0-args.yO2])
        for c in np.where(selected)[0]:
            pv = float(fields[c, 0])
            mix.update(CP.PT_INPUTS, pv, args.T)
            fields[c, 5] = mix.rhomass()
            fields[c, 6] = mix.umass() - offset
            fields[c, 7] = mix.hmass() - offset
        print('PRIME CoolProp-consistent rho/e/h stored')
    except Exception as exc:
        raise SystemExit('priming needs CoolProp: ' + str(exc))
    if False:
        import struct as _struct
        raw = Path(args.table).read_bytes()
        assert raw[:8] == b'RFAMTAB2'
        ny, nt, np_ = _struct.unpack_from('<QQQ', raw, 8)
        off = 32
        Y = np.frombuffer(raw, '<f8', ny, off); off += 8*ny
        T = np.frombuffer(raw, '<f8', nt, off); off += 8*nt
        P = np.frombuffer(raw, '<f8', np_, off); off += 8*np_
        node = np.dtype([('rho','<f8'),('h','<f8'),('e','<f8'),('s','<f8'),
                         ('a2','<f8'),('phase','<i4'),('pad','<i4')])
        nodes = np.frombuffer(raw, node, ny*nt*np_, off).reshape(ny, nt, np_)
        def bracket(axis, v):
            return int(np.clip(np.searchsorted(axis, v)-1, 0, len(axis)-2))
        def interp(yv, tv, pv, key):
            iy, it, ip = bracket(Y,yv), bracket(T,tv), bracket(P,pv)
            fy=(yv-Y[iy])/(Y[iy+1]-Y[iy]); ft=(tv-T[it])/(T[it+1]-T[it]); fp=(pv-P[ip])/(P[ip+1]-P[ip])
            v=0.
            for a in (0,1):
                for b in (0,1):
                    for c in (0,1):
                        w=(fy if a else 1-fy)*(ft if b else 1-ft)*(fp if c else 1-fp)
                        v+=w*nodes[iy+a,it+b,ip+c][key]
            return v
        for c in np.where(selected)[0]:
            pv=float(fields[c,0])
            fields[c,5]=interp(args.yO2, args.T, pv, 'rho')
            fields[c,6]=interp(args.yO2, args.T, pv, 'e')
            fields[c,7]=interp(args.yO2, args.T, pv, 'h')
        print('PRIME table-consistent rho/h stored')
    out = Path(args.out)
    with out.open('wb') as fh:
        fh.write(b'RFAMTETF1')
        fh.write(struct.pack('<QQ', ncells, nsp))
        fh.write(centroids.astype('<f4').tobytes())
        fh.write(fields.astype('<f4').tobytes())
        fh.write(composition.astype('<f4').tobytes())
    print('PRIME wrote', out)


if __name__ == '__main__':
    main()
