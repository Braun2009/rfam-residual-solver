#!/usr/bin/env python3
"""Tetrahedral mesh audit for the A406 fluid grid.

Reads the RFAM binary .rfamtet directly and reports:
  * volume distribution and sliver population
  * per-cell shape quality (radius ratio r_in/r_circ, and edge-ratio)
  * boundary faces by zone (count, area) -- the boundary-condition census
  * 120 degree rotational symmetry about the chamber axis, which the 3-post
    LOX pattern requires. A mesh that fails this cannot support any
    symmetry-based reduction, and a field that fails it means the solution
    broke symmetry numerically.

Usage: python mesh_audit.py <mesh.rfamtet> [--axis 0.0492153,0.04835]
"""
import sys
import struct
import math
import numpy as np

def read_mesh(path):
    with open(path, 'rb') as f:
        magic = f.read(8)
        if magic[:8] != b'RFAMTET1':
            raise SystemExit('bad magic %r' % magic)
        np_, nc, nf = struct.unpack('<QQQ', f.read(24))
        pts = np.frombuffer(f.read(np_ * 24), dtype='<f8').reshape(np_, 3)
        cells = np.frombuffer(f.read(nc * 16), dtype='<i4').reshape(nc, 4)
        vols = np.frombuffer(f.read(nc * 8), dtype='<f8')
        raw = f.read(nf * 28)
        faces = np.frombuffer(raw, dtype=np.dtype([('n', '<i4', 3), ('own', '<i4'),
                                                   ('nb', '<i4'), ('zone', '<i4')]))
    return pts, cells, vols, faces

def cell_quality(pts, cells):
    """Radius ratio 3*r_in/r_circ in [0,1]; 1 is a regular tetrahedron."""
    p = pts[cells]                                  # (nc,4,3)
    e = [p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0],
         p[:, 2] - p[:, 1], p[:, 3] - p[:, 1], p[:, 3] - p[:, 2]]
    edges = np.stack([np.linalg.norm(v, axis=1) for v in e], axis=1)
    emax = edges.max(axis=1)
    emin = edges.min(axis=1)
    vol = np.abs(np.einsum('ij,ij->i', e[0], np.cross(e[1], e[2]))) / 6.0
    area_sum = 0.0
    for a, b in ((1, 2), (1, 3), (2, 3), (0, 2), (0, 3), (0, 1)):
        pass
    # surface area
    surf = np.zeros(len(cells))
    tri = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))
    for i, j, k in tri:
        surf += 0.5 * np.linalg.norm(np.cross(p[:, j] - p[:, i], p[:, k] - p[:, i]), axis=1)
    # Normalised volume ratio: 6*sqrt(2)*V/a_max^3, exactly 1 for a regular
    # tetrahedron of edge a_max and -> 0 as the cell degenerates.
    q = 6.0 * math.sqrt(2.0) * vol / emax ** 3
    return q, emax / np.maximum(emin, 1e-300), vol

def symmetry(pts, axis, period_deg, tol_scale=1.0):
    """Rotate all points about the z axis and measure the nearest-neighbour gap.

    Uses a uniform spatial hash sized from the typical edge so the search is
    O(N) instead of O(N^2); scipy is not available on the run hosts.
    """
    ang = math.radians(360.0 / period_deg)
    c, s = math.cos(ang), math.sin(ang)
    dx = pts[:, 0] - axis[0]
    dy = pts[:, 1] - axis[1]
    rot = np.stack([axis[0] + c * dx - s * dy,
                    axis[1] + s * dx + c * dy,
                    pts[:, 2]], axis=1)

    edge = float(np.median(np.linalg.norm(pts[1:] - pts[:-1], axis=1)))
    h = max(edge * 2.0, 1e-9)
    grid = {}
    keys = np.floor(pts / h).astype(np.int64)
    for i, k in enumerate(map(tuple, keys)):
        grid.setdefault(k, []).append(i)

    gaps = np.empty(0)
    idx = np.linspace(0, len(rot) - 1, num=min(4000, len(rot))).astype(int)
    gaps = np.empty(len(idx))
    maxR = 3
    for n, i in enumerate(idx):
        q = rot[i]
        base = np.floor(q / h).astype(np.int64)
        best = float('inf')
        for r in range(maxR + 1):
            cand = []
            for ix in range(base[0] - r, base[0] + r + 1):
                for iy in range(base[1] - r, base[1] + r + 1):
                    for iz in range(base[2] - r, base[2] + r + 1):
                        v = grid.get((ix, iy, iz))
                        if v:
                            cand.extend(v)
            if cand:
                d = np.linalg.norm(pts[cand] - q, axis=1)
                best = float(d.min())
                if best <= (r - 1) * h:
                    break
        gaps[n] = best if best < float('inf') else (maxR + 1) * h
    return gaps

def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    path = sys.argv[1]
    axis = (0.0492153, 0.04835)
    for i, a in enumerate(sys.argv):
        if a == '--axis':
            axis = tuple(float(v) for v in sys.argv[i + 1].split(','))
    pts, cells, vols, faces = read_mesh(path)
    print('mesh %s' % path)
    print('points=%d cells=%d faces=%d' % (len(pts), len(cells), len(faces)))

    q, ar, vol = cell_quality(pts, cells)
    print('\n-- volume [m3]')
    for p in (0, 1, 50, 99, 100):
        print('  p%-3d %.6e' % (p, np.percentile(vol, p)))
    med = np.median(vol)
    print('  below 1e-3*median: %d cells (%.3e of total volume)'
          % ((vol < 1e-3 * med).sum(), vol[vol < 1e-3 * med].sum() / vol.sum()))

    print('\n-- shape quality (6 sqrt2 V / a_max^3, 1.0 = regular)')
    for p in (0, 1, 5, 50):
        print('  p%-3d %.4f' % (p, np.percentile(q, p)))
    print('  q<0.01: %d   q<0.05: %d   q<0.1: %d'
          % ((q < 0.01).sum(), (q < 0.05).sum(), (q < 0.1).sum()))
    print('  edge ratio p50=%.1f p99=%.1f max=%.1f'
          % (np.percentile(ar, 50), np.percentile(ar, 99), ar.max()))

    print('\n-- boundary census by zone')
    bnd = faces[faces['nb'] < 0]
    p3 = pts[bnd['n']]
    area = 0.5 * np.linalg.norm(np.cross(p3[:, 1] - p3[:, 0], p3[:, 2] - p3[:, 0]), axis=1)
    for z in np.unique(bnd['zone']):
        m = bnd['zone'] == z
        print('  zone %-3d faces %-7d area %.6e m2' % (z, m.sum(), area[m].sum()))

    typical = np.median(np.linalg.norm(pts[cells[:, 1]] - pts[cells[:, 0]], axis=1))
    print('\n-- rotational symmetry about (%.7f,%.7f), typical edge %.3e m'
          % (axis[0], axis[1], typical))
    best = None
    for period in (2, 3, 4, 5, 6, 8, 12):
        g = symmetry(pts, axis, period)
        r = np.percentile(g, 99) / typical
        print('  period %-2d gap/edge p50=%.2f p99=%.2f max=%.2f'
              % (period, np.percentile(g, 50) / typical, r, g.max() / typical))
        if best is None or r < best[1]:
            best = (period, r)
    print('  best period %d (p99/edge = %.2f)' % best)

if __name__ == '__main__':
    main()
