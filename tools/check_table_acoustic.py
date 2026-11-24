#!/usr/bin/env python3
"""Independent check of the RFAMTAB2 acoustic column (a2) against the direct
CoolProp equilibrium sound speed printed by verify_property_table."""
import struct
import sys
from pathlib import Path

import numpy as np

NODE = np.dtype([('rho', '<f8'), ('h', '<f8'), ('e', '<f8'), ('s', '<f8'),
                 ('a2', '<f8'), ('phase', '<i4'), ('pad', '<i4')])  # sizeof=48 (alignment)


def read_table(path):
    raw = Path(path).read_bytes()
    assert raw[:8] == b'RFAMTAB2', 'bad magic'
    ny, nt, np_ = struct.unpack_from('<QQQ', raw, 8)
    off = 32
    Y = np.frombuffer(raw, '<f8', ny, off); off += 8*ny
    T = np.frombuffer(raw, '<f8', nt, off); off += 8*nt
    P = np.frombuffer(raw, '<f8', np_, off); off += 8*np_
    nodes = np.frombuffer(raw, NODE, ny*nt*np_, off).reshape(ny, nt, np_)
    return Y, T, P, nodes


def bracket(axis, v):
    return int(np.clip(np.searchsorted(axis, v) - 1, 0, len(axis) - 2))


def interp_a2(Y, T, P, nodes, yO2, t, p):
    iy, it, ip = bracket(Y, yO2), bracket(T, t), bracket(P, p)
    fy = (yO2 - Y[iy]) / (Y[iy+1] - Y[iy])
    ft = (t - T[it]) / (T[it+1] - T[it])
    fp = (p - P[ip]) / (P[ip+1] - P[ip])
    corners = [(iy+a, it+b, ip+c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
    vals = []
    for (a, b, c) in corners:
        n = nodes[a, b, c]
        if n['phase'] < 0 or n['a2'] <= 0:
            return None
        w = (fy if (a == iy+1) else 1-fy) * (ft if (b == it+1) else 1-ft) * (fp if (c == ip+1) else 1-fp)
        vals.append(w * n['a2'])
    return float(np.sum(vals))


def main():
    table = sys.argv[1] if len(sys.argv) > 1 else 'work/a406_rootcause/property_table_op.bin'
    # (Y, T, p, direct CoolProp equilibrium sound speed from the verification run)
    refs = [(0.99, 150.0, 4.0e6, 131.623), (0.1, 250.0, 4.0e6, 143.222),
            (0.9, 170.0, 4.0e6, 199.168), (0.73, 250.0, 4.0e6, 238.448),
            (0.9162460633, 189.802, 4.09773e6, 221.16)]
    Y, T, P, nodes = read_table(table)
    print('table', table, 'shape', nodes.shape, 'Y', Y[:3], '...', Y[-2:])
    print('phase counts', {int(k): int(v) for k, v in
          zip(*np.unique(nodes['phase'], return_counts=True))})
    print('a2 valid fraction', float((nodes['a2'] > 0).mean()))
    inside = 0; worst = 0.
    for y, t, p, direct in refs:
        if not (Y[0] <= y <= Y[-1] and T[0] <= t <= T[-1] and P[0] <= p <= P[-1]):
            print('REF Y=%.6g T=%.3f p=%.6g -> out of table range' % (y, t, p))
            continue
        inside += 1
        a2 = interp_a2(Y, T, P, nodes, y, t, p)
        if a2 is None:
            print('REF Y=%.6g T=%.3f p=%.6g -> no valid a2 corners' % (y, t, p))
            continue
        a = float(np.sqrt(a2)); rel = abs(a-direct)/direct
        worst = max(worst, rel)
        print('REF Y=%.6g T=%.3f p=%.6g -> table a=%.6g direct a=%.6g rel=%.2e'
              % (y, t, p, a, direct, rel))
    print('in-range refs %d, worst relative error %.2e' % (inside, worst))


if __name__ == '__main__':
    main()
