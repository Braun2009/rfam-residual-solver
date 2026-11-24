#!/usr/bin/env python3
"""Cold-flow probe for the real-geometry (A route) tetra dumps.

Prints board-downstream diagnostics (which include internal passages) and
post-tip diagnostics. Threshold front positions are not physical-time jet
penetration until the actual outlets and a converged flow field are verified.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tetra_field_to_vtu import read_mesh, read_field  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh', required=True)
    ap.add_argument('--field', required=True)
    ap.add_argument('--board-z', type=float, default=0.2709527893066406)
    ap.add_argument('--injector-tip-z', type=float, default=0.2858028)
    ap.add_argument('--axis-x', type=float, default=0.04921534822760748)
    ap.add_argument('--axis-y', type=float, default=0.04835004867381059)
    ap.add_argument('--label', default='')
    args = ap.parse_args()

    _, cells, vol = read_mesh(args.mesh)
    cent, f, y, _ = read_field(args.field)
    p, T, u, v, w, rho, a, _ = (f[:, i] for i in range(8))
    spd = np.sqrt(u*u + v*v + w*w)
    z = cent[:, 2]
    r = np.hypot(cent[:, 0] - args.axis_x, cent[:, 1] - args.axis_y)
    inner = z < args.board_z
    downstream = ~inner  # Internal injector passages protrude above the board.
    post_tip = z >= args.injector_tip_z
    dense = rho > 500.
    o2 = y[:, 0] > 0.95
    fuel = y[:, 0] < 0.2

    print(f'== A-route cold-flow probe {args.label} ==')
    print(f'cells={len(cells)} upstream_board={inner.sum()} ({(vol[inner].sum()*1e6):.4g} cm3) '
          f'board_downstream_including_passages={downstream.sum()} ({(vol[downstream].sum()*1e6):.4g} cm3) '
          f'certain_post_tip={post_tip.sum()} ({(vol[post_tip].sum()*1e6):.4g} cm3)')
    print(f'p [{p.min()/1e6:.4g},{p.max()/1e6:.4g}] MPa  T [{T.min():.5g},{T.max():.5g}] K '
          f'|u| [{spd.min():.4g},{spd.max():.4g}] m/s  Ma_max {max(spd/a):.4g}')
    for name, mask in (('density proxy rho>500, NOT EOS liquid', dense), ('pure O2 Y>0.95', o2), ('fuel Y<0.2', fuel)):
        if mask.any():
            zin, zdown = z[mask & inner], z[mask & downstream]
            print(f'{name}: cells={mask.sum()} mass={float(np.sum(rho[mask]*vol[mask])):.6g} kg '
                  f'upstream_board_cells={len(zin)} board_downstream_cells={len(zdown)} '
                  f'post_tip_cells={int((mask & post_tip).sum())}'
                  + (f' z_upstream [{zin.min():.5g},{zin.max():.5g}]' if len(zin) else ''))
        else:
            print(f'{name}: none')
    # The upper ends of channels are not the board plane.
    for name, mask in (('LOX', o2), ('fuel', fuel)):
        if mask.any():
            print(f'{name} furthest threshold centroid z={z[mask].max():.6g}; '
                  f'CAD tip z={args.injector_tip_z:.6g} (not a physical-time front)')
    slab = downstream & (z < args.board_z + 2.0e-3)
    if slab.any():
        edges = np.array([0, 2, 4, 6, 8, 10, 14, 18, 24, 40]) * 1e-3
        print('board-downstream annulus (INCLUDES internal passages, 2 mm from board):')
        print('  r_mm      cells   p_MPa    T_K     uz_m_s   Y_O2')
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = slab & (r >= lo) & (r < hi)
            if not m.any():
                continue
            print(f'  {1000*(lo+hi)/2:6.1f} {m.sum():7d} {p[m].mean()/1e6:8.4f} '
                  f'{T[m].mean():7.1f} {w[m].mean():9.3f} {y[m,0].mean():8.4f}')
    for label,region in (('board-downstream including passages',downstream),('post-tip',post_tip)):
        if not region.any():continue
        sc = np.where(region)[0]
        top = sc[np.argsort(spd[sc])[-5:][::-1]]
        print(f'fastest {label} cells:')
        for i in top:
            print(f'  z={z[i]:.5g} r={1000*r[i]:.3f} mm |u|={spd[i]:.4g} Ma={spd[i]/a[i]:.4g} '
                  f'p={p[i]/1e6:.4g} MPa T={T[i]:.4g} K rho={rho[i]:.4g} Y_O2={y[i,0]:.4g}')
        volume = vol[region].sum()
        print(f'{label} backflow (w<0) volume fraction: {100.*vol[region][w[region]<0].sum()/volume:.4g}%')
        for thr in (5., 20., 50.):
            print(f'{label} |u|<{thr:g} m/s volume fraction: '
                  f'{100.*vol[region][spd[region]<thr].sum()/volume:.4g}%')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
