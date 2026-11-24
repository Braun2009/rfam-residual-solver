#!/usr/bin/env python3
"""Central-zone kinetic-motion analysis for the 3-pintle A-route field.

The volume sums of rho*u_i*u_j are kinetic second moments (units J),
not momentum fluxes [N]. The region starts at the board and includes the
protruding internal torch passage. Its radial/axial ratio does not prove vector
momentum cancellation at the axis. That requires signed momentum fluxes
through an enclosing control surface and a converged balance.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tetra_field_to_vtu import read_mesh, read_field  # noqa: E402
from coldflow_characteristics import fit_azimuthal_modes  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh', required=True)
    ap.add_argument('--field', required=True)
    ap.add_argument('--board-z', type=float, default=0.2709527893066406)
    ap.add_argument('--axis-x', type=float, default=0.04921534822760748)
    ap.add_argument('--axis-y', type=float, default=0.04835004867381059)
    ap.add_argument('--radius-mm', type=float, default=5.0)
    ap.add_argument('--depth-mm', type=float, default=20.0)
    ap.add_argument('--label', default='')
    args = ap.parse_args()

    _, cells, vol = read_mesh(args.mesh)
    cent, f, y, _ = read_field(args.field)
    p, T, u, v, w, rho, a, _ = (f[:, i] for i in range(8))
    x = cent[:, 0] - args.axis_x
    yy = cent[:, 1] - args.axis_y
    r = np.hypot(x, yy)
    z = cent[:, 2]
    ur = (u * x + v * yy) / np.maximum(r, 1e-12)
    ut = (-u * yy + v * x) / np.maximum(r, 1e-12)
    sel = (r < args.radius_mm * 1e-3) & (z > args.board_z) & (z < args.board_z + args.depth_mm * 1e-3)
    n = int(sel.sum())
    print(f'== board-downstream core diagnostic (includes torch passage) {args.label} ==')
    print(f'zone: r<{args.radius_mm} mm, face..face+{args.depth_mm} mm  ->  cells={n} '
          f'volume={vol[sel].sum()*1e6:.4g} cm3 mass={np.sum(rho[sel]*vol[sel])*1e6:.4g} g')
    if n == 0:
        print('empty board-core selection')
        return 0
    V = vol[sel]
    def M(a_i, a_j):
        return float(np.sum(rho[sel] * a_i[sel] * a_j[sel] * V))
    Mzz, Mrr, Mtt = M(w, w), M(ur, ur), M(ut, ut)
    Mzr, Mzt, Mrt = M(w, ur), M(w, ut), M(ur, ut)
    print(f'volume-integrated kinetic second moments [J]:  Mzz={Mzz:9.3f}  Mrr={Mrr:9.3f}  Mtt={Mtt:9.3f}')
    print(f'                           Mzr={Mzr:9.3f}  Mzt={Mzt:9.3f}  Mrt={Mrt:9.3f}')
    print(f'radial/axial motion ratio (NOT vector cancellation): Mrr/Mzz={Mrr/max(Mzz,1e-30):.3f}  '

          f'off-diagonal magnitude relative to Mzz: '
          f'{max(abs(Mzr),abs(Mzt))/(max(Mzz,1e-30)):.3f}')
    # Azimuthal Fourier modes of the radial kinetic second moment.
    th = np.arctan2(yy[sel], x[sel] - 0 * x[sel])
    q = rho[sel] * ur[sel] * ur[sel] * V          # radial momentum-flux density
    c1 = float(np.sum(q * np.cos(th))); s1 = float(np.sum(q * np.sin(th)))
    c3 = float(np.sum(q * np.cos(3 * th))); s3 = float(np.sum(q * np.sin(3 * th)))
    print(f'azimuthal modes of radial kinetic moment [J]: |m=1|={np.hypot(c1,s1):.4f}  '
          f'|m=3|={np.hypot(c3,s3):.4f} (no signed momentum balance)')
    spd = np.sqrt(u * u + v * v + w * w)
    print(f'speed in zone [m/s]: mean={spd[sel].mean():.2f} std={spd[sel].std():.2f} '
          f'min={spd[sel].min():.2f} max={spd[sel].max():.2f}')
    print(f'Y_O2 in zone: mean={y[sel,0].mean():.4f} std={y[sel,0].std():.4f} '
          f'min={y[sel,0].min():.4f} max={y[sel,0].max():.4f}')
    print(f'T in zone  [K]: mean={T[sel].mean():.1f} std={T[sel].std():.1f} '
          f'min={T[sel].min():.1f} max={T[sel].max():.1f}')
    print(f'p in zone [MPa]: mean={p[sel].mean()/1e6:.4f} std={p[sel].std()/1e6:.4f}')
    # Annulus fit includes an intercept to remove nonuniform mesh-weight bias.
    # Phase modulo 120 degrees can be compared to geometry, but cannot alone
    # establish a steady CFD design verdict.
    ann = (z > args.board_z) & (z < args.board_z + args.depth_mm * 1e-3)         & (r > 0.008) & (r < 0.018)
    if ann.sum() > 100:
        tt = np.arctan2(yy[ann], x[ann])
        wgt = rho[ann] * vol[ann]
        print('board-downstream annulus (includes injector passages, r 8-18 mm, first %.0f mm):' % args.depth_mm)
        for label, field in (('Y_O2', y[ann, 0]), ('u_z', w[ann]), ('T', T[ann]), ('p', p[ann])):
            modes = fit_azimuthal_modes(tt,wgt,field)
            print(f"  {label:5s} mean={modes['mean']:12.5g} |m=1|={modes['m1']:10.4g} "
                  f"|m=3|={modes['m3']:10.4g} (m=3/mean {modes['m3_over_mean']:.4f}, "
                  f"phase {modes['phase_deg']:6.1f} deg; method={modes['method']})")
    # axial velocity profile vs radius (the 'uniform core' shape)
    print('radial profile of the zone:')
    for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5)):
        m = sel & (r >= lo * 1e-3) & (r < hi * 1e-3)
        if m.sum() == 0:
            continue
        print(f'  r {lo}-{hi} mm: cells={m.sum():6d} |w|={abs(w[m]).mean():7.2f} m/s '
              f'url={np.sqrt((ur[m]**2).mean()):7.2f} Y_O2={y[m,0].mean():.4f} T={T[m].mean():.1f} K')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
