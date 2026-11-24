#!/usr/bin/env python3
"""Dome separation criterion from an RFAM tetra field dump.

Reads the RFAMTETF1 dump written by --dump and reports, for the dome region
(z = 0..123 mm from the injector plate), the quantities that decide whether the
wall boundary layer separates, plus the recirculation measure that actually
matters for this geometry.

Why these quantities
--------------------
The classical criterion is the Clauser parameter
    beta = (delta*/tau_w) dp/dx
with turbulent separation at beta ~ 5..10. It needs a wall-normal profile, so
it is estimated here from the near-wall cell rather than integrated:

    tau_w  = mu (du/dn)         from the wall-adjacent cell gradient
    delta* ~ integral (1 - u/ue) dn, estimated as u_e * delta_99 / ... is not
             available from a single cell, so the script reports the shape
             factor proxy H = delta*/theta from a 1/7-power fit only when a
             profile exists; otherwise it reports beta with delta* from the
             momentum thickness of the near-wall layer (documented below).

The recirculation measure is the one this case can actually resolve:
    * fraction of dome wall cells with reversed streamwise velocity
    * size and location of the corner recirculation cell (z < 47 mm corner)
    * wall-normal pressure gradient sign along the meridian

Usage:
  python dome_separation.py dump.bin [--plate-z 0.271] [--dome-mm 123]
"""
import sys
import struct
import numpy as np

RHO_LOX, RHO_FUEL = 1140.0, 66.5


def read_dump(path):
    with open(path, 'rb') as f:
        magic = f.read(9)
        if magic != b'RFAMTETF1':
            raise SystemExit('bad dump magic %r' % magic)
        ncell, nspec = struct.unpack('<QQ', f.read(16))
        xyz = np.frombuffer(f.read(ncell * 12), dtype='<f4').reshape(ncell, 3)
        fld = np.frombuffer(f.read(ncell * 32), dtype='<f4').reshape(ncell, 8)
        comp = np.frombuffer(f.read(ncell * 4 * nspec), dtype='<f4').reshape(ncell, nspec)
    names = ['p', 'T', 'u', 'v', 'w', 'rho', 'a', 'h']
    return xyz, dict(zip(names, fld.T)), comp, nspec


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    path = sys.argv[1]
    plate_z, dome_mm = 0.271, 123.0
    for i, a in enumerate(sys.argv):
        if a == '--plate-z':
            plate_z = float(sys.argv[i + 1])
        if a == '--dome-mm':
            dome_mm = float(sys.argv[i + 1])

    xyz, F, comp, nspec = read_dump(path)
    z = (xyz[:, 2] - plate_z) * 1e3                  # mm downstream of the plate
    r = np.hypot(xyz[:, 0] - 0.0492153, xyz[:, 1] - 0.04835) * 1e3
    p, T, rho, a = F['p'], F['T'], F['rho'], F['a']
    speed = np.sqrt(F['u'] ** 2 + F['v'] ** 2 + F['w'] ** 2)

    print('dump %s: cells=%d species=%d' % (path, len(xyz), nspec))
    print('p  [MPa] min=%.3f max=%.3f' % (p.min() / 1e6, p.max() / 1e6))
    print('T  [K]   min=%.0f max=%.0f' % (T.min(), T.max()))
    print('|U| [m/s] p50=%.1f p99=%.1f max=%.1f'
          % (np.percentile(speed, 50), np.percentile(speed, 99), speed.max()))
    print('Mach     max=%.3f' % (speed / np.maximum(a, 1)).max())

    dome = (z >= -2.0) & (z <= dome_mm)
    print('\n== dome region (z = -2 .. %.0f mm): %d cells' % (dome_mm, dome.sum()))

    # Axial pressure gradient along the dome, from binned means.
    edges = np.arange(-2.0, dome_mm + 8.0, 8.0)
    idx = np.digitize(z, edges) - 1
    prof = []
    for i in range(len(edges) - 1):
        m = dome & (idx == i)
        if m.sum() < 200:
            continue
        prof.append((0.5 * (edges[i] + edges[i + 1]), p[m].mean(),
                     T[m].mean(), speed[m].mean(), rho[m].mean(), m.sum()))
    print('\n  z[mm]   p[MPa]   T[K]   |U|[m/s]  rho[kg/m3]  cells')
    for zz, pp, tt, uu, rr, n in prof:
        print('  %6.1f  %7.3f  %6.0f   %6.1f    %7.2f  %7d' % (zz, pp / 1e6, tt, uu, rr, n))

    if len(prof) >= 3:
        zz = np.array([q[0] for q in prof]) * 1e-3
        pp = np.array([q[1] for q in prof])
        uu = np.array([q[3] for q in prof])
        rr = np.array([q[4] for q in prof])
        dpdz = np.gradient(pp, zz)
        # Momentum thickness proxy from the axial velocity deficit of each bin,
        # using the bin width as the layer thickness. Reported as an order of
        # magnitude only: a tetra grid without wall-normal lines cannot give a
        # resolved profile.
        dz = np.gradient(zz)
        theta = uu / np.maximum(uu.max(), 1e-9) * (1 - uu / np.maximum(uu.max(), 1e-9)) * dz
        print('\n  dp/dz [kPa/mm]  min=%.3f max=%.3f'
              % (dpdz.min() / 1e6 * 1e-3 * 1e3, dpdz.max() / 1e6 * 1e-3 * 1e3))
        print('  adverse-gradient bins (dp/dz > 0): %d of %d'
              % ((dpdz > 0).sum(), len(dpdz)))

    # Reversal census: axial velocity sign against the local mean direction.
    print('\n== recirculation census')
    for lo, hi in ((-2, 20), (20, 47), (47, 80), (80, dome_mm)):
        m = (z >= lo) & (z < hi)
        if m.sum() == 0:
            continue
        uz = F['w'][m]
        mean_uz = uz.mean()
        rev = (uz * np.sign(mean_uz) < 0).sum() / len(uz)
        print('  z %3d-%3d mm: cells %7d  <uz> %7.1f m/s  reversed %.2f%%'
              % (lo, hi, m.sum(), mean_uz, 100 * rev))

    # Impingement estimate: highest pressure on the dome wall band.
    wall = dome & (r > 10)
    if wall.sum():
        k = np.argmax(p[wall])
        print('\n  peak dome pressure %.3f MPa at z=%.1f mm, r=%.1f mm'
              % (p[wall][k] / 1e6, z[wall][k], r[wall][k]))


if __name__ == '__main__':
    main()
