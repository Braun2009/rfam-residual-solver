#!/usr/bin/env python3
"""CoolProp operating envelope for fixed design propellant mass flow.

This is a model-only boundary audit. It does not change CAD, target mdot or
solver settings. Supply p/T are treated as stagnation reservoir conditions and
capacity is the ideal isentropic upper bound through the audited zone area.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import struct
from pathlib import Path

from CoolProp import AbstractState, __version__, __gitrevision__
from CoolProp import CoolProp as CP
import numpy as np


def read_mesh(path: Path):
    raw = path.read_bytes()
    if raw[:8] != b'RFAMTET1':
        raise ValueError('not an RFAMTET1 mesh')
    np_, nc, nf = struct.unpack_from('<QQQ', raw, 8)
    offset = 32
    points = np.frombuffer(raw, '<f8', np_ * 3, offset).reshape(-1, 3)
    offset += np_ * 24
    cells = np.frombuffer(raw, '<i4', nc * 4, offset).reshape(-1, 4)
    offset += nc * 16
    offset += nc * 8
    faces = np.frombuffer(raw, '<i4', nf * 6, offset).reshape(-1, 6)
    tri = points[faces[:, :3]]
    areas = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) * .5
    return raw, faces, areas


def mass_flux(state, p0, t0, p_back, area):
    state.update(CP.PT_INPUTS, p0, t0)
    h0, s0 = state.hmass(), state.smass()
    if p_back >= p0:
        return 0.0
    state.update(CP.PSmass_INPUTS, p_back, s0)
    drop = h0 - state.hmass()
    if drop <= 0.0:
        return 0.0
    return state.rhomass() * math.sqrt(2.0 * drop) * area


def solve_max_backpressure(state, p0, t0, area, target):
    lo = 0.5 * p0
    hi = p0 * (1.0 - 1.0e-12)
    if mass_flux(state, p0, t0, lo, area) < target:
        return None
    # The low-pressure end passes the target and the near-reservoir end does not.
    for _ in range(100):
        mid = .5 * (lo + hi)
        if mass_flux(state, p0, t0, mid, area) >= target:
            lo = mid
        else:
            hi = mid
    return .5 * (lo + hi)


def solve_supply_pressure(state, t0, p_back, area, target, lo, hi):
    def f(p0):
        return mass_flux(state, p0, t0, p_back, area)
    if f(hi) < target:
        return None
    for _ in range(100):
        mid = .5 * (lo + hi)
        if f(mid) >= target:
            hi = mid
        else:
            lo = mid
    return .5 * (lo + hi)


def solve_supply_temperature(state, p0, p_back, area, target, lo, hi):
    def f(t0):
        return mass_flux(state, p0, t0, p_back, area)
    # For this cryogenic propane range, lower reservoir T increases ideal
    # capacity. Keep the monotonic bracket explicit and fail closed otherwise.
    if f(lo) < target or f(hi) > target:
        return None
    for _ in range(100):
        mid = .5 * (lo + hi)
        if f(mid) >= target:
            lo = mid
        else:
            hi = mid
    return .5 * (lo + hi)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--fuel-p', type=float, default=4.5e6)
    ap.add_argument('--fuel-T', type=float, default=453.5)
    ap.add_argument('--fuel-mdot', type=float, default=1.05)
    ap.add_argument('--field-fuel-p', type=float, default=None,
                    help='diagnostic owner/zone backpressure, not a solver input')
    ap.add_argument('--min-supply-p', type=float, default=1.0e6)
    ap.add_argument('--max-supply-p', type=float, default=20.0e6)
    ap.add_argument('--min-supply-T', type=float, default=300.0)
    args = ap.parse_args()
    if not all(math.isfinite(v) and v > 0.0 for v in
               [args.fuel_p, args.fuel_T, args.fuel_mdot, args.min_supply_p,
                args.max_supply_p, args.min_supply_T]):
        ap.error('all thermodynamic inputs must be positive and finite')
    raw, faces, areas = read_mesh(Path(args.mesh))
    rows = []
    state = AbstractState('HEOS', 'Propane')
    target = args.fuel_mdot / 3.0
    for zone in [20, 21, 22]:
        selected = (faces[:, 4] < 0) & (faces[:, 5] == zone)
        if not selected.any():
            raise ValueError(f'missing fuel zone {zone}')
        area = float(areas[selected].sum())
        face = args.field_fuel_p
        max_back = solve_max_backpressure(state, args.fuel_p, args.fuel_T, area, target)
        if face is None:
            face = max_back if max_back is not None else 4.0e6
        cap_at_face = mass_flux(state, args.fuel_p, args.fuel_T, face, area)
        required_p = solve_supply_pressure(state, args.fuel_T, face, area, target,
                                           max(face * (1.0 + 1e-10), args.min_supply_p), args.max_supply_p)
        required_t = solve_supply_temperature(state, args.fuel_p, face, area, target,
                                              args.min_supply_T, args.fuel_T)
        rows.append({
            'zone': zone,
            'area_m2': area,
            'target_mdot_kg_s': target,
            'supply_p_Pa': args.fuel_p,
            'supply_T_K': args.fuel_T,
            'diagnostic_face_p_Pa': face,
            'ideal_capacity_at_supply_kg_s': cap_at_face,
            'capacity_ratio': cap_at_face / target,
            'max_consistent_face_pressure_Pa': max_back,
            'required_supply_pressure_at_face_Pa': required_p,
            'required_supply_temperature_at_face_K': required_t,
            'reachable_under_current_supply': cap_at_face >= target,
        })
    report = {
        'scope': 'boundary operating-point diagnostic; no solver/CAD/design modification',
        'assumption': 'adiabatic no-work stagnation reservoir, ideal isentropic upper bound',
        'coolprop_version': __version__,
        'coolprop_git': __gitrevision__,
        'mesh_sha256': hashlib.sha256(raw).hexdigest(),
        'mesh': str(args.mesh),
        'design_total_fuel_mdot_kg_s': args.fuel_mdot,
        'rows': rows,
        'all_current_points_reachable': all(r['reachable_under_current_supply'] for r in rows),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print('FUEL_OPERATING_ENVELOPE',
          'REACHABLE' if report['all_current_points_reachable'] else 'NOT_REACHABLE')
    for row in rows:
        print('zone', row['zone'],
              'face_MPa=%.6f' % (row['diagnostic_face_p_Pa'] / 1e6),
              'capacity_ratio=%.6f' % row['capacity_ratio'],
              'max_face_MPa=%.6f' % (row['max_consistent_face_pressure_Pa'] / 1e6
                                     if row['max_consistent_face_pressure_Pa'] else float('nan')),
              'required_supply_MPa=%.6f' % (row['required_supply_pressure_at_face_Pa'] / 1e6
                                             if row['required_supply_pressure_at_face_Pa'] else float('nan')))
    print(out)


if __name__ == '__main__':
    main()
