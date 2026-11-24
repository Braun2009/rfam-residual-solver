#!/usr/bin/env python3
"""Isentropic flow capacity of the main inlets under an adiabatic reservoir.

Model-only audit. The driver interprets the configured supply p/T as the
stagnation state of an adiabatic no-work reservoir and prescribes the design
mass flow through the CAD inlet area. For each inlet zone this script computes
the isentropic (maximum) mass flow that the same area can pass when the
discharge static pressure equals the area-weighted average owner pressure of
the saved field, then reports the requested/available ratio. It also inverts
the requirement in three ways (discharge pressure, supply pressure, supply
temperature). Nothing here changes CAD, throat areas, or design flow rates.
"""
import argparse
import hashlib
import json
import math
import struct
from pathlib import Path

import numpy as np
from CoolProp import AbstractState
from CoolProp import CoolProp as CP

from tetra_field_to_vtu import read_mesh, read_field

FLUIDS = {'lox': 'Oxygen', 'fuel': 'Propane'}
SPEED_OF_SOUND_EPS = 1e-9


def isentropic_flow(state, p0, t0, backpressure, area):
    """Return (mdot, velocity, rho, T, s_change, choked_mdot, p_crit)."""
    state.update(CP.PT_INPUTS, p0, t0)
    h0, s0 = state.hmass(), state.smass()
    choked_mdot = None
    p_crit = None

    def flow_at(p):
        try:
            state.update(CP.PSmass_INPUTS, p, s0)
        except Exception:
            return None
        h, rho = state.hmass(), state.rhomass()
        kinetic = 2.0*max(h0 - h, 0.0)
        if kinetic <= 0.0:
            return None
        u = math.sqrt(kinetic)
        return rho*area*u, u, rho, state.T(), state.smass() - s0

    # Locate the sonic point by maximizing the isentropic mass flux.
    lo, hi = 1.0, p0
    for _ in range(200):
        m1, m2 = lo + (hi-lo)/3, hi - (hi-lo)/3
        f1 = flow_at(m1)
        f2 = flow_at(m2)
        if f1 is None or f2 is None:
            break
        if f1[0] < f2[0]:
            lo = m1
        else:
            hi = m2
    p_sonic = 0.5*(lo+hi)
    sonic = flow_at(p_sonic)
    if sonic is not None:
        choked_mdot = sonic[0]
        p_crit = p_sonic
    if backpressure is not None and backpressure >= p_crit:
        result = flow_at(backpressure)
        return (*result[:5], choked_mdot, p_crit)
    return (choked_mdot, sonic[1], sonic[2], sonic[3], sonic[4],
            choked_mdot, p_crit)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', required=True)
    parser.add_argument('--field', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    mesh_path, field_path = Path(args.mesh), Path(args.field)
    points, cells, _ = read_mesh(mesh_path)
    _, fields, _, nspecies = read_field(field_path)
    raw = mesh_path.read_bytes()
    np_, nc, nf = struct.unpack_from('<QQQ', raw, 8)
    face_offset = 32 + np_*24 + nc*24
    if len(raw) != face_offset + nf*24 or len(fields) != len(cells):
        raise ValueError('mesh/field size mismatch')
    faces = np.frombuffer(raw, '<i4', nf*6, face_offset).reshape(-1, 6)
    triangles = points[faces[:, :3]]
    area = .5*np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0],
                                      triangles[:, 2]-triangles[:, 0]), axis=1)

    zone_spec = [('lox', [10, 11, 12]), ('fuel', [20, 21, 22])]

    rows = []
    for family, ids in zone_spec:
        fluid = FLUIDS[family]
        state = AbstractState('HEOS', fluid)
        for zone in ids:
            selected = (faces[:, 4] < 0) & (faces[:, 5] == zone)
            if not selected.any():
                raise ValueError(f'missing inlet zone {zone}')
            patch_area = float(area[selected].sum())
            p_face = float(np.average(fields[faces[selected, 3], 0],
                                      weights=area[selected]))
            rows.append(dict(family=family, fluid=fluid, zone=zone,
                             area_m2=patch_area, p_face_Pa=p_face, state=state))
    report = {'coolprop_version': __import__('CoolProp').__version__,
              'assumption': ('adiabatic no-work reservoir; configured supply p/T '
                             'interpreted as stagnation state; isentropic upper '
                             'bound on realizable mass flow'),
              'scope': 'mass-flow consistency of the inlet boundary model only',
              'mesh_sha256': hashlib.sha256(raw).hexdigest(),
              'field_sha256': hashlib.sha256(field_path.read_bytes()).hexdigest(),
              'rows': []}
    for row in rows:
        row.pop('state', None)
        report['rows'].append(row)
    # Fill in per-family supply conditions from the driver defaults.
    for row in report['rows']:
        row['supply_p_Pa'] = 6e6 if row['family'] == 'lox' else 4.5e6
        row['supply_T_K'] = 133.5 if row['family'] == 'lox' else 453.5
    # 2026-09-17 design table; the Cd model injects Cd*capacity, so the
    # ratio below compares the *design* flow with the isentropic capacity.
    argv_mdot = {'lox': 2.8/3, 'fuel': 1.05/3}
    for row in report['rows']:
        row['target_mdot_kg_s'] = argv_mdot[row['family']]
        row['Cd'] = 0.75
        state = AbstractState('HEOS', row['fluid'])
        mdot, u, rho, t, ds, m_choked, p_crit = isentropic_flow(
            state, row['supply_p_Pa'], row['supply_T_K'], row['p_face_Pa'],
            row['area_m2'])
        row.update(isentropic_mdot_kg_s=mdot, isentropic_velocity_m_s=u,
                   isentropic_rho_kg_m3=rho, isentropic_T_K=t,
                   isentropic_ds_J_kgK=ds, choked_mdot_kg_s=m_choked,
                   critical_pressure_Pa=p_crit,
                   ratio_target_to_isentropic=row['target_mdot_kg_s']/mdot)
        row['isentropic_curve'] = [
            {'p_back_Pa': p, 'mdot_kg_s': isentropic_flow(
                state, row['supply_p_Pa'], row['supply_T_K'], p,
                row['area_m2'])[0]}
            for p in sorted({row['p_face_Pa'], 4.0e6, 3.5e6, 3.0e6,
                             row['critical_pressure_Pa']}, reverse=True)]
        if row['ratio_target_to_isentropic'] > 1.0:
            # pressure that would just pass the target flow isentropically
            lo, hi = 1.0, row['supply_p_Pa']
            for _ in range(80):
                mid = 0.5*(lo+hi)
                f = isentropic_flow(state, row['supply_p_Pa'], row['supply_T_K'],
                                    mid, row['area_m2'])[0]
                if f is None:
                    break
                if f >= row['target_mdot_kg_s']:
                    lo = mid
                else:
                    hi = mid
            row['max_consistent_face_pressure_Pa'] = 0.5*(lo+hi)
            # supply pressure that would just pass the target flow at p_face
            lo, hi = row['supply_p_Pa'], 20e6
            for _ in range(80):
                mid = 0.5*(lo+hi)
                f = isentropic_flow(state, mid, row['supply_T_K'],
                                    row['p_face_Pa'], row['area_m2'])[0]
                if f is None:
                    break
                if f < row['target_mdot_kg_s']:
                    lo = mid
                else:
                    hi = mid
            row['min_consistent_supply_pressure_Pa'] = 0.5*(lo+hi)
            # supply temperature that would just pass the target flow
            lo, hi = max(0.6*row['supply_T_K'], 80.0), row['supply_T_K']
            for _ in range(80):
                mid = 0.5*(lo+hi)
                f = isentropic_flow(state, row['supply_p_Pa'], mid,
                                    row['p_face_Pa'], row['area_m2'])[0]
                if f is None:
                    break
                if f < row['target_mdot_kg_s']:
                    hi = mid
                else:
                    lo = mid
            row['max_consistent_supply_temperature_K'] = 0.5*(lo+hi)
    report['all_targets_isentropically_reachable'] = all(
        row['ratio_target_to_isentropic'] <= 1.0 for row in report['rows'])
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n',
                   encoding='utf-8')
    print('INLET_CAPACITY', 'REACHABLE' if
          report['all_targets_isentropically_reachable'] else 'NOT_REACHABLE')
    for row in report['rows']:
        print('zone', row['zone'], row['family'],
              'ratio=%.4f' % row['ratio_target_to_isentropic'],
              'p_face=%.0f' % row['p_face_Pa'],
              'choked=%.4f' % row['choked_mdot_kg_s'],
              'p_crit=%.0f' % row['critical_pressure_Pa'],
              ('max_face_p=%.0f' % row['max_consistent_face_pressure_Pa'])
              if 'max_consistent_face_pressure_Pa' in row else '')
    print(out)


if __name__ == '__main__':
    main()
