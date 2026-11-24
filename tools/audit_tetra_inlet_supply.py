#!/usr/bin/env python3
"""Audit pure-propellant mass/h0 inlet compatibility using CoolProp 7.2.0.

Supply p/T are interpreted as stagnation reservoir conditions, as in the
current tetra driver. This is a model assumption, not a reinterpretation of
CAD/design inputs. A negative entropy change rejects this boundary model at
the saved field; it does not establish a defect in the injector design.
"""
import argparse
import hashlib
import json
import math
import struct
from pathlib import Path

import CoolProp
from CoolProp import AbstractState, CoolProp as CP
import numpy as np
from tetra_field_to_vtu import read_mesh, read_field


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', required=True)
    parser.add_argument('--field', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--lox-p', type=float, default=6e6)
    parser.add_argument('--lox-T', type=float, default=133.5)
    parser.add_argument('--lox-mdot', type=float, default=2.7021)
    parser.add_argument('--fuel-p', type=float, default=4.5e6)
    parser.add_argument('--fuel-T', type=float, default=453.5)
    parser.add_argument('--fuel-mdot', type=float, default=1.0134)
    args = parser.parse_args()
    if CoolProp.__version__ != '7.2.0':
        raise RuntimeError('Reference audit requires CoolProp 7.2.0')
    mesh_path, field_path = Path(args.mesh), Path(args.field)
    points, cells, _ = read_mesh(mesh_path)
    _, fields, _, nspecies = read_field(field_path)
    raw = mesh_path.read_bytes()
    np_, nc, nf = struct.unpack_from('<QQQ', raw, 8)
    face_offset = 32 + np_*24 + nc*24
    if len(raw) != face_offset + nf*24 or len(fields) != len(cells) or nspecies != 2:
        raise ValueError('Mesh/field size mismatch')
    faces = np.frombuffer(raw, '<i4', nf*6, face_offset).reshape(-1, 6)
    triangles = points[faces[:, :3]]
    area = .5*np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0],
                                      triangles[:, 2]-triangles[:, 0]), axis=1)
    rows = []
    for family, fluid, ids, pressure, temperature, mdot in [
        ('lox', 'Oxygen', [10, 11, 12], args.lox_p, args.lox_T, args.lox_mdot),
        ('fuel', 'Propane', [20, 21, 22], args.fuel_p, args.fuel_T, args.fuel_mdot),
    ]:
        if not all(math.isfinite(v) and v > 0 for v in [pressure, temperature, mdot]):
            raise ValueError('Supply p, T, and mdot must be finite and positive')
        state = AbstractState('HEOS', fluid)
        state.update(CP.PT_INPUTS, pressure, temperature)
        h0, s0 = state.hmass(), state.smass()
        for zone in ids:
            selected = (faces[:, 4] < 0) & (faces[:, 5] == zone)
            if not selected.any():
                raise ValueError(f'Missing inlet zone {zone}')
            patch_area = float(area[selected].sum())
            backpressure = float(np.average(fields[faces[selected, 3], 0], weights=area[selected]))
            mass_flux = (mdot/len(ids))/patch_area
            row = dict(zone=zone, family=family, area_m2=patch_area,
                       target_mdot_kg_s=mdot/len(ids), static_p_Pa=backpressure,
                       supply_p_Pa=pressure, supply_T_K=temperature)
            try:
                trial_t = temperature
                closed = False
                for _ in range(25):
                    state.update(CP.PT_INPUTS, backpressure, trial_t)
                    if state.phase() == CP.iphase_twophase:
                        raise ValueError('Pure inlet PT inversion cannot traverse a two-phase interval')
                    velocity = mass_flux/state.rhomass()
                    residual = state.hmass()+.5*velocity**2-h0
                    if abs(residual) <= 1e-8*max(abs(h0), 1e5):
                        closed = True
                        break
                    dpdt = state.first_partial_deriv(CP.iP, CP.iT, CP.iDmass)
                    dpdrho = state.first_partial_deriv(CP.iP, CP.iDmass, CP.iT)
                    derivative = state.cpmass()+velocity**2*dpdt/(state.rhomass()*dpdrho)
                    if not math.isfinite(derivative) or derivative <= 0:
                        raise ValueError('Invalid h0 residual derivative')
                    trial_t += float(np.clip(-residual/derivative, -.05*trial_t, .05*trial_t))
                if not closed:
                    raise ValueError('h0 inversion did not converge')
                ds = state.smass()-s0
                compatible = backpressure < pressure and ds >= -1e-7*max(abs(s0), 1.)
                row.update(static_T_K=state.T(), rho_kg_m3=state.rhomass(),
                           velocity_m_s=velocity, h0_error_J_kg=residual,
                           supply_entropy_J_kgK=s0, entropy_change_J_kgK=ds,
                           compatible=compatible)
            except Exception as exc:
                row.update(compatible=False, error=str(exc))
            rows.append(row)
    report = {
        'coolprop_version': CoolProp.__version__, 'coolprop_git': CoolProp.__gitrevision__,
        'assumption': 'adiabatic no-work reservoir; supplied p/T interpreted as stagnation conditions',
        'scope': 'six main pure-propellant inlet zones; torch is not audited',
        'mesh_sha256': hashlib.sha256(raw).hexdigest(),
        'field_sha256': hashlib.sha256(field_path.read_bytes()).hexdigest(),
        'all_main_inlets_compatible': all(row['compatible'] for row in rows), 'rows': rows,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('INLET_SUPPLY_AUDIT', 'COMPATIBLE' if report['all_main_inlets_compatible'] else 'INCOMPATIBLE')
    for row in rows:
        print('zone', row['zone'], 'delta_s_J_kgK', row.get('entropy_change_J_kgK'),
              'compatible', row['compatible'])
    print(out)


if __name__ == '__main__':
    main()
