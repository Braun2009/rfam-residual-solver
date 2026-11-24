#!/usr/bin/env python3
"""Quantitative cold-flow report from an RFAMTETF1 field dump.

Produces axial (station-averaged) and radial (per-station) profiles of
velocity, Mach, pressure, temperature and species mass fractions, so the
injector near-field can be judged without a viewer.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np

from tetra_field_to_vtu import read_mesh, read_field, verify_phase_pair


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', required=True)
    parser.add_argument('--field', required=True)
    parser.add_argument('--out-prefix', required=True)
    parser.add_argument('--board-z', type=float, default=0.2709527893066406,
                        help='CAD board plane in metres; NOT the injector exit')
    parser.add_argument('--injector-tip-z', type=float, default=0.2858028,
                        help='highest solid injector tip; below this Z includes internal passages')
    parser.add_argument('--axis-x', type=float, default=0.04921534822760748)
    parser.add_argument('--axis-y', type=float, default=0.04835004867381059)
    parser.add_argument('--species', default='O2,C3H8')
    args = parser.parse_args()

    points, cells, volumes = read_mesh(args.mesh)
    centroids, fields, composition, nspecies = read_field(args.field)
    names = [s.strip() for s in args.species.split(',')]
    if len(cells)!=len(fields) or len(names)!=nspecies:
        raise ValueError('mesh, field or species sizes differ')
    centre = np.array([args.axis_x,args.axis_y])
    offset_xy = centroids[:, :2].astype(np.float64) - centre
    radial = np.linalg.norm(offset_xy, axis=1)
    radial_direction = offset_xy/np.maximum(radial[:,None],1e-30)
    radial_velocity = np.sum(fields[:,2:4]*radial_direction,axis=1)
    azimuthal_velocity = -fields[:,2]*radial_direction[:,1]+fields[:,3]*radial_direction[:,0]
    axial_velocity = fields[:,4]
    status_file=Path(args.field+'.status.json')
    status=json.loads(status_file.read_text()) if status_file.exists() else {'converged':False,'production_validated':False,'detail':'missing solver status'}
    Path(args.out_prefix).parent.mkdir(parents=True,exist_ok=True)
    print('solver status:',status)
    if not status.get('converged',False):print('UNCONVERGED_DIAGNOSTIC_FIELD: no ignition or injector-performance conclusion')
    if status.get('boundary_state_valid') is False:
        print('REJECTED_BOUNDARY_STATE: archived field does not satisfy the current boundary model')
    axial = (centroids[:, 2] - args.board_z) * 1e3      # mm from the board
    speed = np.linalg.norm(fields[:, 2:5], axis=1)
    mach = speed / np.maximum(fields[:, 6], 1e-30)

    total_volume = float(volumes.sum())
    mass = float((fields[:, 5].astype(np.float64) * volumes).sum())
    print('mesh cells', len(cells), 'total volume m3', total_volume, 'mass kg', mass)
    print(f'p  [{fields[:,0].min():.4g}, {fields[:,0].max():.4g}] Pa')
    print(f'T  [{fields[:,1].min():.4g}, {fields[:,1].max():.4g}] K')
    print(f'|u|[{speed.min():.4g}, {speed.max():.4g}] m/s  Mach[{mach.min():.4g}, {mach.max():.4g}]')
    print(f'rho[{fields[:,5].min():.4g}, {fields[:,5].max():.4g}] kg/m3')

    # ---- axial station averages (volume weighted) ----
    edges = np.linspace(-2.0, 410.0, 207)
    axial_rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        use = (axial >= lo) & (axial < hi)
        if use.sum() == 0:
            continue
        w = volumes[use]
        wt = w.sum()
        row = {
            'axial_mm': 0.5 * (lo + hi),
            'may_include_injector_passages': bool(args.board_z + lo*1e-3 < args.injector_tip_z),
            'cells': int(use.sum()),
            'area_mm2': wt / ((hi - lo) * 1e-3) * 1e6,
            'p_MPa': float((fields[use, 0] * w).sum() / wt / 1e6),
            'T_K': float((fields[use, 1] * w).sum() / wt),
            'uz_m_s': float((axial_velocity[use] * w).sum() / wt),
            'ur_m_s': float((radial_velocity[use] * w).sum() / wt),
            'utheta_m_s': float((azimuthal_velocity[use] * w).sum() / wt),
            'reverse_volume_fraction':float(w[axial_velocity[use]<0].sum()/wt),
            'slow_lt_5m_s_volume_fraction':float(w[speed[use]<5].sum()/wt),
            'speed_m_s': float((speed[use] * w).sum() / wt),
            'Mach': float((mach[use] * w).sum() / wt),
            'rho': float((fields[use, 5] * w).sum() / wt),
        }
        for index, name in enumerate(names):
            row['Y_' + name] = float((composition[use, index] * w).sum() / wt)
        axial_rows.append(row)
    with open(args.out_prefix + '_axial.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(axial_rows[0].keys()))
        writer.writeheader()
        writer.writerows(axial_rows)
    print('wrote', args.out_prefix + '_axial.csv', len(axial_rows), 'stations')

    # ---- radial profiles at injector-relevant stations ----
    stations = [1.0, 3.0, 6.0, 10.0, 15.0, 25.0, 40.0, 60.0, 100.0]
    radial_rows = []
    for station in stations:
        use = (axial >= station - 0.5) & (axial < station + 0.5)
        if use.sum() == 0:
            continue
        r = radial[use]
        edges_r = np.linspace(0.0, 0.070, 36)
        for lo, hi in zip(edges_r[:-1], edges_r[1:]):
            band = (r >= lo) & (r < hi)
            if band.sum() == 0:
                continue
            w = volumes[use][band]
            wt = w.sum()
            row = {
                'axial_mm': station,
                'may_include_injector_passages': bool(args.board_z+(station-0.5)*1e-3 < args.injector_tip_z),
                'r_mm': 0.5 * (lo + hi) * 1e3,
                'cells': int(band.sum()),
                'p_MPa': float((fields[use, 0][band] * w).sum() / wt / 1e6),
                'T_K': float((fields[use, 1][band] * w).sum() / wt),
                'uz_m_s': float((axial_velocity[use][band] * w).sum() / wt),
                'ur_m_s': float((radial_velocity[use][band] * w).sum() / wt),
                'utheta_m_s': float((azimuthal_velocity[use][band] * w).sum() / wt),
                'speed_m_s': float((speed[use][band] * w).sum() / wt),
                'Mach': float((mach[use][band] * w).sum() / wt),
            }
            for index, name in enumerate(names):
                row['Y_' + name] = float((composition[use, index][band] * w).sum() / wt)
            radial_rows.append(row)
    with open(args.out_prefix + '_radial.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(radial_rows[0].keys()) if radial_rows else ['axial_mm','r_mm','cells'])
        writer.writeheader()
        writer.writerows(radial_rows)
    print('wrote', args.out_prefix + '_radial.csv', len(radial_rows), 'rows')

    summary={**status,'cells':len(cells),'volume_m3':total_volume,'mass_kg':mass,
             'coordinates':'SI, CAD Z downstream','axis_xy_m':centre.tolist(),
             'board_z_m':args.board_z,'injector_tip_z_m':args.injector_tip_z,
             'board_downstream_is_chamber_only':False,
             'p_range_Pa':[float(fields[:,0].min()),float(fields[:,0].max())],
             'T_range_K':[float(fields[:,1].min()),float(fields[:,1].max())],
             'speed_range_m_s':[float(speed.min()),float(speed.max())],
             'Mach_max':float(mach.max()),
             'reverse_volume_fraction':float(volumes[axial_velocity<0].sum()/total_volume),
             'slow_lt_5m_s_volume_fraction':float(volumes[speed<5].sum()/total_volume),
             'species':{name:{'min':float(composition[:,i].min()),'max':float(composition[:,i].max()),
                              'inventory_kg':float((volumes*fields[:,5]*composition[:,i]).sum())} for i,name in enumerate(names)}}
    boundary_path=Path(args.field+'.boundary.csv')
    if boundary_path.exists():
        budgets=list(csv.DictReader(boundary_path.open(newline='')))
        inflow=sum(float(row['target_inflow_kg_s']) for row in budgets)
        net=sum(float(row['outward_mass_kg_s']) for row in budgets)
        summary['boundary_budget']={'inlet_target_kg_s':inflow,
            'net_outward_mass_kg_s':net,'relative_mass_imbalance':abs(net)/inflow if inflow>0 else None,
            'net_outward_species_kg_s':{name:sum(float(row['outward_'+name+'_kg_s']) for row in budgets) for name in names},
            'net_outward_sensible_total_energy_W':sum(float(row['outward_energy_W']) for row in budgets)}
    phase_path=Path(args.field+'.phase.csv')
    if phase_path.exists():
        verify_phase_pair(args.field,phase_path)
        phase=np.genfromtxt(phase_path,delimiter=',',names=True)
        if len(phase)!=len(cells) or not np.array_equal(phase['cell'],np.arange(len(cells))):
            raise ValueError('phase sidecar and mesh cell order differ')
        two_phase=phase['region']==5
        summary['two_phase_cells']=int(two_phase.sum())
        summary['two_phase_volume_fraction']=float(volumes[two_phase].sum()/total_volume)
        summary['projected_cells']=int(np.count_nonzero(phase['projected']))
        summary['two_phase_vapor_mass_fraction_range']=([float(phase['vapor_mass_fraction'][two_phase].min()),
            float(phase['vapor_mass_fraction'][two_phase].max())] if two_phase.any() else None)
    Path(args.out_prefix+'_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False),encoding='utf-8')
    print('wrote',args.out_prefix+'_summary.json')


if __name__ == '__main__':
    main()
