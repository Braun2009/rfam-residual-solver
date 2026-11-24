#!/usr/bin/env python3
"""Cold-flow field characteristics for the A-route (real injector) dumps.

Produces the consolidated characteristic set asked for by the injector
verification: what the flow field actually looks like, with every number
traceable to the dump and the mesh.

Sections
  * zones: injector block vs chamber (cells, volume, mass).
  * streams: liquid-LOX and pure-fuel inventory, the front positions (how far
    each propellant has advanced), and how much of each has entered the chamber.
  * mixing: O2 mass-fraction moments and the unmixedness theta = <(1-Y)^2>-style
    segregation measure in the chamber, plus the fuel-rich / O2-rich fractions.
  * recirculation: backflow (w<0) and slow-zone volume fractions in the chamber.
  * pattern: azimuthal m=1/m=3 modes of the face-plane annulus and their phases
    (the three-pintle signature and the swirl rotation).
  * central zone: volume-integrated kinetic second moments (units J), not a
    momentum-flux tensor or proof of vector momentum cancellation.
  * equilibrium phase: actual EOS phase-sidecar region==5, when present.
    Cryogenic/dense thresholds are reported separately and are NOT phase.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tetra_field_to_vtu import read_mesh, read_field, verify_phase_pair  # noqa: E402


def fit_azimuthal_modes(theta: np.ndarray, weight: np.ndarray, field: np.ndarray) -> dict:
    """Weighted angular fit with an intercept, so a constant field has zero modes."""
    if len(theta) != len(weight) or len(theta) != len(field) or len(theta) < 5:
        raise ValueError('azimuthal fit needs at least five matching samples')
    w = np.asarray(weight, dtype=np.float64)
    values = np.asarray(field, dtype=np.float64)
    if not (np.isfinite(w).all() and (w >= 0).all() and w.sum() > 0
            and np.isfinite(values).all() and np.isfinite(theta).all()):
        raise ValueError('invalid azimuthal weights or field')
    design = np.column_stack((np.ones(len(theta)), np.cos(theta), np.sin(theta),
                              np.cos(3*theta), np.sin(3*theta)))
    root = np.sqrt(w / w.sum())
    weighted = design * root[:, None]
    if np.linalg.cond(weighted.T @ weighted) > 1e8:
        raise ValueError('azimuthal mesh sampling is rank deficient')
    mean, c1, s1, c3, s3 = np.linalg.lstsq(weighted, values*root, rcond=None)[0]
    m1, m3 = float(np.hypot(c1, s1)), float(np.hypot(c3, s3))
    return dict(mean=float(mean), m1=m1, m3=m3,
                m3_over_mean=m3/abs(mean) if mean else 0.,
                m3_over_m1=m3/max(m1,1e-30),
                phase_deg=float(np.degrees(np.arctan2(s3,c3))/3.0 % 120.0),
                method='weighted_least_squares_intercept_m1_m3')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh', required=True)
    ap.add_argument('--field', required=True)
    ap.add_argument('--board-z', type=float, default=0.2709527893066406)
    ap.add_argument('--injector-tip-z', type=float, default=0.2858028,
                    help='maximum injector solid tip Z from audited CAD (metres)')
    ap.add_argument('--axis-x', type=float, default=0.04921534822760748)
    ap.add_argument('--axis-y', type=float, default=0.04835004867381059)
    ap.add_argument('--label', default='')
    ap.add_argument('--json', default=None)
    ap.add_argument('--phase', default=None, help='EOS phase CSV (defaults to FIELD.phase.csv if present)')
    args = ap.parse_args()

    _, cells, vol = read_mesh(args.mesh)
    cent, f, y, _ = read_field(args.field)
    p, T, u, v, w, rho = (f[:, i] for i in range(6))
    spd = np.sqrt(u*u + v*v + w*w)
    z = cent[:, 2]
    r = np.hypot(cent[:, 0] - args.axis_x, cent[:, 1] - args.axis_y)
    inner = z < args.board_z
    downstream = ~inner  # Includes fluid INSIDE protruding injector passages.
    post_tip = z >= args.injector_tip_z
    volD = float(vol[downstream].sum())
    volTip = float(vol[post_tip].sum())
    Y = y[:, 0]

    def frac(mask, base=downstream):
        return 100.0*float(vol[mask & base].sum())/max(float(vol[base].sum()),1e-30)

    out = dict(label=args.label, cells=int(len(cells)),
               injector_cells=int(inner.sum()),board_downstream_cells=int(downstream.sum()),
               post_tip_cells=int(post_tip.sum()),
               injector_volume_cm3=float(vol[inner].sum()*1e6),
               board_downstream_volume_cm3=volD*1e6,
               post_tip_volume_cm3=volTip*1e6,
               board_z_m=args.board_z,injector_tip_z_m=args.injector_tip_z,
               p_MPa=[float(p.min()/1e6), float(p.max()/1e6)],
               T_K=[float(T.min()), float(T.max())],
               speed_ms=[float(spd.min()), float(spd.max())],
               mach_max=float(np.max(spd/f[:, 6])),
               rho_kg_m3=[float(rho.min()), float(rho.max())])
    # A density threshold is only a dense-fluid proxy, not a liquid phase tag.
    liq = rho > 500.
    o2pure = Y > 0.95
    fuel = Y < 0.2
    out['streams'] = dict(
        density_proxy_definition='rho>500 kg/m3 (not an EOS liquid classification)',
        dense_cells=int(liq.sum()), dense_mass_mg=float(np.sum(rho[liq]*vol[liq])*1e6),
        dense_in_board_downstream_cells=int((liq & downstream).sum()),
        dense_post_tip_cells=int((liq & post_tip).sum()),
        pure_o2_in_board_downstream_cells=int((o2pure & downstream).sum()),
        pure_o2_post_tip_cells=int((o2pure & post_tip).sum()),
        fuel_in_board_downstream_cells=int((fuel & downstream).sum()),
        fuel_post_tip_cells=int((fuel & post_tip).sum()),
        o2_max_z=float(z[o2pure].max()) if o2pure.any() else None,
        fuel_max_z=float(z[fuel].max()) if fuel.any() else None,
        o2_post_tip_max_z=float(z[o2pure & post_tip].max()) if (o2pure & post_tip).any() else None,
        fuel_post_tip_max_z=float(z[fuel & post_tip].max()) if (fuel & post_tip).any() else None)
    # Board-downstream includes injector internals; neither mask implies steady mixing.
    def mixing(mask):
        if not mask.any():
            return dict(cells=0,Y_mean=None,Y_std=None,Y_min=None,Y_max=None,
                        fuel_rich_frac=None,o2_rich_frac=None,nonbackground_volume_frac=None)
        total=float(vol[mask].sum())
        ym=float(np.sum(Y[mask]*vol[mask])/total)
        return dict(cells=int(mask.sum()),Y_mean=ym,
                    Y_std=float(np.sqrt(np.sum((Y[mask]-ym)**2*vol[mask])/total)),
                    Y_min=float(Y[mask].min()),Y_max=float(Y[mask].max()),
                    fuel_rich_frac=frac(Y<0.4,mask),o2_rich_frac=frac(Y>0.95,mask),
                    nonbackground_volume_frac=frac(np.abs(Y-0.784)>0.1,mask))
    out['board_downstream_mixing']=mixing(downstream)
    out['post_tip_mixing']=mixing(post_tip)
    # Recirculation/slow zones are diagnostic until the field is stationary.
    def recirculation(mask):
        if not mask.any():
            return dict(cells=0,backflow_volume_frac=None,slow_lt_5_frac=None,
                        slow_lt_20_frac=None,mach_gt_0p3_frac=None)
        return dict(cells=int(mask.sum()),backflow_volume_frac=frac(w<0.,mask),
                    slow_lt_5_frac=frac(spd<5.,mask),slow_lt_20_frac=frac(spd<20.,mask),
                    mach_gt_0p3_frac=frac(spd/np.maximum(f[:,6],1e-30)>0.3,mask))
    out['board_downstream_recirculation']=recirculation(downstream)
    out['post_tip_recirculation']=recirculation(post_tip)

    # Thresholds identify cold/dense fluid, not an EOS two-phase region.
    dense = (rho > 500.) | ((T < 300.) & (Y > 0.7))
    out['cryogenic_or_dense'] = dict(cells=int(dense.sum()),
        board_downstream_volume_percent=frac(dense),
        board_downstream_cells=int((dense & downstream).sum()),
        post_tip_cells=int((dense & post_tip).sum()),
        definition='rho>500 kg/m3 or (T<300 K and Y_O2>0.7)')
    phase_path = Path(args.phase) if args.phase else Path(args.field + '.phase.csv')
    if phase_path.exists():
        verify_phase_pair(args.field,phase_path)
        phase = np.genfromtxt(phase_path,delimiter=',',skip_header=1,usecols=(0,1),dtype=np.int32)
        if len(phase)!=len(cells) or not np.array_equal(phase[:,0],np.arange(len(cells))):
            raise ValueError('phase sidecar/mesh cell order differs')
        region=phase[:,1]
        tp=region==5
        liquid=region==2
        out['eos_liquid'] = dict(source=str(phase_path),definition='EOS region==2',
            cells=int(liquid.sum()),board_downstream_cells=int((liquid & downstream).sum()),
            post_tip_cells=int((liquid & post_tip).sum()))
        out['two_phase'] = dict(source=str(phase_path),definition='EOS region==5',
            cells=int(tp.sum()),board_downstream_cells=int((tp & downstream).sum()),
            post_tip_cells=int((tp & post_tip).sum()),
            total_volume_percent=100.*float(vol[tp].sum())/float(vol.sum()),
            board_downstream_volume_percent=frac(tp),
            post_tip_volume_percent=frac(tp,post_tip))
    else:
        out['two_phase'] = dict(source=None,definition='EOS sidecar unavailable; not inferred',
            cells=None,board_downstream_cells=None,post_tip_cells=None,
            total_volume_percent=None,board_downstream_volume_percent=None,
            post_tip_volume_percent=None)

    # This annulus includes internal passages above the plate; not a chamber-only pattern.
    ann = downstream & (z < args.board_z + 20e-3) & (r > 8e-3) & (r < 18e-3)
    if ann.sum() > 100:
        th = np.arctan2(cent[ann, 1] - args.axis_y, cent[ann, 0] - args.axis_x)
        wgt = rho[ann]*vol[ann]
        modes = {}
        for name, field in (('Y_O2', Y[ann]), ('u_z', w[ann]), ('T', T[ann]), ('p', p[ann])):
            modes[name] = fit_azimuthal_modes(th,wgt,field)
        out['board_annulus_pattern'] = modes
    # This volume also includes torch passage fluid until the top of the tip.
    zone = downstream & (r < 5e-3) & (z < args.board_z + 20e-3)
    if zone.sum() > 100:
        V = vol[zone]
        ur = (u[zone]*(cent[zone, 0]-args.axis_x) + v[zone]*(cent[zone, 1]-args.axis_y))/np.maximum(r[zone], 1e-30)
        ut = (-u[zone]*(cent[zone, 1]-args.axis_y) + v[zone]*(cent[zone, 0]-args.axis_x))/np.maximum(r[zone], 1e-30)
        Mzz = float(np.sum(rho[zone]*w[zone]*w[zone]*V))
        Mrr = float(np.sum(rho[zone]*ur*ur*V))
        Mtt = float(np.sum(rho[zone]*ut*ut*V))
        out['board_core_motion'] = dict(cells=int(zone.sum()),
                              Mzz=Mzz, Mrr=Mrr, Mtt=Mtt,
                              moment_units='J = integral rho*u_i^2 dV',
                              radial_to_axial_motion_ratio=Mrr/max(Mzz, 1e-30),
                              momentum_cancellation_verified=False,
                              speed_mean=float(spd[zone].mean()), speed_std=float(spd[zone].std()),
                              radial_speed_rms=float(np.sqrt(np.mean(ur**2))),
                              Y_mean=float(Y[zone].mean()), Y_std=float(Y[zone].std()),
                              Y_min=float(Y[zone].min()),
                              T_std=float(T[zone].std()), p_std_MPa=float(p[zone].std()/1e6))
    print(f"=== cold-flow characteristics {args.label} ===")
    print(f"zones: upstream of board {out['injector_cells']} cells ({out['injector_volume_cm3']:.4g} cm3), "
          f"board-downstream INCLUDING passages {out['board_downstream_cells']} cells "
          f"({out['board_downstream_volume_cm3']:.4g} cm3), "
          f"certain post-tip {out['post_tip_cells']} cells ({out['post_tip_volume_cm3']:.4g} cm3)")
    print(f"state: p {out['p_MPa'][0]:.3f}-{out['p_MPa'][1]:.3f} MPa, "
          f"T {out['T_K'][0]:.1f}-{out['T_K'][1]:.1f} K, |u| {out['speed_ms'][0]:.2f}-{out['speed_ms'][1]:.1f} m/s, "
          f"Mach<= {out['mach_max']:.3f}, rho {out['rho_kg_m3'][0]:.1f}-{out['rho_kg_m3'][1]:.1f}")
    s = out['streams']
    print(f"density proxy (rho>500, NOT liquid phase): {s['dense_cells']} cells / {s['dense_mass_mg']:.3f} mg "
          f"(board-downstream {s['dense_in_board_downstream_cells']}, post-tip {s['dense_post_tip_cells']}), "
          f"pure O2 board-downstream/post-tip {s['pure_o2_in_board_downstream_cells']}/{s['pure_o2_post_tip_cells']}, "
          f"fuel-rich {s['fuel_in_board_downstream_cells']}/{s['fuel_post_tip_cells']}")
    def format_front(value):
        return f'{value:.6g}' if value is not None else 'not present'
    print(f"         maximum pure-stream z: O2={format_front(s['o2_max_z'])}, "
          f"fuel={format_front(s['fuel_max_z'])}; tip={args.injector_tip_z:.6g}, "
          f"O2 post-tip={format_front(s['o2_post_tip_max_z'])}")
    for label in ('board_downstream','post_tip'):
        m = out[label+'_mixing']
        rc = out[label+'_recirculation']
        if m['cells']==0:
            print(f'{label}: empty region')
        else:
            print(f"{label} (UNCONVERGED): Y_O2 {m['Y_min']:.3f}-{m['Y_max']:.3f} "
                  f"(mean {m['Y_mean']:.4f}, std {m['Y_std']:.4f}); "
                  f"backflow {rc['backflow_volume_frac']:.3f}%, "
                  f"<20 m/s {rc['slow_lt_20_frac']:.3f}%")
    phase = out['two_phase']
    if phase['source']:
        print(f"EOS two-phase: {phase['cells']} cells, {phase['total_volume_percent']:.6f}% total volume "
              f"(board-downstream {phase['board_downstream_cells']} cells / "
              f"{phase['board_downstream_volume_percent']:.6f}% volume; "
              f"post-tip {phase['post_tip_cells']} cells)")
    else:
        print('EOS two-phase: unavailable (phase sidecar missing)')
    dense = out['cryogenic_or_dense']
    print(f"cryogenic/dense threshold (NOT phase): {dense['cells']} cells "
          f"(board-downstream {dense['board_downstream_cells']}, post-tip {dense['post_tip_cells']})")

    if 'board_annulus_pattern' in out:
        a = out['board_annulus_pattern']
        print("board-downstream annulus INCLUDING passages (r 8-18 mm, first 20 mm):")
        for k, v in a.items():
            print(f"  {k:5s} mean={v['mean']:12.5g} |m=1|={v['m1']:10.4g} |m=3|={v['m3']:10.4g} "
                  f"(m3/mean {v['m3_over_mean']:.4f}, m3/m1 {v['m3_over_m1']:.1f}, phase {v['phase_deg']:.1f} deg)")
    if 'board_core_motion' in out:
        c = out['board_core_motion']
        print(f"board-core INCLUDING torch passage: volume-integrated kinetic moments [J] "
              f"Mzz={c['Mzz']:.4f} Mrr={c['Mrr']:.5f} Mtt={c['Mtt']:.5f} "
              f"-> Mrr/Mzz={c['radial_to_axial_motion_ratio']:.4f} "
              '(NOT vector momentum cancellation)')
        print(f"  speed {c['speed_mean']:.2f}+-{c['speed_std']:.2f} m/s (radial rms {c['radial_speed_rms']:.2f}); "
              f"Y_O2 {c['Y_mean']:.4f}+-{c['Y_std']:.4f} (min {c['Y_min']:.4f}); "
              f"sigmaT={c['T_std']:.2f} K, sigmaP={c['p_std_MPa']:.4f} MPa")
    if args.json:
        json.dump(out, open(args.json, 'w'), indent=1)
        print('wrote', args.json)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
