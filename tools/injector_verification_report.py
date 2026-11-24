#!/usr/bin/env python3
"""Consolidated injector verification record (Vortex 3-pintle, design v5).

Runs the design-faithful screening, maps it onto the published LOX/propane test
envelope, and -- when a cold-flow dump is given -- attaches the CFD evidence
(injector-channel priming state, composition front positions, near-face
velocity samples, and cold/dense proxies). It does not certify a converged
field, a CAD exit flux or an EOS two-phase inventory without a paired sidecar.

The point is traceability: every number carries its source (design file, mesh,
CoolProp, public test report, or a specific dump), so the "is this injector
usable?" verdict can be audited and re-run.
"""
from __future__ import annotations
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import injector_screening as S  # noqa: E402
from tetra_field_to_vtu import read_mesh, read_field  # noqa: E402

# Published LOX/propane test envelope (NASA NTRS 19810006605 / 19840003267)
ENVELOPE = dict(source='NASA NTRS 19810006605 (Aerojet/NASA-JSC NAS 9-15724, 1980) '
                        '+ 19840003267 (1983)',
                pc_bar=(8.6, 103.0), fuel_v_ms=(14.6, 215.0), fuel_T_K=(119.0, 343.0),
                orifice_mm=(0.61, 1.14), tests=127 + 77)


def coldflow_evidence(mesh_path: str, field_path: str, board_z: float) -> dict:
    _, cells, vol = read_mesh(mesh_path)
    cent, f, y, _ = read_field(field_path)
    p, T, u, v, w, rho, a, h = (f[:, i] for i in range(8))
    z = cent[:, 2]
    inner = z < board_z
    liq = rho > 500.
    o2 = y[:, 0] > 0.95
    fuel = y[:, 0] < 0.2
    out = dict(cells=int(len(cells)),
               density_proxy_definition='rho>500 kg/m3; not an EOS phase',
               liquid_cells=int(liq.sum()),
               liquid_mass_kg=float(np.sum(rho[liq] * vol[liq])),
               pure_o2_front_z=float(z[o2].max()) if o2.any() else None,
               fuel_cells=int(fuel.sum()),
               fuel_mass_kg=float(np.sum(rho[fuel] * vol[fuel])),
               fuel_front_z=float(z[fuel].max()) if fuel.any() else None,
               board_z=board_z,
               two_phase_like_cells=int(np.sum((T < 300.) & (y[:, 0] > 0.7))),
               two_phase_like_definition='T<300 K and Y_O2>0.7; not EOS phase')
    status_path=Path(field_path+'.status.json')
    status=json.loads(status_path.read_text()) if status_path.exists() else None
    out['solver_converged']=status.get('converged') if status else None
    out['production_validated']=status.get('production_validated') if status else None
    phase_path=Path(field_path+'.phase.csv')
    if phase_path.exists():
        phase=np.genfromtxt(phase_path,delimiter=',',skip_header=1,usecols=(0,1),dtype=np.int32)
        if len(phase)!=len(cells) or not np.array_equal(phase[:,0],np.arange(len(cells))):
            raise ValueError('phase sidecar/mesh cell order differs')
        out['eos_two_phase_cells']=int(np.count_nonzero(phase[:,1]==5))
        out['eos_two_phase_source']=str(phase_path)
    else:
        out['eos_two_phase_cells']=None
        out['eos_two_phase_source']=None
    # Board-downstream near-field includes protruding injector passages.
    near_board = (~inner) & (z < board_z + 2e-3)
    if near_board.any():
        out['board_downstream_includes_passages']=True
        out['face_p_MPa'] = float(p[near_board].mean() / 1e6)
        out['face_T_K'] = float(T[near_board].mean())
        out['face_YO2'] = float(y[near_board, 0].mean())
        out['face_uz_ms'] = float(w[near_board].mean())
    post_tip=z>=0.2858028
    out['post_tip_pure_O2_cells']=int(((y[:,0]>0.95)&post_tip).sum())
    out['post_tip_fuel_rich_cells']=int(((y[:,0]<0.2)&post_tip).sum())
    # near-exit stream velocities
    for lab, msk in (('lox', (y[:, 0] > 0.9) & (z > board_z - 0.26e-3) & inner),
                     ('fuel', (y[:, 0] < 0.1) & (z > board_z - 0.55e-3) & inner)):
        if msk.sum() > 3:
            out[f'{lab}_exit_cells'] = int(msk.sum())
            out[f'{lab}_exit_rho'] = float(rho[msk].mean())
            out[f'{lab}_exit_uz'] = float(w[msk].mean())
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh', default='work/mesh_rebuild/delivery/mesh/fluid_volume.rfamtet')
    ap.add_argument('--field', default=None, help='optional cold-flow dump for the CFD evidence')
    ap.add_argument('--json', default='work/a406_rootcause/injector_verification.json')
    ap.add_argument('--cd', type=float, default=0.75)
    args = ap.parse_args()

    d = S.DESIGN
    r = S.screen(args.cd, d['lox_cant_deg'], 5.0)
    record = dict(
        design_source='docs/design/vortex.txt v5 (2026-08-21)',
        geometry=dict(n_pintles=d['n_pintles'], r_bolt_mm=d['r_bolt'] * 1e3,
                      d_pintle_mm=d['d_pintle'] * 1e3,
                      lox_slot_mm=d['lox_slot_h'] * 1e3,
                      fuel_orifice_d_mm=d['fuel_hole_d'] * 1e3,
                      fuel_orifice_count=d['fuel_holes'],
                      torch_nozzle_mm=d['torch_nozzle_d'] * 1e3,
                      lox_exit_area_mm2=d['n_pintles'] * math.pi * d['d_pintle'] * d['lox_slot_h'] * 1e6,
                      fuel_exit_area_mm2=d['fuel_holes'] * math.pi * (d['fuel_hole_d'] / 2) ** 2 * 1e6),
        hydraulics=dict(v_lox_ms=r['v_o'], v_fuel_ms=r['v_f'], tmr=r['tmr'],
                        dp_lox_MPa=r['dp_o'] / 1e6, dp_fuel_MPa=r['dp_f'] / 1e6,
                        dp_lox_frac_pc=r['dp_o'] / d['pc'], dp_fuel_frac_pc=r['dp_f'] / d['pc'],
                        cd_assumed=args.cd),
        atomization=dict(we_hot=r['we_hot'], we_cold=r['we_cold'], ohnesorge=r['oh'],
                         sigma_N_m=r['sigma'], fuel_phase='supercritical gas (Tc=370 K < 450 K)'),
        cancellation=dict(per_pintle_lox_radial_N=r['j_lox_r'], per_pintle_fuel_axial_N=r['j_fuel_z'],
                          residual_N_for_5pct=r['j_res'], residual_ratio=r['res_ratio'],
                          lateral_drift_ms=r['v_lateral'],
                          sheet_travel_to_axis_ms=d['r_bolt'] / r['v_o'] * 1e3,
                          sheet_breakup_ms=d['l_spray'] / r['v_o'] * 1e3),
        ignition=dict(torch_nozzle_mm=d['torch_nozzle_d'] * 1e3,
                      residence_ms=r['tau_hot'] * 1e3, chamber_hot_mass_g=r['fill_mass'] * 1e3,
                      accumulation_g_per_ms=3.85),
        checks=[dict(name=n, value=v, criterion=c, pass_=ok) for n, v, c, ok in r['checks']],
        envelope=ENVELOPE,
    )
    if args.field:
        record['coldflow'] = coldflow_evidence(args.mesh, args.field, d and 0.2709527893066406)
    json.dump(record, open(args.json, 'w'), indent=1)
    npass = sum(1 for c in record['checks'] if c['pass_'])
    print(f'=== injector verification record ({record["design_source"]}) ===')
    print(f'geometry: {d["n_pintles"]}x pintle R={record["geometry"]["r_bolt_mm"]:.1f} mm, '
          f'LOX slot {record["geometry"]["lox_slot_mm"]:.2f} mm '
          f'({record["geometry"]["lox_exit_area_mm2"]:.2f} mm2 exit), '
          f'fuel {record["geometry"]["fuel_orifice_count"]}x dia '
          f'{record["geometry"]["fuel_orifice_d_mm"]:.2f} mm '
          f'({record["geometry"]["fuel_exit_area_mm2"]:.2f} mm2 exit), '
          f'torch {record["geometry"]["torch_nozzle_mm"]:.2f} mm')
    print(f'hydraulics: v_LOX={r["v_o"]:.1f} v_fuel={r["v_f"]:.1f} m/s  TMR={r["tmr"]:.2f}  '
          f'dP={r["dp_o"]/1e6:.2f}/{r["dp_f"]/1e6:.2f} MPa ({100*r["dp_o"]/d["pc"]:.0f}%/'
          f'{100*r["dp_f"]/d["pc"]:.0f}% Pc) at Cd={args.cd}')
    print(f'atomization: We={r["we_hot"]:.0f} (hot) Oh={r["oh"]:.1e} sigma={r["sigma"]*1e3:.2f} mN/m; '
          f'{record["atomization"]["fuel_phase"]}')
    print(f'cancellation: per pintle LOX radial {r["j_lox_r"]:.1f} N vs fuel axial {r["j_fuel_z"]:.1f} N; '
          f'sheet reaches axis in {record["cancellation"]["sheet_travel_to_axis_ms"]:.3f} ms vs '
          f'breakup {record["cancellation"]["sheet_breakup_ms"]:.3f} ms')
    print(f'ignition: torch {record["ignition"]["torch_nozzle_mm"]:.2f} mm, residence '
          f'{record["ignition"]["residence_ms"]:.2f} ms, chamber holds '
          f'{record["ignition"]["chamber_hot_mass_g"]:.1f} g (hot) vs 3.85 g/ms inflow')
    if 'coldflow' in record:
        cf = record['coldflow']
        print(f'cold-flow diagnostic: cells={cf["cells"]} rho>500 proxy={cf["liquid_cells"]} '
              f'({cf["liquid_mass_kg"]*1e6:.3f} mg), max pure-O2 centroid z={cf["pure_o2_front_z"]:.6f} '
              f'(board {cf["board_z"]:.6f}, NOT chamber penetration), '
              f'max fuel-rich centroid z={cf["fuel_front_z"]:.6f}, '
              f'cold/O2-rich proxy {cf["two_phase_like_cells"]}, EOS two phase '
              f'{cf["eos_two_phase_cells"] if cf["eos_two_phase_cells"] is not None else "unavailable"}; '
              f'converged={cf["solver_converged"]}')
        if 'lox_exit_uz' in cf:
            print(f'  near-face samples (NOT integrated CAD exit flux): LOX rho={cf["lox_exit_rho"]:.1f} '
                  f'uz={cf["lox_exit_uz"]:.1f} m/s; fuel rho={cf.get("fuel_exit_rho", float("nan")):.1f} '
                  f'uz={cf.get("fuel_exit_uz", float("nan")):.1f} m/s')
        if 'face_p_MPa' in cf:
            print(f'  board-downstream near face (INCLUDES channels): p={cf["face_p_MPa"]:.4f} MPa '
                  f'T={cf["face_T_K"]:.1f} K Y_O2={cf["face_YO2"]:.4f} '
                  f'uz={cf["face_uz_ms"]:.1f} m/s; pure O2/fuel-rich post-tip cells '
                  f'{cf["post_tip_pure_O2_cells"]}/{cf["post_tip_fuel_rich_cells"]}')
    print(f'\nscreening (Cd={args.cd} assumed, not measured): {npass}/{len(record["checks"])} '
          f'design-point checks; NOT a CFD or ignition acceptance. wrote {args.json}')
    for c in record['checks']:
        if not c['pass_']:
            print(f'  MARGINAL: {c["name"]} = {c["value"]}  ({c["criterion"]})')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
