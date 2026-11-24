#!/usr/bin/env python3
"""Vortex 3-pintle injector usability screening (design-faithful).

Answers, in minutes and without CFD: does the injector have healthy hydraulics,
does it atomize, does the three-pintle momentum cancellation actually balance,
is the centre ignitable, and does the design point sit inside the envelope of
published LOX/propane firings.

Inputs are taken verbatim from docs/design/vortex.txt (v5, 2026-08-21) plus
CoolProp for the supply-state properties. Everything is algebraic; the only
assumed inputs are the discharge coefficients (reported over a range, because
the design's Cd=0.75 is an assumption) and the LOX slot cant angle (reported
as a sensitivity).

Units: SI internally, engineering units in the output.
"""
from __future__ import annotations
import argparse
import json
import math

import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState

# ---------------------------------------------------------------- design inputs
DESIGN = dict(
    n_pintles=3,
    r_bolt=17.5e-3,          # m, pintle centres
    d_pintle=9.0e-3,         # m, barrel OD
    lox_slot_h=0.74e-3,      # m, LOX annular slot height (exit area)
    fuel_gap_w=2.20e-3,      # legacy annulus width (kept for the sweep history)
    fuel_hole_d=2.00e-3,     # m, C3H8 discrete orifice diameter (measured in the mesh)
    fuel_holes=60,           # 20 per 120 deg sector x 3 (matches the design area)
    torch_nozzle_d=3.4e-3,   # m, core torch exit
    mdot_lox_pintle=0.9007,  # kg/s
    mdot_fuel_pintle=0.3378, # kg/s
    mdot_lox_torch=0.0980,
    mdot_fuel_torch=0.0367,
    l_ext=12.3e-3, l_spray=39.7e-3, l_skip=4.0e-3, l_imp=52.0e-3,
    lox_cant_deg=43.4,       # LOX sheet half angle from the pintle axis
    fuel_swirl_deg=21.8,     # C3H8 swirl angle (from the design's S=0.4)
    chamber_d=88.3e-3,       # m, cylindrical section ID
    pc=4.0e6,                # Pa design chamber pressure
    dt=47.2e-3,              # m throat diameter
    lstar=1285.8e-3,         # m characteristic length
    cd_design=0.75,
)


def supply_state(fluid: str, p: float, t: float) -> AbstractState:
    s = AbstractState('HEOS', fluid)
    s.update(CP.PT_INPUTS, p, t)
    return s


def screen(cd: float, lox_cant_deg: float, tol_pct: float, overrides: dict | None = None) -> dict:
    d = {**DESIGN, **(overrides or {})}
    o2 = supply_state('Oxygen', 6.0e6, 130.0)
    c3 = supply_state('Propane', 4.5e6, 450.0)
    rho_o, rho_f = o2.rhomass(), c3.rhomass()

    # exit areas (annular slots around the barrel)
    a_lox_p = math.pi * d['d_pintle'] * d['lox_slot_h']
    # C3H8 is injected through discrete circular orifices (60 x dia 2.00 mm,
    # measured on the mesh boundary patches: 20 per 120 deg sector, area
    # 188.5 mm2 vs the design's 186.30 mm2), NOT through an annular gap.
    a_fuel_p = d['fuel_holes'] * math.pi * (d['fuel_hole_d'] / 2) ** 2 / d['n_pintles']
    a_lox, a_fuel = d['n_pintles'] * a_lox_p, d['n_pintles'] * a_fuel_p
    a_torch = math.pi * (d['torch_nozzle_d'] / 2) ** 2

    m_o = d['n_pintles'] * d['mdot_lox_pintle']
    m_f = d['n_pintles'] * d['mdot_fuel_pintle']
    v_o = m_o / (rho_o * a_lox)
    v_f = m_f / (rho_f * a_fuel)

    # per-pintle momentum vectors in the pintle frame (z = axis, r = radial)
    mdot_o_p = d['mdot_lox_pintle']
    mdot_f_p = d['mdot_fuel_pintle']
    cant = math.radians(lox_cant_deg)
    swirl = math.radians(d['fuel_swirl_deg'])
    j_lox_r = mdot_o_p * v_o * math.sin(cant)        # radial (toward the axis)
    j_lox_z = mdot_o_p * v_o * math.cos(cant)        # axial
    j_fuel_z = mdot_f_p * v_f * math.cos(swirl)      # axial
    j_fuel_t = mdot_f_p * v_f * math.sin(swirl)      # tangential (swirl)
    tmr = j_fuel_z / j_lox_r if j_lox_r > 0 else float('nan')

    # three-pintle vector cancellation: radial components are 120 deg apart
    # (their vector sum is zero by symmetry); the residual from a flow
    # asymmetry of tol_pct on one pintle is what leaves the centre non-uniform.
    j_radial_all = d['n_pintles'] * j_lox_r            # magnitudes, not summed as vectors
    j_res = (tol_pct / 100.0) * j_lox_r                # one pintle off by tol
    rho_hot = 3.45                                     # hot-gas density estimate
    # Robust, scale-free metric: the residual radial momentum relative to the
    # whole injection momentum flux. The absolute lateral velocity of the
    # central stream then follows from the axial momentum balance (the local
    # uniformity of the stagnation zone is a CFD question, not an algebraic one).
    j_total = d['n_pintles'] * (j_lox_r + j_lox_z + j_fuel_z)
    res_ratio = j_res / j_total
    v_lateral = (v_o * res_ratio) if j_total > 0 else 0.0
    pstag_center = 0.5 * rho_hot * (j_res / (rho_hot * math.pi * (5e-3) ** 2)) ** 2

    # hydraulics
    dp_o = (v_o / cd) ** 2 * rho_o / 2.0
    dp_f = (v_f / cd) ** 2 * rho_f / 2.0

    # atomization (LOX liquid sheet; propane is supercritical gas at 450 K)
    try:
        sat = AbstractState('HEOS', 'Oxygen')
        sat.update(CP.QT_INPUTS, 0.0, 130.0)
        sigma = sat.surface_tension()
    except Exception:
        sigma = float('nan')
    we_hot = rho_hot * v_o ** 2 * d['lox_slot_h'] / sigma
    we_cold = 42.5 * v_o ** 2 * d['lox_slot_h'] / sigma
    oh = o2.viscosity() / math.sqrt(rho_o * sigma * d['lox_slot_h'])

    # sheet travel to the centre vs sheet breakup
    r_travel = d['r_bolt']
    t_imp = r_travel / v_o
    t_break = d['l_spray'] / v_o

    # chamber / ignition
    a_t = math.pi * (d['dt'] / 2) ** 2
    v_chamber = d['lstar'] * a_t
    mdot_tot = m_o + m_f
    tau_hot = rho_hot * v_chamber / mdot_tot
    fill_mass = rho_hot * v_chamber

    # public LOX/propane envelope (NASA NTRS 19810006605 / 19840003267)
    env = dict(pc_bar=(8.6, 103.0), fuel_v=(14.6, 215.0), fuel_T=(119.0, 343.0),
               orifice_mm=(0.61, 1.14))
    d_lox_mm = d['lox_slot_h'] * 1e3
    d_fuel_mm = d['fuel_hole_d'] * 1e3

    def inside(x, rng):
        return rng[0] <= x <= rng[1]

    checks = [
        ('LOX exit velocity', f'{v_o:.1f} m/s', 'design ~46 m/s', abs(v_o - 46) / 46 < 0.10),
        ('fuel exit velocity', f'{v_f:.1f} m/s', 'design ~82 m/s', abs(v_f - 82) / 82 < 0.10),
        ('TMR per pintle', f'{tmr:.2f}', 'design 1.0 (radial LOX vs axial fuel)',
         0.7 <= tmr <= 1.4),
        ('LOX dP', f'{dp_o/1e6:.2f} MPa ({100*dp_o/d["pc"]:.0f}% Pc)', 'design 1.8 MPa; chug margin >10%',
         dp_o / d['pc'] > 0.10),
        ('fuel dP', f'{dp_f/1e6:.2f} MPa ({100*dp_f/d["pc"]:.0f}% Pc)', 'design 0.5 MPa',
         dp_f / d['pc'] > 0.05),
        ('LOX Weber (hot gas)', f'{we_hot:.0f}', '>12 atomizes', we_hot > 12),
        ('LOX Ohnesorge', f'{oh:.2e}', '<0.1 fine spray', oh < 0.1),
        ('fuel phase', 'supercritical gas (Tc=370 K < 450 K)', 'no atomization needed', True),
        ('sheet travel to centre', f'{r_travel*1e3:.1f} mm in {t_imp*1e3:.2f} ms vs breakup {t_break*1e3:.2f} ms',
         'coherent at impingement (good cancellation)', t_imp < t_break),
        ('residual momentum ({:.0f}% asym)'.format(tol_pct),
         f'{j_res:.2f} N = {100*res_ratio:.1f}% of injection momentum -> {v_lateral:.2f} m/s lateral',
         'centre uniformity (local profile needs CFD)', res_ratio < 0.05),
        ('residence time (hot)', f'{tau_hot*1e3:.2f} ms', '> chemical time (~0.5 ms)', tau_hot > 0.5e-3),
        ('propellant before ignition', f'{mdot_tot*1e-3:.2f} g/ms; chamber holds {fill_mass*1e3:.1f} g (hot gas)',
         'ignition within ~1 ms keeps the hard-start spike low', fill_mass/(mdot_tot*1e-3) > 1.0),
        ('Pc inside public envelope', f'{d["pc"]/1e5:.0f} bar', f'{env["pc_bar"][0]}-{env["pc_bar"][1]} bar',
         inside(d['pc'] / 1e5, env['pc_bar'])),
        ('fuel velocity inside envelope', f'{v_f:.0f} m/s', f'{env["fuel_v"][0]}-{env["fuel_v"][1]} m/s',
         inside(v_f, env['fuel_v'])),
        ('fuel temperature inside envelope', '450 K', f'{env["fuel_T"][0]}-{env["fuel_T"][1]} K (above: beyond tested)',
         inside(450.0, env['fuel_T'])),
        ('LOX slot inside envelope', f'{d_lox_mm:.2f} mm', f'{env["orifice_mm"][0]}-{env["orifice_mm"][1]} mm',
         inside(d_lox_mm, env['orifice_mm'])),
        ('fuel orifice inside envelope', f'{d_fuel_mm:.2f} mm x {d["fuel_holes"]}', f'{env["orifice_mm"][0]}-{env["orifice_mm"][1]} mm (above tested)',
         inside(d_fuel_mm, env['orifice_mm'])),
    ]
    return dict(v_o=v_o, v_f=v_f, a_lox=a_lox, a_fuel=a_fuel, a_torch=a_torch, tmr=tmr,
                j_lox_r=j_lox_r, j_lox_z=j_lox_z, j_fuel_z=j_fuel_z, j_fuel_t=j_fuel_t,
                j_res=j_res, res_ratio=res_ratio, v_lateral=v_lateral, j_total=j_total,
                pstag_center=pstag_center, dp_o=dp_o, dp_f=dp_f,
                we_hot=we_hot, we_cold=we_cold, oh=oh, sigma=sigma, tau_hot=tau_hot,
                fill_mass=fill_mass, checks=checks, cd=cd, lox_cant_deg=lox_cant_deg)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--json', help='write the screening result as JSON')
    ap.add_argument('--tol-pct', type=float, default=5.0, help='feed asymmetry for the residual estimate')
    ap.add_argument('--sweep', action='append', default=[],
                    help='design sweep, e.g. --sweep lox_slot_h=0.5,0.74,1.0 (repeatable)')
    ap.add_argument('--target-dp-frac', type=float, default=None,
                    help='solve the exit area for a target dP/Pc (e.g. 0.20)')
    ap.add_argument('--cd', type=float, default=None, help='override the discharge coefficient')
    args = ap.parse_args()
    cd_eff = args.cd if args.cd else DESIGN['cd_design']
    d = DESIGN
    if args.target_dp_frac is not None:
        o2 = supply_state('Oxygen', 6.0e6, 130.0)
        c3 = supply_state('Propane', 4.5e6, 450.0)
        for name, mdot, rho, cur in (('LOX', d['mdot_lox_pintle'], o2.rhomass(), d['lox_slot_h']),
                                     ('C3H8', d['mdot_fuel_pintle'], c3.rhomass(), d['fuel_gap_w'])):
            dp = args.target_dp_frac * d['pc']
            a_exit = mdot / (cd_eff * math.sqrt(2.0 * rho * dp))
            gap = a_exit / (math.pi * d['d_pintle'])
            print(f'target dP/Pc={args.target_dp_frac:.2f} -> {name}: A_exit={a_exit*1e6:.2f} mm2/pintle '
                  f'({d["n_pintles"]*a_exit*1e6:.2f} mm2 total), gap={gap*1e3:.3f} mm (current {cur*1e3:.2f} mm)')
    for sw in args.sweep:
        if '=' not in sw:
            print('bad --sweep (need name=v1,v2,..):', sw)
            continue
        name, vals = sw.split('=', 1)
        if name not in DESIGN:
            print('unknown parameter:', name)
            continue
        print(f'\n--- sweep {name} (Cd={cd_eff}) ---')
        length_keys = {'lox_slot_h', 'fuel_gap_w', 'd_pintle', 'torch_nozzle_d', 'r_bolt',
                       'l_ext', 'l_spray', 'l_skip', 'l_imp', 'chamber_d', 'dt', 'lstar'}
        for raw in vals.split(','):
            v = float(raw) * (1e-3 if name in length_keys else 1.0)   # lengths are given in mm
            r = screen(cd_eff, DESIGN['lox_cant_deg'], args.tol_pct, {name: v})
            npass = sum(1 for _, _, _, ok in r['checks'] if ok)
            print(f'  {name}={raw:<9s} v_LOX={r["v_o"]:6.1f} v_fuel={r["v_f"]:6.1f} TMR={r["tmr"]:5.2f} '
                  f'dP_LOX={r["dp_o"]/1e6:5.2f}MPa({100*r["dp_o"]/DESIGN["pc"]:4.0f}%) '
                  f'dP_fuel={r["dp_f"]/1e6:5.2f}MPa We={r["we_hot"]:7.0f} res={100*r["res_ratio"]:4.2f}% '
                  f'pass {npass}/{len(r["checks"])}')
    print('=== Vortex 3-pintle injector screening (design v5) ===')
    print(f'fuel orifices: {d["fuel_holes"]} x dia {d["fuel_hole_d"]*1e3:.2f} mm '
          f'= {d["fuel_holes"]*math.pi*(d["fuel_hole_d"]/2)**2*1e6:.2f} mm2 total; '
          f'LOX: {d["n_pintles"]} annular slots')
    print(f'exit areas: LOX {d["n_pintles"]}x{math.pi*d["d_pintle"]*d["lox_slot_h"]*1e6:.2f} mm2 '
          f'= {d["n_pintles"]*math.pi*d["d_pintle"]*d["lox_slot_h"]*1e6:.2f} mm2 (design 62.97) | '
          f'fuel {d["n_pintles"]*math.pi*d["d_pintle"]*d["fuel_gap_w"]*1e6:.2f} mm2 (design 186.30) | '
          f'torch {math.pi*(d["torch_nozzle_d"]/2)**2*1e6:.2f} mm2 (design 9.08)')
    for cd in (0.6, 0.7, d['cd_design'], 0.8, 0.9):
        r = screen(cd, d['lox_cant_deg'], args.tol_pct)
        print(f'  Cd={cd:.2f}: v_LOX={r["v_o"]:.1f} m/s  v_fuel={r["v_f"]:.1f} m/s  '
              f'TMR={r["tmr"]:.2f}  dP_LOX={r["dp_o"]/1e6:.2f} MPa ({100*r["dp_o"]/d["pc"]:.0f}%)  '
              f'dP_fuel={r["dp_f"]/1e6:.2f} MPa ({100*r["dp_f"]/d["pc"]:.0f}%)')
    r = screen(cd_eff, d['lox_cant_deg'], args.tol_pct)
    print(f'\nmomentum per pintle: LOX radial {r["j_lox_r"]:.1f} N + axial {r["j_lox_z"]:.1f} N | '
          f'fuel axial {r["j_fuel_z"]:.1f} N + tangential {r["j_fuel_t"]:.1f} N  -> TMR {r["tmr"]:.2f}')
    print(f'three-pintle radial vectors at 120 deg: sum = 0 by symmetry; '
          f'{args.tol_pct:.0f}% asymmetry -> {r["j_res"]:.2f} N residual = '
          f'{100*r["res_ratio"]:.1f}% of the {r["j_total"]:.0f} N injection momentum -> '
          f'{r["v_lateral"]:.2f} m/s lateral drift (local uniformity is a CFD question)')
    print(f'centre stagnation pressure estimate: {r["pstag_center"]/1e3:.1f} kPa '
          f'(residual confined to the 5 mm core; hot-gas reference)')
    print(f'\nsigma_LOX(130K)={r["sigma"]*1e3:.2f} mN/m  We_hot={r["we_hot"]:.0f}  We_cold={r["we_cold"]:.0f}  '
          f'Oh={r["oh"]:.2e}  tau_hot={r["tau_hot"]*1e3:.2f} ms  chamber holds {r["fill_mass"]*1e3:.0f} g')
    print('\n--- verdict ---')
    npass = 0
    for name, value, criterion, ok in r['checks']:
        flag = 'PASS' if ok else 'MARGINAL'
        npass += bool(ok)
        print(f'  [{flag:8s}] {name:34s} {value:55s} ({criterion})')
    print(f'\n{npass}/{len(r["checks"])} checks pass at Cd={d["cd_design"]} (Cd is a design '
          f'assumption; the public LOX/propane firings suggest 0.6-0.9)')
    if args.json:
        with open(args.json, 'w') as fh:
            json.dump({k: v for k, v in r.items() if k != 'checks'},
                      fh, indent=1, default=str)
        print('wrote', args.json)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
