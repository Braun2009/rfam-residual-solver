#!/usr/bin/env python3
"""Design-consistency audit of the injector feed data (2026-09-17 input).

Inputs (user-supplied design table):
  LOX   T=130 K, p=6.0 MPa, rho=934 kg/m3, mu=8.58e-5 Pa s, mdot=2.8  kg/s
  C3H8  T=450 K, p=4.5 MPa, rho=66.6 kg/m3, mu=1.40e-5 Pa s, mdot=1.05 kg/s
        (total 1.5 = 1.05 main combustion + 0.45 liquid film)
  both Cd = 0.75

The geometric areas are the CAD injector-plate inlet patches used by the B
mesh (three LOX zones, three C3H8 zones, no film slots). For each propellant
this script:
  1. verifies rho/h/s/a against CoolProp 7.2.0 HEOS at the stated p/T,
  2. computes the isentropic mass-flow curve through the geometric area,
  3. reports the mass flow at a grid of discharge pressures with Cd=1 and
     Cd=0.75, and
  4. inverts Cd*mdot_isen(p2)=design mdot to get the injector-face pressure
     implied by each propellant.
No CAD, throat, or design-flow value is modified.
"""
import argparse
import json
from pathlib import Path

from CoolProp import AbstractState
from CoolProp import CoolProp as CP
from CoolProp.CoolProp import PhaseSI

FLUIDS = {
    'lox': dict(fluid='Oxygen', T=130.0, p=6.0e6, rho=934.0, mu=8.58e-5,
                mdot=2.8, area=3*18.28284446549178e-6, Cd=0.75),
    'fuel': dict(fluid='Propane', T=450.0, p=4.5e6, rho=66.6, mu=1.40e-5,
                 mdot=1.05, area=3*60.47183163113513e-6, Cd=0.75),
}


def isentropic_flow(state, p0, t0, p_back, area):
    """mdot through `area` when the static state is on the supply isentrope."""
    state.update(CP.PT_INPUTS, p0, t0)
    h0, s0 = state.hmass(), state.smass()
    state.update(CP.PSmass_INPUTS, p_back, s0)
    h, rho = state.hmass(), state.rhomass()
    kinetic = 2.0*max(h0 - h, 0.0)
    if kinetic <= 0.0:
        return None
    u = kinetic**0.5
    return dict(p_back_Pa=p_back, mdot_kg_s=rho*area*u, rho=rho, T=state.T(),
                u_isen=u, s_change=state.smass() - s0)


def choked_state(state, p0, t0):
    """Isentropic sonic point by maximising the mass flux."""
    state.update(CP.PT_INPUTS, p0, t0)
    h0, s0 = state.hmass(), state.smass()
    cache = {}

    def flux(p):
        if p in cache:
            return cache[p]
        try:
            state.update(CP.PSmass_INPUTS, p, s0)
            kinetic = 2.0*max(h0 - state.hmass(), 0.0)
            value = state.rhomass()*max(kinetic, 0.0)**0.5 if kinetic > 0 else 0.0
        except Exception:
            value = None
        cache[p] = value
        return value

    lo, hi = 1.0, p0
    for _ in range(200):
        m1, m2 = lo+(hi-lo)/3, hi-(hi-lo)/3
        f1, f2 = flux(m1), flux(m2)
        if f1 is None or f2 is None:
            break
        if f1 < f2:
            lo = m1
        else:
            hi = m2
    p = 0.5*(lo+hi)
    state.update(CP.PSmass_INPUTS, p, s0)
    return dict(p_crit_Pa=p, mdot_per_area=flux(p), T=state.T(),
                rho=state.rhomass(), s0=s0)


def solve_backpressure(state, p0, t0, area, target, Cd):
    """Find p_back with Cd*mdot_isen(p_back)=target (subsonic branch)."""
    lo, hi = 1.0, p0
    for _ in range(100):
        mid = 0.5*(lo+hi)
        flow = isentropic_flow(state, p0, t0, mid, area)
        if flow is None:
            hi = mid
            continue
        value = Cd*flow['mdot_kg_s']
        if value > target:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo+hi)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', default='work/a406_rootcause/injector_design_point.json')
    parser.add_argument('--grid', default='4.4,4.3,4.2,4.1,4.0,3.9,3.8,3.7,3.6,3.4,3.2,3.0')
    args = parser.parse_args()
    rows = []
    for name, spec in FLUIDS.items():
        state = AbstractState('HEOS', spec['fluid'])
        state.update(CP.PT_INPUTS, spec['p'], spec['T'])
        supply = dict(fluid=spec['fluid'], p_Pa=spec['p'], T_K=spec['T'],
                      rho_kg_m3=state.rhomass(), rho_table=spec['rho'],
                      rho_rel_error=(state.rhomass()-spec['rho'])/spec['rho'],
                      mu_Pa_s=state.viscosity(), mu_table=spec['mu'],
                      h_J_kg=state.hmass(), s_J_kgK=state.smass(),
                      a_m_s=state.speed_sound(), phase=str(state.phase()),
                      area_m2=spec['area'], mdot_design=spec['mdot'],
                      Cd=spec['Cd'])
        choked = choked_state(state, spec['p'], spec['T'])
        supply['choked'] = choked
        curve = [isentropic_flow(state, spec['p'], spec['T'], float(p)*1e6,
                                 spec['area']) for p in args.grid.split(',')]
        supply['curve'] = [dict(c, mdot_Cd=spec['Cd']*c['mdot_kg_s'])
                           for c in curve if c]
        supply['implied_face_pressure_Pa'] = solve_backpressure(
            state, spec['p'], spec['T'], spec['area'], spec['mdot'], spec['Cd'])
        supply['implied_face_pressure_no_Cd_Pa'] = solve_backpressure(
            state, spec['p'], spec['T'], spec['area'], spec['mdot'], 1.0)
        supply['mdot_if_Cd1_at_implied_p'] = isentropic_flow(
            state, spec['p'], spec['T'],
            supply['implied_face_pressure_Pa'], spec['area'])['mdot_kg_s']
        rows.append(supply)
        print('%s %s rho=%.4f (table %.1f, rel %.2e) mu=%.3e (table %.2e) '
              'a=%.1f phase=%s' % (name, spec['fluid'], supply['rho_kg_m3'],
              spec['rho'], supply['rho_rel_error'], supply['mu_Pa_s'],
              spec['mu'], supply['a_m_s'], supply['phase']))
        print('   choked p=%.3f MPa  mdot/A=%.1f  T=%.2f' % (
            choked['p_crit_Pa']/1e6, choked['mdot_per_area'], choked['T']))
        for c in supply['curve']:
            print('   p_back=%.2f MPa  mdot_isen=%.4f  mdot_Cd=%.4f  '
                  'u_isen=%.1f  T=%.2f' % (c['p_back_Pa']/1e6,
                  c['mdot_kg_s'], spec['Cd']*c['mdot_kg_s'], c['u_isen'], c['T']))
        print('   implied face pressure: Cd=%.2f -> %.4f MPa ; Cd=1 -> %.4f MPa' % (
            spec['Cd'], supply['implied_face_pressure_Pa']/1e6,
            supply['implied_face_pressure_no_Cd_Pa']/1e6))
    report = {'source': 'user-supplied injector feed table 2026-09-17',
              'assumption': 'supply p/T are stagnation conditions; steady '
                            'isentropic acceleration to the face static pressure; '
                            'Cd multiplies the isentropic mass flux',
              'note': '0.45 kg/s C3H8 film does not pass through the modelled '
                      'plate zones and is excluded',
              'rows': rows,
              'implied_pressure_spread_Pa': abs(
                  rows[0]['implied_face_pressure_Pa']
                  - rows[1]['implied_face_pressure_Pa'])}
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(args.out)


if __name__ == '__main__':
    main()
