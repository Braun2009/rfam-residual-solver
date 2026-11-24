#!/usr/bin/env python3
"""WITHDRAWN: historical, internally inconsistent operating-point estimate.

The throat curve keeps design Y and h0 while the inlet flows vary with
pressure, changing both Y and h0. No self-consistent coupling is solved.
Preserved only for audit; see docs/coldflow_diagnosis_recheck.md.

Historical calculation: operating-point balance from HEM capacity samples.

Input Gmax values are produced by tools/diag_hem_throat.py (CoolProp 7.2.0
equilibrium isentropes, mass-averaged propellant stagnation enthalpy, real
nozzle throat r=23.6 mm). Injector flows are the design orifice fits:
  LOX  per patch: 0.9007 kg/s at delta_p=1.8 MPa (supply 6.0 MPa)
  C3H8 per patch: 0.3378 kg/s at delta_p=0.5 MPa (supply 4.5 MPa)
Both are scaled by sqrt(delta_p/design delta_p). This is a model comparison,
not a design verdict.
"""
import argparse
import json
import math
from pathlib import Path

THROAT_RADIUS = 23.6e-3
THROAT_AREA = math.pi * THROAT_RADIUS ** 2
# CoolProp HEM capacity points (pressure Pa, Gmax kg/m2/s), see diag_hem_throat.txt
CAPACITY = [(4e6, 12048.0), (3e6, 9042.0), (2e6, 6037.0), (1.5e6, 4536.0), (1e6, 3033.0)]


def capacity_per_pa():
    return sum(flux / pressure for pressure, flux in CAPACITY) / len(CAPACITY)


def throat_flow(pressure):
    return capacity_per_pa() * THROAT_AREA * pressure


def injector_flow(pressure):
    return 2.7021 * math.sqrt(max(0.0, (6e6 - pressure) / 1.8e6)) \
        + 1.0134 * math.sqrt(max(0.0, (4.5e6 - pressure) / 0.5e6))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--historical-unvalidated', action='store_true',
                        help='reproduce withdrawn arithmetic without treating it as self-consistent')
    args = parser.parse_args()
    if not args.historical_unvalidated:
        parser.error('WITHDRAWN: inlet composition/enthalpy and throat model are inconsistent')
    print('WITHDRAWN_INCONSISTENT_ESTIMATE: no operating point or hardware recommendation')
    slope = capacity_per_pa()
    low, high = 0.3e6, 4e6
    for _ in range(100):
        middle = 0.5 * (low + high)
        if injector_flow(middle) > throat_flow(middle):
            low = middle
        else:
            high = middle
    equilibrium = 0.5 * (low + high)

    def velocities(pressure):
        lox = math.sqrt(max(0.0, (6e6 - pressure) / 1.8e6))
        fuel = math.sqrt(max(0.0, (4.5e6 - pressure) / 0.5e6))
        return {'lox_m_s': 47.19 * lox, 'fuel_axial_m_s': 86.03 * fuel,
                'VR': (47.19 * lox) / (86.03 * fuel)}

    restrictor_area = 3.7155 / (slope * 4.0e6)
    report = {
        'status': 'WITHDRAWN_INCONSISTENT_ESTIMATE',
        'validated': False,
        'design_use_allowed': False,
        'known_error': 'Throat curve has fixed Y/h0; pressure-dependent inlet flows require different Y/h0.',
        'reference': 'docs/coldflow_diagnosis_recheck.md',
        'coolprop': '7.2.0',
        'throat_area_m2': THROAT_AREA,
        'hem_capacity_slope_kg_m2s_per_Pa': slope,
        'design_point': {'chamber_p_Pa': 4.0e6, 'injector_flow_kg_s': injector_flow(4.0e6),
                         'throat_capacity_kg_s': throat_flow(4.0e6),
                         'hot_gas_c_star_m_s': 4.0e6 * THROAT_AREA / 3.7155,
                         'velocities': velocities(4.0e6)},
        'fixed_flow_choke_pressure_Pa': 3.7155 / (slope * THROAT_AREA),
        'self_consistent_cold_point': {'chamber_p_Pa': equilibrium,
                                       'mass_flow_kg_s': injector_flow(equilibrium),
                                       'velocities': velocities(equilibrium)},
        'cold_restrictor_for_design_point': {
            'area_m2': restrictor_area,
            'radius_mm': math.sqrt(restrictor_area / math.pi) * 1e3,
            'diameter_mm': 2.0 * math.sqrt(restrictor_area / math.pi) * 1e3},
    }
    out = Path(__file__).resolve().parents[1] / 'work/coldflow_B_verified/diagnosis_recheck/historical_balance_unvalidated.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('HEM capacity slope: %.4f kg/(m2 s Pa)  =>  %.1f kg/s per MPa at the real throat'
          % (slope, slope * THROAT_AREA * 1e6))
    print('design point 4.000 MPa: injector %.3f kg/s, throat capacity %.3f kg/s, hot-gas c* %.0f m/s'
          % (report['design_point']['injector_flow_kg_s'], report['design_point']['throat_capacity_kg_s'],
             report['design_point']['hot_gas_c_star_m_s']))
    print('UNVALIDATED fixed-flow linear extrapolation: %.3f MPa, NOT verified chamber pressure'
          % (report['fixed_flow_choke_pressure_Pa'] / 1e6))
    point = report['self_consistent_cold_point']
    print('WITHDRAWN intersection (not self-consistent): %.3f MPa, %.3f kg/s, scaled LOX %.1f m/s, fuel %.1f m/s, VR %.3f (reference %.3f)'
          % (point['chamber_p_Pa'] / 1e6, point['mass_flow_kg_s'],
             point['velocities']['lox_m_s'], point['velocities']['fuel_axial_m_s'],
             point['velocities']['VR'], 0.5485))
    restrictor = report['cold_restrictor_for_design_point']
    print('WITHDRAWN area conversion: d=%.1f mm (reference throat d=%.1f mm); NOT a design recommendation'
          % (restrictor['diameter_mm'], 2 * THROAT_RADIUS * 1e3))
    print('wrote', out)


if __name__ == '__main__':
    main()
