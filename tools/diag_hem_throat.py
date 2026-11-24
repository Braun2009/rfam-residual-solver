#!/usr/bin/env python3
"""WITHDRAWN: historical unvalidated HEM capacity estimate.

This script skips flash failures, interpolates properties without closure
verification, and has no grid/sonic-maximum convergence check. Its output
must not be used as a validated operating point or hardware recommendation.
See docs/coldflow_diagnosis_recheck.md. Preserved for audit reproduction only.
"""
import argparse
import math
import sys

from CoolProp import AbstractState, CoolProp as CP

Y = [2.7021 / 3.7155, 1.0134 / 3.7155]


def make_state():
    state = AbstractState('HEOS', 'Oxygen&Propane')
    state.set_mass_fractions(Y)
    return state


def evaluate(state, pressure, temperature):
    try:
        state.update(CP.PT_INPUTS, pressure, temperature)
        return state.hmass(), state.smass(), state.rhomass()
    except Exception:
        return None


def interpolate(samples, key, target):
    for a, b in zip(samples[:-1], samples[1:]):
        if (a[key] - target) * (b[key] - target) <= 0 and a[key] != b[key]:
            weight = (target - a[key]) / (b[key] - a[key])
            return tuple(a[i] + weight * (b[i] - a[i]) for i in range(len(a)))
    return None


def chamber_state(state, pressure, enthalpy):
    samples = []
    for index in range(361):
        temperature = 140.0 + index
        row = evaluate(state, pressure, temperature)
        if row is not None:
            samples.append((temperature,) + row)
    return interpolate(samples, 1, enthalpy)


def throat(state, pressure0, enthalpy0):
    chamber = chamber_state(state, pressure0, enthalpy0)
    if chamber is None:
        return None
    entropy = chamber[2]
    previous_temperature = chamber[0]
    pressure = pressure0 * 0.99
    maximum = None
    profile = []
    while pressure > 0.08e6:
        samples = []
        window = 6.0
        for widen in (1.0, 3.0, 8.0):
            narrow = previous_temperature - window * widen
            wide = previous_temperature + window * widen
            step = max(0.05, 0.25 * widen)
            for index in range(int((wide - narrow) / step) + 1):
                temperature = narrow + index * step
                row = evaluate(state, pressure, temperature)
                if row is not None:
                    samples.append((temperature,) + row)
            found = interpolate(samples, 2, entropy)
            if found is not None:
                break
        if found is None:
            break
        temperature, h, s, rho = found
        previous_temperature = temperature
        speed = math.sqrt(max(0.0, 2.0 * (enthalpy0 - h)))
        flux = rho * speed
        profile.append((flux, pressure, rho, speed, temperature, h))
        if maximum is None or flux > maximum[0]:
            maximum = profile[-1]
        pressure *= 0.975
    return entropy, maximum, profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--historical-unvalidated', action='store_true',
                        help='explicitly reproduce the withdrawn calculation, not a design result')
    args = parser.parse_args()
    if not args.historical_unvalidated:
        parser.error('WITHDRAWN calculation; see docs/coldflow_diagnosis_recheck.md')
    print('WITHDRAWN_UNVALIDATED_ESTIMATE: not a flow/phase/design conclusion', flush=True)
    lox = AbstractState('HEOS', 'Oxygen')
    lox.update(CP.PT_INPUTS, 6e6, 133.5)
    fuel = AbstractState('HEOS', 'Propane')
    fuel.update(CP.PT_INPUTS, 4.5e6, 453.5)
    mass_flow_lox, mass_flow_fuel = 2.7021, 1.0134
    enthalpy0 = (mass_flow_lox * lox.hmass() + mass_flow_fuel * fuel.hmass()) / (mass_flow_lox + mass_flow_fuel)
    print('h0mix=%.1f J/kg  (CoolProp absolute enthalpy)' % enthalpy0)
    for throat_radius in [23.6e-3]:
        area = math.pi * throat_radius ** 2
        print('throat area=%.6f m2' % area)
        for pressure0 in [4e6, 3e6, 2e6, 1.5e6, 1e6]:
            state = make_state()
            result = throat(state, pressure0, enthalpy0)
            injected = 2.702 * math.sqrt(max(0.0, (6e6 - pressure0) / 1.8e6)) \
                + 1.0134 * math.sqrt(max(0.0, (4.5e6 - pressure0) / 0.5e6))
            if result is None:
                print('p0=%8.3g  no isentropic state found' % pressure0)
                continue
            entropy, maximum, profile = result
            if maximum is None:
                print('UNVALIDATED p0=%8.3g no usable flux samples' % pressure0)
                continue
            flux, sample_p, rho, sample_speed, temperature, h = maximum
            print('UNVALIDATED p0=%8.3g s=%.1f sampled_Gmax=%8.0f sampled_mdot=%7.3f kg/s '
                  'sample_p=%.3g rho=%.2f speed_not_verified_sound=%.1f T=%.1f scaled_inlet_mdot=%7.3f'
                  % (pressure0, entropy, flux, flux * area, sample_p, rho, sample_speed,
                     temperature, injected))
            sys.stdout.flush()


if __name__ == '__main__':
    main()
