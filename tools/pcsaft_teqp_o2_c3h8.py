#!/usr/bin/env python3
"""Predictive PC-SAFT (teqp) for O2 + C3H8 — the two-phase uncertainty band.

Why: LOX + hydrocarbon mixtures are auto-oxidation/explosion systems, so no
binary VLE/density data and no fitted departure function exist (CoolProp/GERG use
ideal mixing, F=0, for this pair). The predictive route is a molecular EOS with a
*defaulted* binary interaction: PC-SAFT with kij=0. Pure-component parameters are
fitted to CoolProp pure-fluid data only, so no hazardous experiment is involved.

The PC-SAFT evaluation uses teqp (an independent, well-tested implementation).
This script
  1. validates the literature propane parameters against CoolProp,
  2. fits (m, sigma, eps/k) for O2 to CoolProp saturation densities,
  3. traces the O2/propane isobar with kij=0 and compares the bubble/dew
     temperatures with CoolProp's ideal-mixing HEOS at the design compositions.

Output: the temperature difference at the operating pressure is the model
uncertainty band for the near-critical two-phase region.
"""
from __future__ import annotations
import argparse
import json
import math

import numpy as np
from scipy.optimize import least_squares

import teqp
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState

# literature PC-SAFT parameters (Gross & Sadowski 2001)
LIT = {'Propane': dict(m=2.002, sigma=3.618, eps=208.11)}


def make_pure(m, sigma, eps):
    sc = teqp.SAFTCoeffs()
    sc.m = m
    sc.sigma_Angstrom = sigma
    sc.epsilon_over_k = eps
    return teqp.PCSAFTEOS([sc])


def make_mixture(p1, p2, kij=None):
    scs = []
    for p in (p1, p2):
        sc = teqp.SAFTCoeffs()
        sc.m, sc.sigma_Angstrom, sc.epsilon_over_k = p['m'], p['sigma'], p['eps']
        scs.append(sc)
    if kij is None:
        return teqp.PCSAFTEOS(scs)
    # kmat is a symmetric matrix with zeros on the diagonal
    kmat = np.zeros((2, 2))
    kmat[0, 1] = kmat[1, 0] = kij
    return teqp.PCSAFTEOS(scs, kmat)


def coolprop_sat(fluid, temps):
    st = AbstractState('HEOS', fluid)
    out = []
    for t in temps:
        try:
            st.update(CP.QT_INPUTS, 0.0, t)
            rl = st.rhomass() / st.molar_mass()   # mol/m3
            st.update(CP.QT_INPUTS, 1.0, t)
            rv = st.rhomass() / st.molar_mass()
            out.append((t, rl, rv))
        except Exception:
            pass
    return out


def fit_o2(tmin=100.0, tmax=150.0, n=9):
    data = coolprop_sat('Oxygen', np.linspace(tmin, tmax, n))

    def residuals(theta):
        m, sigma, eps = theta
        model = make_pure(m, sigma, eps)
        res = []
        for t, rl, rv in data:
            try:
                out = model.pure_VLE_T(t, rl, max(rv, 1e-3), 100, None)
                res += [(out[0] - rl) / rl, (out[1] - rv) / rv]
            except Exception:
                res += [1.0, 1.0]
        return np.array(res)

    sol = least_squares(residuals, [1.2, 3.1, 105.0],
                        bounds=([0.5, 2.0, 60.0], [3.0, 4.5, 250.0]),
                        xtol=1e-10, ftol=1e-10)
    m, sigma, eps = sol.x
    model = make_pure(m, sigma, eps)
    errs_l, errs_v = [], []
    for t, rl, rv in data:
        out = model.pure_VLE_T(t, rl, max(rv, 1e-3), 100, None)
        errs_l.append(abs(out[0] - rl) / rl)
        errs_v.append(abs(out[1] - rv) / rv)
    return dict(m=float(m), sigma=float(sigma), eps=float(eps),
                max_rhol_err=float(max(errs_l)), max_rhov_err=float(max(errs_v)),
                temperatures=[d[0] for d in data])


def coolprop_bulb(p, y_mass):
    s = AbstractState('HEOS', 'Oxygen&Propane')
    s.set_mass_fractions([y_mass, 1.0 - y_mass])
    out = {}
    for q, lab in ((0.0, 'bubble'), (1.0, 'dew')):
        try:
            s.update(CP.PQ_INPUTS, p, q)
            out[lab] = s.T()
        except Exception:
            out[lab] = float('nan')
    return out


def mass_to_mole(y_mass, mw={'O2': 31.9988e-3, 'C3H8': 44.0956e-3}):
    w_o, w_c = y_mass, 1.0 - y_mass
    n_o, n_c = w_o / mw['O2'], w_c / mw['C3H8']
    return np.array([n_o / (n_o + n_c), n_c / (n_o + n_c)])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--json', default='work/a406_rootcause/pcsaft_teqp.json')
    ap.add_argument('--kij', type=float, default=0.0)
    args = ap.parse_args()

    print('== 1. literature propane parameters vs CoolProp ==')
    p = LIT['Propane']
    model_c3 = make_pure(p['m'], p['sigma'], p['eps'])
    for t in (230.0, 250.0, 300.0, 350.0):
        st = AbstractState('HEOS', 'Propane'); st.update(CP.QT_INPUTS, 0.0, t)
        rl_cp = st.rhomass() / st.molar_mass()
        out = model_c3.pure_VLE_T(t, rl_cp, 100.0, 100, None)
        print(f'  T={t:5.0f} K  rho_L: PC-SAFT {out[0]:8.1f} vs CoolProp {rl_cp:8.1f} mol/m3 '
              f'({100*(out[0]-rl_cp)/rl_cp:+.2f}%)')

    print('\n== 2. fit O2 (m, sigma, eps/k) to CoolProp saturation densities ==')
    o2 = fit_o2()
    print(f"  O2: m={o2['m']:.4f} sigma={o2['sigma']:.4f} A eps/k={o2['eps']:.2f} K "
          f"(max rho_L err {100*o2['max_rhol_err']:.2f}%, max rho_V err {100*o2['max_rhov_err']:.2f}%)")

    print(f"\n== 3. mixture isobar (kij={args.kij}) vs CoolProp ideal-mixing HEOS ==")
    model = make_mixture(o2, LIT['Propane'], args.kij if args.kij else None)
    print(' Y_O2,mass  x_O2,mol   p[MPa]   T_bub_PCSAFT  T_bub_CoolProp   dT    T_dew_PCSAFT  T_dew_CoolProp   dT')
    rows = []
    for y in (0.7273, 0.8728, 0.93):
        x = mass_to_mole(y)
        for pp in (2.0e6, 3.0e6, 4.0e6):
            cp = coolprop_bulb(pp, y)
            got = {}
            for lab, xspec in (('bubble', x), ('dew', x)):
                try:
                    r = model.mixture_VLE_px(pp, xspec, cp.get(lab, 250.0) or 250.0,
                                             np.array([x[0] * 5000.0, x[1] * 5000.0]),
                                             np.array([x[0] * 100.0, x[1] * 100.0]))
                    t_pcsaft = float(r.T) if hasattr(r, 'T') else float(r[0])
                except Exception:
                    t_pcsaft = float('nan')
                got[lab] = t_pcsaft
            db = got['bubble'] - cp['bubble'] if np.isfinite(got['bubble']) and np.isfinite(cp['bubble']) else float('nan')
            dd = got['dew'] - cp['dew'] if np.isfinite(got['dew']) and np.isfinite(cp['dew']) else float('nan')
            print(f' {y:9.4f} {x[0]:9.4f} {pp/1e6:8.1f} {got["bubble"]:12.2f} {cp["bubble"]:14.2f} {db:+7.2f} '
                  f'{got["dew"]:12.2f} {cp["dew"]:14.2f} {dd:+7.2f}')
            rows.append(dict(y=y, x=float(x[0]), p=pp, t_bub_pcsaft=got['bubble'],
                             t_bub_cp=cp['bubble'], t_dew_pcsaft=got['dew'], t_dew_cp=cp['dew']))
    json.dump(dict(o2=o2, propane=LIT['Propane'], kij=args.kij, rows=rows),
              open(args.json, 'w'), indent=1)
    print('wrote', args.json)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
