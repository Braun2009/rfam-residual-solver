#!/usr/bin/env python3
"""Predictive PC-SAFT model for the O2 + C3H8 mixture (uses no mixture data).

LOX + hydrocarbon mixtures are auto-oxidation/explosion systems: nobody measures
their VLE or density, and no public library carries a binary departure function
for the pair (CoolProp/GERG treat O2-C3H8 as ideal mixing, F=0). The predictive
route is a molecular EOS whose binary interaction is *defaulted*: PC-SAFT with
kij = 0. Its pure-component parameters are fitted here to CoolProp's pure-fluid
reference data only (safe and abundant), so no hazardous experiment is involved.

Steps: implement PC-SAFT residual Helmholtz -> fit (m, sigma, eps/k) for O2 and
C3H8 to CoolProp saturation data -> validate the pure-fluid fit -> compute the
mixture bubble/dew line at the injector pressure with kij=0 and compare with
CoolProp's ideal-mixing HEOS (the model-uncertainty band for the two-phase zone).
"""
from __future__ import annotations
import argparse
import json
import math

import numpy as np
from scipy.optimize import least_squares, brentq

import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState

A_CONST = np.array([0.9105631445, 0.6361281449, 2.6861347891, -26.547362491,
                    97.759208784, -159.59154087, 91.297774084])
B_CONST = np.array([-0.3084016918, -0.0906148351, 0.4527842806, 0.5962700728,
                    -1.7241829131, -4.1302112531, 13.776631870])
N_A = 6.02214076e23
R_GAS = 8.314462618


def a_res(T, rho_molar, x, m, sigma, eps_k, kij=0.0):
    """Residual Helmholtz energy per molecule / kT. rho_molar in mol/m3."""
    rho_n = rho_molar * N_A * 1e-30          # 1/A^3 (sigma given in Angstrom)
    d = sigma * (1.0 - 0.12 * np.exp(-3.0 * eps_k / T))
    mbar = float(np.dot(x, m))
    zeta = np.array([math.pi / 6.0 * rho_n * float(np.sum(x * m * d ** n)) for n in range(4)])
    z0, z1, z2, z3 = zeta
    eta = z3
    if eta <= 0.0:
        return 0.0
    if eta >= 1.0:
        return float('inf')
    a_hs = (1.0 / z0) * (3.0 * z1 * z2 / (1.0 - z3)
                         + z2 ** 3 / (z3 * (1.0 - z3) ** 2)
                         + (z2 ** 3 / z3 ** 2 - z0) * math.log(max(1.0 - z3, 1e-300)))
    a_hc = mbar * a_hs
    for i in range(len(m)):
        dd = d[i]
        g_ii = (1.0 / (1.0 - z3) + (dd / 2.0) * 3.0 * z2 / (1.0 - z3) ** 2
                + (dd / 2.0) ** 2 * 2.0 * z2 ** 2 / (1.0 - z3) ** 3)
        a_hc -= x[i] * (m[i] - 1.0) * math.log(max(g_ii, 1e-300))
    s = np.sqrt(eps_k[:, None] * eps_k[None, :]) * (1.0 - kij)
    sig = 0.5 * (sigma[:, None] + sigma[None, :])
    mm = m[:, None] * m[None, :]
    xx = x[:, None] * x[None, :]
    m2es3 = float(np.sum(xx * mm * (s / T) * sig ** 3))
    m2e2s3 = float(np.sum(xx * mm * (s / T) ** 2 * sig ** 3))
    i1 = float(np.sum(A_CONST * eta ** np.arange(7)))
    i2 = float(np.sum(B_CONST * eta ** np.arange(7)))
    c1 = 1.0 / (1.0 + mbar * (8.0 * eta - 2.0 * eta ** 2) / (1.0 - eta) ** 4
                + (1.0 - mbar) * (20.0 * eta - 27.0 * eta ** 2 + 12.0 * eta ** 3 - 2.0 * eta ** 4)
                / ((1.0 - eta) ** 2 * (2.0 - eta) ** 2))
    a_disp = -2.0 * math.pi * rho_n * i1 * m2es3 - math.pi * rho_n * mbar * c1 * i2 * m2e2s3
    return a_hc + a_disp


def z_factor(T, rho_molar, x, m, sigma, eps_k, kij=0.0):
    h = 1e-7 * max(rho_molar, 1e-9)
    f = lambda r: a_res(T, r, x, m, sigma, eps_k, kij)
    return 1.0 + rho_molar * (f(rho_molar + h) - f(rho_molar - h)) / (2.0 * h)


def pressure(T, rho_molar, x, m, sigma, eps_k, kij=0.0):
    return rho_molar * R_GAS * T * z_factor(T, rho_molar, x, m, sigma, eps_k, kij)


def rho_from_p(T, p, x, m, sigma, eps_k, kij=0.0, branch='auto', rho_guess=None):
    """Density roots of p(T, rho) = p; returns the stable branch requested."""
    def g(r):
        return pressure(T, r, x, m, sigma, eps_k, kij) - p
    grid = np.concatenate([np.linspace(1.0, 500.0, 200), np.linspace(500.0, 40000.0, 400)])
    roots = []
    prev = g(grid[0])
    for r in grid[1:]:
        cur = g(r)
        if cur == 0.0:
            roots.append(r)
        elif prev * cur < 0.0:
            try:
                roots.append(brentq(g, grid[np.where(grid == r)[0][0] - 1], r, xtol=1e-12))
            except Exception:
                pass
        prev = cur
    if not roots:
        raise RuntimeError('no density root')
    if branch == 'liquid':
        return max(roots)
    if branch == 'vapour':
        return min(roots)
    if rho_guess is None:
        return max(roots)
    return min(roots, key=lambda r: abs(r - rho_guess))


def lnphi(T, p, x, m, sigma, eps_k, kij=0.0, branch='auto', rho_guess=None):
    """ln fugacity coefficients (residual part) via numerical composition derivative."""
    x = np.asarray(x, dtype=float)
    x = x / x.sum()
    rho = rho_from_p(T, p, x, m, sigma, eps_k, kij, branch, rho_guess)
    a = a_res(T, rho, x, m, sigma, eps_k, kij)
    z = z_factor(T, rho, x, m, sigma, eps_k, kij)
    grad = np.zeros(len(x))
    h = 1e-6
    for i in range(len(x)):
        xp = x.copy()
        xp[i] += h
        xp = xp / xp.sum()
        grad[i] = (a_res(T, rho, xp, m, sigma, eps_k, kij) - a) / h
    return a + grad - float(np.dot(x, grad)) + (z - 1.0) - math.log(z), rho, z


def fit_pure(fluid, tmin, tmax, n, x0, bounds):
    st = AbstractState('HEOS', fluid)
    temps, psat, rhol = [], [], []
    for t in np.linspace(tmin, tmax, n):
        try:
            st.update(CP.QT_INPUTS, 0.0, t)
            temps.append(t); psat.append(st.p()); rhol.append(st.rhomass())
        except Exception:
            pass
    temps = np.array(temps); psat = np.array(psat); rhol = np.array(rhol)
    mw = st.molar_mass() * 1e-3   # kg/mol

    def residuals(theta):
        m, sigma, eps = theta
        mv, sv, ev = np.array([m]), np.array([sigma]), np.array([eps])
        out = []
        for t, p, rl in zip(temps, psat, rhol):
            try:
                lnp_l, _, _ = lnphi(t, p, [1.0], mv, sv, ev, branch='liquid')
                lnp_v, _, _ = lnphi(t, p, [1.0], mv, sv, ev, branch='vapour')
                out.append((lnp_l[0] - lnp_v[0]))
                rho_l = rho_from_p(t, p, [1.0], mv, sv, ev, branch='liquid')
                out.append((rho_l * mw - rl) / rl)
            except Exception:
                out.append(1.0)
        return np.array(out)

    sol = least_squares(residuals, x0, bounds=bounds, xtol=1e-12, ftol=1e-12)
    m, sigma, eps = sol.x
    mv, sv, ev = np.array([m]), np.array([sigma]), np.array([eps])
    # quality: saturation pressure reproduced by the equal-fugacity condition
    perr, derr = [], []
    for t, p, rl in zip(temps, psat, rhol):
        try:
            lnp_l, _, _ = lnphi(t, p, [1.0], mv, sv, ev, branch='liquid')
            lnp_v, _, _ = lnphi(t, p, [1.0], mv, sv, ev, branch='vapour')
            # p_sat from the model: p*exp(lnphi_l - lnphi_v)
            perr.append(math.exp(lnp_l[0] - lnp_v[0]) - 1.0)
            rho_l = rho_from_p(t, p, [1.0], mv, sv, ev, branch='liquid')
            derr.append((rho_l * mw - rl) / rl)
        except Exception:
            pass
    return dict(fluid=fluid, m=float(m), sigma=float(sigma), eps_k=float(eps),
                tmin=float(tmin), tmax=float(tmax), n=len(temps),
                max_psat_rel=float(np.max(np.abs(perr))) if perr else float('nan'),
                max_rholiq_rel=float(np.max(np.abs(derr))) if derr else float('nan'))


def load_params(path_or_fits, fluid):
    if isinstance(path_or_fits, dict):
        f = path_or_fits[fluid]
    else:
        f = json.load(open(path_or_fits))[fluid]
    return f['m'], f['sigma'], f['eps_k']


def mixture_x(y_oxygen_mass):
    """Mass fraction -> mole fraction for O2/C3H8."""
    mw = {'O2': 31.9988e-3, 'C3H8': 44.0956e-3}
    w_o, w_c = y_oxygen_mass, 1.0 - y_oxygen_mass
    n_o, n_c = w_o / mw['O2'], w_c / mw['C3H8']
    return np.array([n_o / (n_o + n_c), n_c / (n_o + n_c)])


def mixture_params(fits):
    m = np.array([fits['Oxygen']['m'], fits['Propane']['m']])
    sigma = np.array([fits['Oxygen']['sigma'], fits['Propane']['sigma']])
    eps = np.array([fits['Oxygen']['eps_k'], fits['Propane']['eps_k']])
    return m, sigma, eps


def bubble_T(p, x, m, sigma, eps, t_lo=110.0, t_hi=380.0):
    def f(t):
        try:
            l_l, _, _ = lnphi(t, p, x, m, sigma, eps, branch='liquid')
            l_v, _, _ = lnphi(t, p, x, m, sigma, eps, branch='vapour')
            k = np.exp(l_l - l_v)
            return float(np.sum(x * k) - 1.0)
        except Exception:
            return float('nan')
    grid = np.linspace(t_lo, t_hi, 40)
    vals = [f(t) for t in grid]
    for i in range(len(vals) - 1):
        a, b = vals[i], vals[i + 1]
        if np.isfinite(a) and np.isfinite(b) and a * b < 0:
            return brentq(f, grid[i], grid[i + 1], xtol=1e-6)
    return float('nan')


def dew_T(p, y, m, sigma, eps, t_lo=110.0, t_hi=380.0):
    def f(t):
        try:
            l_l, _, _ = lnphi(t, p, y, m, sigma, eps, branch='liquid')
            l_v, _, _ = lnphi(t, p, y, m, sigma, eps, branch='vapour')
            k = np.exp(l_l - l_v)
            return float(np.sum(y / k) - 1.0)
        except Exception:
            return float('nan')
    grid = np.linspace(t_lo, t_hi, 40)
    vals = [f(t) for t in grid]
    for i in range(len(vals) - 1):
        a, b = vals[i], vals[i + 1]
        if np.isfinite(a) and np.isfinite(b) and a * b < 0:
            return brentq(f, grid[i], grid[i + 1], xtol=1e-6)
    return float('nan')


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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--json', default='work/a406_rootcause/pcsaft_o2_c3h8.json')
    ap.add_argument('--fit-only', action='store_true')
    args = ap.parse_args()
    fits = {}
    fits['Oxygen'] = fit_pure('Oxygen', 95.0, 150.0, 10, [1.2, 3.0, 110.0],
                              ([0.5, 2.0, 60.0], [3.0, 4.5, 250.0]))
    fits['Propane'] = fit_pure('Propane', 230.0, 360.0, 10, [2.0, 3.6, 210.0],
                               ([1.0, 2.5, 120.0], [4.0, 5.5, 400.0]))
    for k, v in fits.items():
        print(f"{k}: m={v['m']:.4f} sigma={v['sigma']:.4f} A eps/k={v['eps_k']:.2f} K  "
              f"({v['n']} states, {v['tmin']:.0f}-{v['tmax']:.0f} K)  "
              f"max|p_sat err|={v['max_psat_rel']*100:.2f}%  max|rho_l err|={v['max_rholiq_rel']*100:.2f}%")
    json.dump(fits, open(args.json, 'w'), indent=1)
    print('wrote', args.json)
    if args.fit_only:
        return 0
    m, sigma, eps = mixture_params(fits)
    print('\nPC-SAFT kij=0 mixture vs CoolProp HEOS (ideal mixing) — bubble/dew:')
    print(' Y_O2,mass  x_O2,mol   p[MPa]  T_bub_PCSAFT  T_bub_CoolProp  dT   T_dew_PCSAFT  T_dew_CoolProp  dT')
    for y in (0.5, 0.727, 0.872, 0.93, 0.95):
        x = mixture_x(y)
        for p in (2.0e6, 3.0e6, 4.0e6):
            tb = bubble_T(p, x, m, sigma, eps)
            td = dew_T(p, x, m, sigma, eps)
            cp = coolprop_bulb(p, y)
            fmt = lambda v: f'{v:12.2f}' if np.isfinite(v) else '         ---'
            db = f'{tb-cp["bubble"]:5.1f}' if np.isfinite(tb) and np.isfinite(cp['bubble']) else '  ---'
            dd = f'{td-cp["dew"]:5.1f}' if np.isfinite(td) and np.isfinite(cp['dew']) else '  ---'
            print(f' {y:9.3f} {x[0]:9.4f} {p/1e6:7.1f} {fmt(tb)} {fmt(cp["bubble"])} {db} '
                  f'{fmt(td)} {fmt(cp["dew"])} {dd}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
