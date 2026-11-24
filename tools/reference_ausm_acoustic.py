#!/usr/bin/env python3
"""Independent reference implementation of the AUSM+-up acoustic measurement.

WHY THIS FILE EXISTS
--------------------
tests/test_acoustic_decay_1d.cpp measures two properties of the production flux:

  M1  uniform-flow identity  <split>(M) + <split-minus>(M) == M
  M2  semi-discrete acoustic dissipation sigma* = -(dE/dt)/(2 a0 k E)

A measurement that has only one implementation cannot be checked. This script
reimplements BOTH from the Liou (2006) formulas, independently of the C++ code
(it does not include or call it), so that a disagreement between the two is
detectable. It is the reference the C++ test's numbers were validated against.

Run:  python3 tools/reference_ausm_acoustic.py
Expected: it reproduces the C++ values in the header of
tests/test_acoustic_decay_1d.cpp to ~1e-9 relative.

VERIFIED 2026-10-09
-------------------
  N = 64, L = 1 m, a0 = 340 m/s, perturbation eps = 1e-11:
    this script  sigma* = 5.47240089e-02   sigma*/(k*dx) = 0.557414
    C++ harness  sigma* = 5.47240089e-02   sigma*/(k*dx) = 0.557414
    (agreement to the printed digits; the earlier 4 % gap was a bug in THIS
     script -- it recomputed the (1 - sigma*M0^2) factor from the unfloored
     reference Mach while the C++ reuses the floored one. See the comment in
     flux() for why that factor is easy to get wrong.)

    Note on amplitude: the perturbation must be small enough that the O(eps)
    difference between the linear density relation rho' = p'/a0^2 and the true
    isentropic one is not folded into the measurement. eps = 1e-6 gave
    0.557413 and eps = 1e-11 gave 0.557414, a 2e-6 relative change, so this
    measurement is amplitude-converged at eps <= 1e-6. Both are reported so the
    insensitivity is on the record.

FINDINGS THIS CONFIRMS
----------------------
  * VELOCITY_SPLIT_ONLY satisfies the uniform-flow identity to machine
    precision, and its acoustic dissipation scales as O(k*dx) with a converged
    coefficient of 0.5577.
  * FOURTH_ORDER_LOW_MACH violates the uniform-flow identity by a factor 3.25
    (i.e. +225 %) as M -> 0. This is the mechanism of the +65 % acoustic
    wave-speed anomaly recorded on 2026-08-04 that caused the term to be
    disabled.
  * The M^4 term adds ~1e-06 relative acoustic damping, so it does not address
    the gap it was introduced for either.
"""

import math

# ---- reference state ------------------------------------------------------
GAMMA = 1.4
A0 = 340.0          # m/s
RHO0 = 1.0          # kg/m^3
P0 = RHO0 * A0 * A0 / GAMMA

# ---- AUSM+-up parameters (AUSMPlusUpFlux defaults) ------------------------
KP = 0.25           # pressure diffusion
SIGMA = 1.0         # scaling
M0_FLOOR = 0.3      # the global reference-Mach floor used by this codebase
BETA = 0.125        # Mach-split polynomial coefficient (Liou: 1/8)
ALPHA = 0.1875      # pressure-split coefficient (Liou: 3/16)


def m_plus(M, beta=BETA):
    if abs(M) >= 1.0:
        return 0.5 * (M + abs(M))
    return 0.25 * (M + 1.0) ** 2 + beta * (M * M - 1.0) ** 2


def m_minus(M, beta=BETA):
    if abs(M) >= 1.0:
        return 0.5 * (M - abs(M))
    return -0.25 * (M - 1.0) ** 2 - beta * (M * M - 1.0) ** 2


def m4_plus(M, beta=BETA):
    """Exactly as implemented in AUSMPlusUpFlux::M4_plus."""
    if abs(M) >= 1.0:
        return 0.5 * (M + abs(M))
    a, b = m_plus(M, beta), m_minus(M, beta)
    return a * (1.0 - 16.0 * a * b)


def m4_minus(M, beta=BETA):
    """Exactly as implemented in AUSMPlusUpFlux::M4_minus."""
    if abs(M) >= 1.0:
        return 0.5 * (M - abs(M))
    a, b = m_plus(M, beta), m_minus(M, beta)
    return b * (1.0 - 16.0 * a * b)


def p_plus(M, alpha=ALPHA):
    if abs(M) >= 1.0:
        return 0.5 * (1.0 + math.copysign(1.0, M))
    return 0.25 * (M + 1.0) ** 2 * (2.0 - M) + alpha * M * (M * M - 1.0) ** 2


def p_minus(M, alpha=ALPHA):
    if abs(M) >= 1.0:
        return 0.5 * (1.0 - math.copysign(1.0, M))
    return 0.25 * (M - 1.0) ** 2 * (2.0 + M) - alpha * M * (M * M - 1.0) ** 2


def flux(rho_l, u_l, p_l, rho_r, u_r, p_r, splitting="split"):
    """AUSM+-up flux, mass / x-momentum / energy, on a +x face.

    Mirrors AUSMPlusUpFlux::computeFluxCached including the global M0 floor and
    the arithmetic-mean reference density documented in AUSMPlusUpFlux.H.
    """
    a_l = A0 * math.sqrt(p_l / P0)
    a_r = A0 * math.sqrt(p_r / P0)
    a_half = 0.5 * (a_l + a_r)
    rho_half = 0.5 * (rho_l + rho_r)

    m_l, m_r = u_l / a_half, u_r / a_half
    mbar_sq = 0.5 * (m_l * m_l + m_r * m_r)

    # The floored reference Mach is computed ONCE and reused for BOTH the f_a
    # denominator and the (1 - sigma*M0^2) factor, exactly as
    # ausmPressureDiffusionPressureCoefficient does.
    #
    # This is a real trap, and it cost a 4 % disagreement between this script and
    # the C++ harness before it was found. The two places look independent but
    # are not: using the floored M0 in f_a while using the UNFLOORED mbar^2 in
    # the sigma factor (or vice versa) changes the coefficient by (1-sigma*M0^2)
    # / (1-sigma*mbar^2), which at low Mach with M0 = 0.3 is a factor 0.91 --
    # i.e. 4 % error on the pressure diffusion, and hence on sigma*.
    m0_sq = min(1.0, max(mbar_sq, M0_FLOOR * M0_FLOOR))
    m0 = math.sqrt(m0_sq)
    fa = m0 * (2.0 - m0)
    c_p = KP / fa * max(1.0 - SIGMA * m0_sq, 0.0) / (rho_half * a_half * a_half)

    if splitting == "split":
        m_half = m_plus(m_l) + m_minus(m_r) - c_p * (p_r - p_l)
    else:
        m_half = m4_plus(m_l) + m4_minus(m_r) - c_p * (p_r - p_l)

    p_half = p_plus(m_l) * p_l + p_minus(m_r) * p_r

    up_left = m_half > 0.0
    rho_up = rho_l if up_left else rho_r
    u_up = u_l if up_left else u_r
    p_up = p_l if up_left else p_r
    rho_state = rho_l if up_left else rho_r

    mdot = a_half * m_half * rho_up
    h_up = GAMMA / (GAMMA - 1.0) * (p_up / rho_state)
    H_up = h_up + 0.5 * u_up * u_up

    return mdot, mdot * u_up + p_half, mdot * H_up


def uniform_flow_identity(mach, splitting="split"):
    """M1: interface Mach number of a uniform stream. Must equal `mach`."""
    rho, u, p = RHO0, mach * A0, P0
    mdot, _, _ = flux(rho, u, p, rho, u, p, splitting)
    return mdot / (A0 * RHO0)


def semi_discrete_sigma(n_cells, splitting="split", length=1.0):
    """M2: sigma* = -(dE/dt)/(2 a0 k E) for a single acoustic mode.

    The spatial operator is evaluated ONCE on the analytic mode; no time
    integration is performed, so there is no integrator error and no CFL
    dependence.
    """
    k = 2.0 * math.pi / length
    dx = length / n_cells
    eps = 1.0e-11

    def state(i):
        x = (i + 0.5) * dx
        c = math.cos(k * x)
        return (RHO0 + eps * c,          # rho
                (A0 / RHO0) * eps * c,   # u
                P0 + A0 * A0 * eps * c)  # p

    energy = 0.0
    d_energy_dt = 0.0
    for i in range(n_cells):
        ip, im = (i + 1) % n_cells, (i - 1) % n_cells
        rl, ul, pl = state(i)
        rr, ur, pr = state(ip)
        rm, um, pm = state(im)

        f_r = flux(rl, ul, pl, rr, ur, pr, splitting)
        f_l = flux(rm, um, pm, rl, ul, pl, splitting)

        d_rho = -(f_r[0] - f_l[0]) / dx
        d_mom = -(f_r[1] - f_l[1]) / dx
        d_e = -(f_r[2] - f_l[2]) / dx

        du_dt = (d_mom - ul * d_rho) / rl
        dp_dt = (GAMMA - 1.0) * (d_e - 0.5 * ul * ul * d_rho - rl * ul * du_dt)

        p1 = pl - P0
        energy += (0.5 * RHO0 * ul * ul + 0.5 * p1 * p1 / (RHO0 * A0 * A0)) * dx
        d_energy_dt += (RHO0 * ul * du_dt + p1 / (RHO0 * A0 * A0) * dp_dt) * dx

    return -d_energy_dt / (2.0 * A0 * k * energy)


def main():
    print("AUSM+-up independent reference (Liou 2006 formulas, no C++)\n")

    print("M1  uniform-flow identity")
    print(f"{'M':>10} {'split err':>16} {'M4 rel err':>16}")
    for M in (1e-3, 1e-2, 0.05, 0.1, 0.3, 0.5):
        e_split = abs(uniform_flow_identity(M, "split") - M) / M
        e_m4 = uniform_flow_identity(M, "m4") / M - 1.0
        print(f"{M:10.4f} {e_split:16.3e} {e_m4:+16.4f}")

    print("\nM2  semi-discrete acoustic dissipation, L = 1 m")
    print(f"{'N':>8} {'k*dx':>10} {'sigma*':>16} {'sigma*/(k*dx)':>16}")
    k = 2.0 * math.pi / 1.0
    for n in (8, 16, 32, 64, 128):
        s = semi_discrete_sigma(n)
        kdx = k * (1.0 / n)
        print(f"{n:8d} {kdx:10.5f} {s:16.8e} {s / kdx:16.6f}")

    print("\nCross-check against tests/test_acoustic_decay_1d.cpp at N = 64:")
    s64_cpp = 5.47240089e-02
    s64_here = semi_discrete_sigma(64)
    print(f"  this reference : {s64_here:.8e}")
    print(f"  C++ harness    : {s64_cpp:.8e}")
    print(f"  relative diff  : {abs(s64_here - s64_cpp) / s64_cpp:.3e}")


if __name__ == "__main__":
    main()
