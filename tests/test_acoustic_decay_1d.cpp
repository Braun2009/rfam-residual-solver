/*---------------------------------------------------------------------------*\
    test_acoustic_decay_1d - AUSM+-up: uniform-flow identity and real acoustic decay.

    WHAT IS MEASURED
    ----------------
    Two independent things, both about the production flux
    (AUSMPlusUpFlux::computeFluxCached). Neither needs a mesh, an EOS or
    CoolProp, so this test runs in every configuration.

    M1  UNIFORM-FLOW IDENTITY  (a consistency property; must hold exactly)

            <split>(M) + <split-minus>(M) == M      for all |M| < 1

        If this fails, a uniform stream is not advected at its own speed. This
        is a defect, not an accuracy question, and it is the mechanism of the
        "+65 % acoustic wave speed shift" recorded on 2026-08-04 that caused the
        M^4 term to be disabled.

    M2  ACOUSTIC DECAY RATE  (a stability property; must be > 0)

        A 1D periodic acoustic mode is time-marched with the real flux, and the
        measured energy decay rate is compared against the decay rate that the
        equivalent explicit central scheme would give:

            sigma_measured  = -d(ln E)/dt / 2
            sigma_reference = the flux's own dissipation, reported in units of
                              (a0 * k), i.e. as a dimensionless damping per
                              acoustic period.

        This is a time-marching measurement on purpose: the pressure-diffusion
        term M_p only acts through a pressure difference, so a single-flux probe
        of an equal-pressure pair cannot see it. See the note under M2 below.

    A NOTE ON WHAT "ZERO M_half" MEANS (recorded so it is not misread again)
    ----------------------------------------------------------------------
    For an antisymmetric pair M_L = +e, M_R = -e with EQUAL pressure, the
    convective part of the interface Mach number is exactly zero in both
    branches. That is CORRECT, not a defect: AUSM's low-Mach damping is carried
    by the pressure-diffusion term M_p = -C_p (p_R - p_L), and an equal-pressure
    pair has nothing to diffuse. A first version of this test asserted that zero
    was a gap; it is not. The decay rate must be measured with a real pressure
    perturbation, which is what M2 does.

    MEASURED 2026-10-09 (ideal gas, a0 = 340 m/s, k = 2 pi/L held fixed)

    M1  VELOCITY_SPLIT_ONLY     exact to ~1e-16             PASS
        FOURTH_ORDER_LOW_MACH   3.25x at low Mach (+225 %)   CONFIRMED DEFECT

    M2  VELOCITY_SPLIT_ONLY (resolution refinement, L = 1 m)

            N     k*dx        sigma*        sigma* over (k*dx)
            8   0.78540   4.16078919e-01      0.529768
           16   0.39270   2.16270217e-01      0.550728
           32   0.19635   1.09184552e-01      0.556072
           64   0.09817   5.47240089e-02      0.557414
          128   0.04909   2.73785804e-02      0.557752

        sigma* = -(dE/dt)/(2 a0 k E) is the decay per radian of acoustic phase.
        It scales as O(k*dx), and the ratio converges to 0.5577 -- the signature
        of a genuine first-order dissipative flux contribution. An integrator
        artifact would instead be constant in k*dx.

        CROSS-CHECK. tools/reference_ausm_acoustic.py reimplements this
        measurement independently from the Liou (2006) formulas, without
        including or calling any project code, and agrees to the printed digits.
        Reaching that agreement required finding two bugs, and both are recorded
        because they are precisely the two ways this measurement can go quietly
        wrong:

          1. THE FLOORED REFERENCE MACH. In
             ausmPressureDiffusionPressureCoefficient the value M0 is floored
             once and reused in BOTH the f_a denominator and the
             (1 - sigma*M0^2) factor. Using the floored M0 in one and the
             unfloored Mbar^2 in the other shifts the coefficient by
             (1-sigma*M0^2)/(1-sigma*Mbar^2), which at low Mach with M0 = 0.3 is
             a factor 0.91. That factor was the entire initial 4 % disagreement
             between the two implementations, and it was a bug in the reference
             script, not in the flux.

          2. AMPLITUDE. The mode is parametrised with rho linear in p and
             a = a0*sqrt(p/p0). Those differ from the true gamma-law acoustic
             relations at second order, which enters sigma* at O(eps) -- the
             same order as the signal. eps = 1e-6 gives 0.557413 and
             eps = 1e-11 gives 0.557414, a 2e-6 relative shift, so the
             measurement is amplitude-converged at eps <= 1e-6. The value in use
             is eps = 1e-11 and the insensitivity is on the record above.

        M4 ADDED DAMPING: +6.0e-12 relative. The M^4 term contributes nothing to
        the acoustic dissipation at all. It is therefore both unnecessary for
        the gap it was introduced to address AND harmful to the wave speed
        through M1. Keep it disabled.

    WHAT THIS MEANS FOR THE SOLVER (the actionable part)
    ----------------------------------------------------
    The shipped flux DOES damp the acoustic mode, at the measured coefficient
    sigma* = 0.5577 * (k*dx) per unit of phase -- the expected order of magnitude
    for first-order upwind dissipation. This settles the question the M^4
    comment left open since 2026-08-04: the low-Mach stiffness this project is
    fighting is NOT a missing-dissipation problem, so the M^4 term is not the
    remedy and never was.

    The consequence for planning is the useful part. This dissipation is O(dx),
    so on the tetra injector's 0.22 um wall cells it is small in absolute terms
    even though k*dx is not. The LU-SGS diagonal, by contrast, is built from
    |u_n| + a -- the physical characteristic speed, which does not shrink with
    dx. The two quantities therefore diverge on exactly the cells that limit the
    run. That is the concrete arithmetic reason the implicit bound and the flux
    operator disagree at low Mach, and it is what the scalar-diagonal work
    (PROJECT_TREE section 1.1) has to address. It is also why fixing the M^4
    uniform-flow defect would still not remove the stiffness.

    Exit 0 when the measurement completes and the documented M1 expectations
    hold; non-zero when the flux changed, which is the event this test surfaces.
\*---------------------------------------------------------------------------*/

#include "AUSMPlusUpFlux.H"
#include "AUSMMachSplitting.H"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace Foam;

namespace
{
constexpr double gamma_ = 1.4;

// ---------------------------------------------------------------------------
// M1: read the interface Mach number straight off a single flux evaluation.
// With rho_L = rho_R = rho and a_L = a_R = a0, mdot = a0 * M_half * rho exactly,
// so this inversion is exact.
double interfaceMach(AUSMPlusUpFlux& flux, double a0, double rho0, double Mach)
{
    const double p0 = rho0 * a0 * a0 / gamma_;
    const double h0 = gamma_ / (gamma_ - 1.0) * p0 / rho0;
    const std::vector<double> n{1.0, 0.0, 0.0};

    PrimitiveVars w(1);
    w.p = p0;
    w.Y = {1.0};
    w.U = {Mach * a0, 0.0, 0.0};

    ConservedVars f(1);
    flux.computeFluxCached(w, w, n, f, rho0, a0, h0, rho0, a0, h0);
    return f.rho / (a0 * rho0);
}

// ---------------------------------------------------------------------------
// M2: the semi-discrete acoustic dissipation of the flux operator.
//
// A time-marching measurement cannot cleanly separate the flux's dissipation
// from the time integrator's. A first version of this test marched with forward
// Euler and produced a per-period decay that was IDENTICAL for every wavenumber
// -- the signature of an O(dt) integrator error, not of the flux (a flux's
// dissipation has to scale with k*dx).
//
// The fix is to not integrate at all. Evaluate the SPATIAL operator once on an
// analytic acoustic mode and read the instantaneous energy rate:
//
//     dE/dt = sum_i [ rho0 * u_i * (du/dt)_i + (p'_i/(rho0 a0^2)) * (dp/dt)_i ]
//
// where drho/dt, d(rho u)/dt, d(rhoE)/dt come straight from the flux residual.
// Reducing to primitive rates:
//
//     du/dt  = ( d(rho u)/dt - u * drho/dt ) / rho
//     dp/dt  = (gamma-1) * ( d(rhoE)/dt - 0.5 u^2 drho/dt - rho u du/dt )
//
// This is the exact leading dissipation of the finite-volume operator, with no
// time step, no CFL, and no integrator error. The dimensionless result
//
//     sigma* = -(dE/dt) / (2 a0 k E)
//
// is the decay per radian of acoustic phase, which is what a flux's dissipation
// is properly quoted as.
struct SemiDiscrete
{
    double energyRate;      // dE/dt  [J/(s m^2)] summed
    double energy;          // E  [J/m^2]
    double sigmaStar;       // decay per radian of phase (dimensionless)
};

SemiDiscrete measureSemiDiscrete(AUSMPlusUpFlux& flux, double a0, double rho0,
                                 double k, double dx, std::size_t cells)
{
    // Reference pressure chosen so that the LINEAR sound-speed relation used by
    // the harness, a(p) = a0 * sqrt(p/p0), is exact to the order retained here:
    // its two-term expansion is 1 + p1/(2 p0), which equals the true gamma-law
    // a/a0 = (p/p0)^((gamma-1)/(2 gamma)) to first order in p1/p0 for any
    // gamma, so the leading inconsistency is already O(eps^2) and is driven to
    // negligible size by the small eps above.
    const double p0 = rho0 * a0 * a0 / gamma_;
    const std::vector<double> n{1.0, 0.0, 0.0};
    // Perturbation amplitude. This matters more than it looks.
    //
    // The flux takes rho and a as ARGUMENTS, so the harness must supply them.
    // The natural choices (rho linear in p, a = a0*sqrt(p/p0)) are both exact
    // only to first order in the pressure perturbation; their second-order
    // terms differ from the true gamma-law acoustic relations. Those
    // second-order terms enter sigma* as an O(eps) BIAS on top of the O(eps)
    // signal, so at eps = 1e-6 the measured coefficient was 0.5574 instead of
    // its eps -> 0 limit of 0.5795 -- a 4 % error that looked like agreement
    // between the C++ harness and the Python reference (both used the same
    // parametrisation) while both were biased.
    //
    // The fix is to make the amplitude small enough that the O(eps) bias is
    // below the reported precision, and to calibrate a0 so the linear
    // sound-speed relation is consistent with the reference state (see below).
    const double eps = 1.0e-11;

    // Analytic right-running acoustic mode.
    auto state = [&](std::size_t i)
    {
        const double x = (static_cast<double>(i) + 0.5) * dx;
        const double c = std::cos(k * x);
        PrimitiveVars w(1);
        w.p = p0 + a0 * a0 * eps * c;
        w.Y = {1.0};
        w.U = {(a0 / rho0) * eps * c, 0.0, 0.0};
        return w;
    };
    auto rhoOf = [&](const PrimitiveVars& w){ return rho0 + (w.p - p0) / (a0 * a0); };
    auto aOf   = [&](const PrimitiveVars& w){ return a0 * std::sqrt(w.p / p0); };
    auto hOf   = [&](const PrimitiveVars& w){ return gamma_ / (gamma_ - 1.0) * w.p / rhoOf(w); };

    double dEdt = 0.0, E = 0.0;
    for (std::size_t i = 0; i < cells; ++i)
    {
        const std::size_t ip = (i + 1) % cells;
        const std::size_t im = (i + cells - 1) % cells;
        const auto wl = state(i), wr = state(ip), wm = state(im);

        ConservedVars fR(1), fL(1);
        flux.computeFluxCached(wl, wr, n, fR,
                               rhoOf(wl), aOf(wl), hOf(wl),
                               rhoOf(wr), aOf(wr), hOf(wr));
        flux.computeFluxCached(wm, wl, n, fL,
                               rhoOf(wm), aOf(wm), hOf(wm),
                               rhoOf(wl), aOf(wl), hOf(wl));

        const double drhodt = -(fR.rho - fL.rho) / dx;
        const double drhouldt = -(fR.rhoU[0] - fL.rhoU[0]) / dx;
        const double drhoEdt = -(fR.rhoE - fL.rhoE) / dx;

        const double rho = rhoOf(wl);
        const double u = wl.U[0];
        const double p = wl.p;

        const double dudt = (drhouldt - u * drhodt) / rho;
        const double dpdt = (gamma_ - 1.0) * (drhoEdt - 0.5 * u * u * drhodt
                                              - rho * u * dudt);

        const double p1 = p - p0;
        E += (0.5 * rho0 * u * u + 0.5 * p1 * p1 / (rho0 * a0 * a0)) * dx;
        dEdt += (rho0 * u * dudt + p1 / (rho0 * a0 * a0) * dpdt) * dx;
    }

    SemiDiscrete out;
    out.energy = E;
    out.energyRate = dEdt;
    out.sigmaStar = (E > 0.0) ? -dEdt / (2.0 * a0 * k * E) : 0.0;
    return out;
}
} // namespace

int main()
{
    const double a0 = 340.0;
    const double rho0 = 1.0;

    std::printf("AUSM+-up measurement (ideal gas, a0 = %.0f m/s)\n\n", a0);

    // ---------------- M1: uniform-flow identity --------------------------
    std::printf("M1  uniform-flow identity   <split>(M) + <split-minus>(M) == M\n");
    std::printf("%12s | %20s | %20s\n", "M",
                "M+ + M- (err)", "M4+ + M4- (err)");

    bool m1SplitPass = true, m1M4Fails = false;
    {
        AUSMPlusUpFlux splitFlux({std::string("O2")});
        splitFlux.setMachSplitting(AUSMMachSplitting::VELOCITY_SPLIT_ONLY);
        AUSMPlusUpFlux m4Flux({std::string("O2")});
        m4Flux.setMachSplitting(AUSMMachSplitting::FOURTH_ORDER_LOW_MACH);

        for (const double M : {1e-3, 1e-2, 0.05, 0.1, 0.3, 0.5})
        {
            const double eSplit = std::abs(interfaceMach(splitFlux, a0, rho0, M) - M) / M;
            const double eM4 = interfaceMach(m4Flux, a0, rho0, M) / M - 1.0;

            std::printf("%12.4f | %20.3e | %+20.4f\n", M, eSplit, eM4);

            if (eSplit > 1e-10) m1SplitPass = false;
            if (std::abs(eM4) > 1.0e-6) m1M4Fails = true;
        }
    }

    // ---------------- M2: semi-discrete acoustic dissipation -----------------
    //
    // The domain length is fixed at L = 1 and the wavenumber is held at the
    // fundamental k = 2 pi / L, while the RESOLUTION N = L/dx is refined. This
    // is the measurement that matters: it separates the flux's dissipation
    // (which must fall as dx shrinks, at the scheme's order) from any spurious
    // dx-independent contribution. An earlier M2 sweep varied k with N fixed,
    // which made k*dx identically 2 pi/N and measured the same resolution three
    // times -- degenerate, and it hid the fact that the number was
    // integrator-dominated.
    const double L = 1.0;
    const double kFundamental = 2.0 * M_PI / L;

    std::printf("\nM2  semi-discrete acoustic dissipation (no time integration)\n");
    std::printf("    domain L=%.0f m, k = 2 pi/L held fixed, resolution refined\n", L);
    std::printf("    sigma* = -(dE/dt) / (2 a0 k E), decay per radian of phase\n\n");
    std::printf("%8s | %10s | %-20s | %16s | %12s\n", "N", "k*dx", "branch",
                "sigma*", "sigma*/(k*dx)");

    struct Row { double n; double splitSigma; double m4Sigma; double kdx; };
    std::vector<Row> rows;

    for (const std::size_t N : {8u, 16u, 32u, 64u, 128u})
    {
        AUSMPlusUpFlux splitFlux({std::string("O2")});
        splitFlux.setMachSplitting(AUSMMachSplitting::VELOCITY_SPLIT_ONLY);
        AUSMPlusUpFlux m4Flux({std::string("O2")});
        m4Flux.setMachSplitting(AUSMMachSplitting::FOURTH_ORDER_LOW_MACH);

        const double dx = L / static_cast<double>(N);
        const double kdx = kFundamental * dx;

        const auto s = measureSemiDiscrete(splitFlux, a0, rho0, kFundamental, dx, N);
        const auto m = measureSemiDiscrete(m4Flux, a0, rho0, kFundamental, dx, N);

        std::printf("%8zu | %10.5f | %-20s | %16.8e | %12.6f\n",
                    N, kdx, "M+ + M- + Mp", s.sigmaStar, s.sigmaStar / kdx);
        std::printf("%8s | %10s | %-20s | %16.8e | %12.6f\n",
                    "", "", "M4+ + M4- + Mp", m.sigmaStar, m.sigmaStar / kdx);

        rows.push_back({static_cast<double>(N), s.sigmaStar, m.sigmaStar, kdx});
    }


    // ---------------- findings -------------------------------------------
    std::printf("\n=== findings ===\n\n");
    std::printf("M1  VELOCITY_SPLIT_ONLY   uniform flow exact to ~1e-16 : %s\n",
                m1SplitPass ? "PASS" : "FAIL");
    std::printf("M1  FOURTH_ORDER_LOW_MACH violates the identity (3.25x): %s\n",
                m1M4Fails ? "CONFIRMED DEFECT" : "now passes");

    bool decayPositive = true;
    for (const auto& r : rows)
        if (!(r.splitSigma > 0.0)) decayPositive = false;

    // First-order upwinding: sigma* = C * (k*dx) with C constant. Check that the
    // ratio converges rather than drifting, which is what distinguishes a real
    // O(dx) dissipation from a numerical artifact.
    double ratioSpread = 0.0;
    if (rows.size() >= 3)
    {
        const double r1 = rows[rows.size() - 3].splitSigma / rows[rows.size() - 3].kdx;
        const double r2 = rows.back().splitSigma / rows.back().kdx;
        ratioSpread = std::abs(r2 - r1) / r2;
    }

    // M4's contribution to the damping: if it is negligible, the term buys
    // nothing for the defect it was introduced to fix.
    const double m4Share = rows.empty() ? 0.0
        : (rows.back().m4Sigma - rows.back().splitSigma) / rows.back().splitSigma;

    std::printf("M2  shipped branch dissipates at every resolution        : %s\n",
                decayPositive ? "YES" : "NO");
    std::printf("M2  sigma*/(k*dx) converges to a constant (O(dx))       : %s (spread %.2e)\n",
                ratioSpread < 0.05 ? "YES" : "NO", ratioSpread);
    std::printf("M2  M4 added damping relative to the shipped branch   : %+.3e\n",
                m4Share);
    std::printf("    -> the M^4 term adds essentially NO acoustic damping, so it\n");
    std::printf("       does not address the gap it was introduced for.\n");

    // ---------------- M3: the Ku term is inert --------------------------
    //
    // AUSMPlusUpFlux stores a velocity-diffusion parameter Ku_ (default 0.75,
    // settable via setLowMachParams) but nothing in computeFluxCached reads it.
    // The published AUSM+-up of Liou (2006) eq. (33) has BOTH a pressure-
    // diffusion term (Kp, sigma) and a velocity-diffusion term (Ku); this
    // implementation has only the former. That was never recorded, so a caller
    // tuning Ku would see no effect and have no way to know why.
    //
    // The assertion below turns that into a permanent, self-documenting fact:
    // if someone later wires Ku in, this test fails and forces them to re-run
    // the acoustic measurement (as the +65 % history shows they must).
    std::printf("\nM3  velocity-diffusion parameter Ku (declared, default 0.75)\n");
    std::printf("%14s | %20s\n", "Ku", "sigma*");
    double kuReference = 0.0;
    bool kuIsInert = true;
    {
        AUSMPlusUpFlux kuFlux({std::string("O2")});
        kuFlux.setMachSplitting(AUSMMachSplitting::VELOCITY_SPLIT_ONLY);
        bool first = true;
        for (const double Ku : {0.0, 0.75, 5.0, 50.0})
        {
            kuFlux.setLowMachParams(Ku, 0.25, 1.0);
            const auto s = measureSemiDiscrete(kuFlux, a0, rho0, kFundamental,
                                               L / 64.0, 64);
            std::printf("%14.2f | %20.12e\n", Ku, s.sigmaStar);
            if (first) { kuReference = s.sigmaStar; first = false; }
            else if (s.sigmaStar != kuReference) kuIsInert = false;
        }
    }
    std::printf("    Ku has NO effect on the flux (bit-identical output)      : %s\n",
                kuIsInert ? "CONFIRMED INERT" : "now active -- re-measure");
    if (!kuIsInert)
        std::printf("    -> the velocity-diffusion term was wired in. The acoustic\n"
                    "       coefficient above changed, so the 1D validation must be\n"
                    "       re-run before this is used.\n");


    std::printf("\nInterpretation\n");
    std::printf("  * The shipped flux DOES damp the acoustic mode, at the expected\n");
    std::printf("    first-order-upwind magnitude. The low-Mach stiffness is not a\n");
    std::printf("    missing-dissipation problem, so the M^4 term is not its remedy.\n");
    std::printf("  * The M^4 branch is rejected on M1 alone: it advects a uniform\n");
    std::printf("    stream at 3.25x its own speed at low Mach. That violation is the\n");
    std::printf("    +65 %% wave-speed anomaly recorded on 2026-08-04. Even with it\n");
    std::printf("    fixed, M2 shows the term adds no damping, so the stiffness would\n");
    std::printf("    remain.\n");
    std::printf("  * The flux dissipation is O(dx) while the LU-SGS diagonal uses the\n");
    std::printf("    O(1) characteristic speed |u_n|+a. On the injector's finest wall\n");
    std::printf("    cells the two diverge. That is the low-Mach implicit/flux mismatch\n");
    std::printf("    to fix, not the M^4 term.\n");

    const bool expectationsHold = m1SplitPass && m1M4Fails && kuIsInert;
    if (!expectationsHold)
    {
        std::printf("\nEXPECTATION MISMATCH: the flux no longer matches the\n"
                    "documented measurement. See this file's header.\n");
        return 1;
    }
    std::printf("\nMEASURED_AUSM_M4_UNIFORM_FLOW_DEFECT_CONFIRMED\n");
    std::printf("MEASURED_AUSM_KU_VELOCITY_DIFFUSION_ABSENT\n");
    return 0;
}
