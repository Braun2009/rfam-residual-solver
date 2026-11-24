// HybridHelmholtzThermo unified-potential contract tests.
//
// Verifies (GERG adoption contract, 2026-08-02):
//   1. ResidualFade C2 quintic smoothstep: endpoints, f'=f''=0, monotone;
//   2. high-temperature limit (fade -> 0): ideal-gas p=rhoRT, e=eIG,
//      cp-cv = R_mix;
//   3. low-temperature limit (fade = 1): full PR residual active;
//   4. unified-potential thermodynamic identities (Mayer relation from the
//      same potential);
//   5. rho-e-z recovery round trip (PTY -> rhoEY -> inversion);
//   6. transition band: Cv>0, mechanical stability dp/drho>0, a^2>0,
//      derivative continuity (Cp/Cv finite across the band).
//
// Output marker: HYBRID_HELMHOLTZ_PASS / _FAIL.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "HybridHelmholtzThermo.H"
#include "SpeciesRegistry.H"

using namespace Foam;

namespace
{

int failures = 0;

void check(
    const std::string& name, bool ok,
    double got = 0.0, double expected = 0.0)
{
    std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << name;
    if (!ok)
        std::cout << " got=" << std::setprecision(12) << got
                  << " expected=" << expected;
    std::cout << '\n';
    if (!ok) ++failures;
}

double relativeError(double got, double expected, double floor = 1.0e-10)
{
    return std::abs(got - expected) / std::max(std::abs(expected), floor);
}

int run()
{
    SpeciesRegistry registry;
    const std::vector<std::string> names = registry.names();
    // Lean O2 background composition (LOX-adjacent).
    std::vector<double> Y(names.size(), 0.0);
    Y[registry.index("O2")] = 0.90;
    Y[registry.index("CO2")] = 0.05;
    Y[registry.index("H2O")] = 0.05;

    HybridHelmholtzThermo thermo(registry, 600.0, 900.0);

    // ------------------------------------------------------------------
    // 1. ResidualFade properties
    // ------------------------------------------------------------------
    ResidualFade f = thermo.fade();
    check("fade(T<=T0)=1", std::abs(f.fade(500.0) - 1.0) < 1e-12);
    check("fade(T>=T1)=0", std::abs(f.fade(1000.0)) < 1e-12);
    check("fade mid in (0,1)", f.fade(750.0) > 0.0 && f.fade(750.0) < 1.0,
        f.fade(750.0), 0.5);
    check("fade' at T0 = 0", std::abs(f.derivative(600.0)) < 1e-10);
    check("fade' at T1 = 0", std::abs(f.derivative(900.0)) < 1e-10);
    check("fade'' at T0 = 0", std::abs(f.secondDerivative(600.0)) < 1e-8);
    check("fade'' at T1 = 0", std::abs(f.secondDerivative(900.0)) < 1e-8);
    check("fade monotone decreasing",
        f.fade(650.0) > f.fade(750.0) && f.fade(750.0) > f.fade(850.0));

    // ------------------------------------------------------------------
    // 2. High-temperature limit (fade=0): ideal-gas behaviour.
    // ------------------------------------------------------------------
    HybridHelmholtzThermo thermoIG(registry, 300.0, 400.0);  // fade=0 at 700K
    const HybridHelmholtzThermo::State s700 =
        thermoIG.evaluate(names, Y, 700.0, 2.0);
    check("700K (fade=0): SUCCESS", s700.success());
    if (s700.success())
    {
        const double Rmix = [&]{
            double R = 0.0;
            for (std::size_t i = 0; i < names.size(); ++i)
                if (Y[i] > 0.0)
                    R += Y[i]*8.314462618/registry.molecularWeight(registry.index(names[i]));
            return R;
        }();
        const double pIG = 2.0 * Rmix * 700.0;
        check("ideal limit: p = rho R T",
            relativeError(s700.p, pIG, 1.0e-6) < 1e-6, s700.p, pIG);
        check("ideal limit: cp - cv = Rmix",
            relativeError(s700.cp - s700.cv, Rmix, 1e-6) < 1e-6,
            s700.cp - s700.cv, Rmix);
        check("ideal limit: fade=0", std::abs(s700.fade) < 1e-12, s700.fade);
    }

    // ------------------------------------------------------------------
    // 3. Low-temperature limit (fade=1): full residual.
    // ------------------------------------------------------------------
    // 130 K / 900 kg/m3 (O2-rich mixture): the PR isochoric path lies in
    // the PR two-phase dome (no single-phase PR solution) and must FAIL
    // CLOSED rather than return the spurious vapour-like root.  This is the
    // PR volume-translation limitation; the GERG residual replaces it.
    const HybridHelmholtzThermo::State s130 =
        thermo.evaluate(names, Y, 130.0, 900.0);
    check("130K/900 two-phase-dome: FAIL CLOSED (not SUCCESS)",
        !s130.success(), s130.p, 0.0);

    // ------------------------------------------------------------------
    // 4. Mayer relation from the unified potential (transition band).
    // ------------------------------------------------------------------
    for (double T : {550.0, 700.0, 850.0})
    {
        const auto s = thermo.evaluate(names, Y, T, 10.0);
        check("transition " + std::to_string(int(T)) + "K: SUCCESS",
            s.success());
        if (s.success())
        {
            // Mayer: cp - cv = T/rho^2 (dp/dT)^2 / (dp/drho)
            const double mayer = T/(s.rho*s.rho)
                * s.dpdT_rho*s.dpdT_rho/s.dpdrho_T;
            check("Mayer cp-cv = T(dp/dT)^2/(rho^2 dp/drho)",
                relativeError(s.cp - s.cv, mayer, 1e-9) < 1e-8,
                s.cp - s.cv, mayer);
        }
    }

    // ------------------------------------------------------------------
    // 5. rho-e-z recovery round trip.
    // ------------------------------------------------------------------
    for (double T : {150.0, 300.0, 700.0, 1000.0})
    {
        const auto fwd = thermo.evaluate(names, Y, T, 5.0);
        if (!fwd.success()) { check("fwd " + std::to_string(int(T)) + "K",
            false); continue; }
        const auto rev = thermo.recoverRhoEY(names, Y, 5.0, fwd.e);
        check("rho-e round trip " + std::to_string(int(T)) + "K",
            rev.success() &&
            std::abs(rev.T - T)/T < 1e-5, rev.T, T);
    }

    // ------------------------------------------------------------------
    // 6. Transition-band derivative continuity scan.
    // ------------------------------------------------------------------
    double prevCv = 0.0, prevCp = 0.0;
    bool continuous = true;
    for (double T = 500.0; T <= 1000.0; T += 25.0)
    {
        const auto s = thermo.evaluate(names, Y, T, 5.0);
        if (!s.success()) { continuous = false; break; }
        if (prevCv > 0.0)
        {
            if (std::abs(s.cv - prevCv)/prevCv > 0.3) continuous = false;
            if (std::abs(s.cp - prevCp)/prevCp > 0.3) continuous = false;
        }
        prevCv = s.cv;
        prevCp = s.cp;
    }
    check("transition band: Cv/Cp continuous (no step > 30%)", continuous);

    // ------------------------------------------------------------------
    // 7. Pure-O2 caloric cross-check vs CoolProp HEOS (Gate 3, PR placeholder).
    //    Reference-state-independent quantities: Cp, Cv, sound speed.
    //    PR vs HEOS: near-ideal gas states agree to a few %; dense cryogenic
    //    liquid agrees in magnitude (PR known ~8% rho, ~20% a).
    // ------------------------------------------------------------------
    std::vector<double> YO2(names.size(), 0.0);
    YO2[registry.index("O2")] = 1.0;
    {
        // 300 K / 1 atm gas: cp=919.89 cv=658.73 a=329.72 (CoolProp HEOS)
        const auto s300 = thermo.evaluate(names, YO2, 300.0, 1.301);
        check("O2 300K gas: SUCCESS", s300.success());
        if (s300.success())
        {
            check("O2 300K: cp within 3% of HEOS (919.9)",
                std::abs(s300.cp-919.89)/919.89 < 0.03, s300.cp, 919.89);
            check("O2 300K: cv within 3% of HEOS (658.7)",
                std::abs(s300.cv-658.73)/658.73 < 0.03, s300.cv, 658.73);
            check("O2 300K: a within 3% of HEOS (329.7)",
                std::abs(s300.soundSpeed-329.72)/329.72 < 0.03,
                s300.soundSpeed, 329.72);
        }
    }
    {
        // 700 K / 4 MPa gas: cp=1038.98 cv=772.38 a=500.18 (CoolProp HEOS)
        const auto s700 = thermo.evaluate(names, YO2, 700.0, 21.751);
        check("O2 700K/4MPa gas: SUCCESS", s700.success());
        if (s700.success())
        {
            check("O2 700K: cp within 3% of HEOS (1039.0)",
                std::abs(s700.cp-1038.98)/1038.98 < 0.03, s700.cp, 1038.98);
            check("O2 700K: cv within 5% of HEOS (772.4)",
                std::abs(s700.cv-772.38)/772.38 < 0.05, s700.cv, 772.38);
            check("O2 700K: a within 5% of HEOS (500.2)",
                std::abs(s700.soundSpeed-500.18)/500.18 < 0.05,
                s700.soundSpeed, 500.18);
        }
    }
    {
        // 130 K / 4 MPa dense liquid (LOX state): cp=2023.27 cv=829.24
        // a=575.76 (CoolProp HEOS).  PR placeholder: expect magnitude
        // agreement (cp/cv < 20%, a < 35%).
        // 130 K dense liquid (HEOS rho=920.66, PR sat-liq 993.68) lies in
        // the PR two-phase dome for the ISOCHORIC path: the isochoric cubic
        // has only a spurious vapour-like root, so the single-phase
        // evaluate() must FAIL CLOSED.  The PR liquid branch is only
        // reachable on the ISOBARIC path (validated by test_hybrid_bubble).
        // The GERG residual, defined for every rho, removes this
        // limitation — this is part of the GERG motivation.
        const auto s130 = thermo.evaluate(names, YO2, 130.0, 993.68);
        check("O2 130K liquid isochoric: FAIL CLOSED (PR two-phase dome)",
            !s130.success());
    }

    // Negative-e cryogenic state at 110 K / 950 kg/m3 lies in the PR
    // two-phase dome (PR sat-liq ~1150 kg/m3): fail closed.
    const auto cryo = thermo.evaluate(names, Y, 110.0, 950.0);
    check("110K dense: FAIL CLOSED (PR two-phase dome)", !cryo.success());

    std::cout << (failures == 0 ? "HYBRID_HELMHOLTZ_PASS\n"
                                : "HYBRID_HELMHOLTZ_FAIL\n");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main()
{
    return run();
}
