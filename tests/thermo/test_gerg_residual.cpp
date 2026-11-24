// HybridHelmholtzThermo with the GERG-2008 residual table.
//
// Validates the GERG-residual integration:
//   1. With the GERG table registered, the dense liquid branch (where the
//      PR isochoric path fails closed) becomes a valid single-phase state:
//      O2 130 K / 920.66 kg/m3 (LOX injector state) must SUCCEED.
//   2. Caloric quantities (Cp, Cv, sound speed) must approach the CoolProp
//      HEOS reference at the dense liquid state (GERG is the HEOS-family
//      reference, so agreement should be far better than the PR placeholder).
//   3. rho-e recovery works in the dense liquid branch.
//   4. Without the table, the same dense state still fails closed (PR
//      placeholder limitation).
//
// Output marker: GERG_RESIDUAL_PASS / _FAIL.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "GERGResidualTable.H"
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
    std::vector<double> YO2(names.size(), 0.0);
    YO2[registry.index("O2")] = 1.0;

    // Load the GERG residual table.
    static GERGResidualTable table;
    if (!table.load("tests/reference/gerg_residual_table_v3.csv"))
    {
        check("GERG table load", false);
        std::cout << "GERG_RESIDUAL_FAIL\n";
        return 1;
    }
    check("GERG table load", table.has("O2") && table.has("C3H8"));

    HybridHelmholtzThermo thermo(registry, 600.0, 900.0);

    // ------------------------------------------------------------------
    // 1. Dense liquid without the table: PR placeholder fails closed.
    // ------------------------------------------------------------------
    {
        const auto s = thermo.evaluate(names, YO2, 130.0, 920.66);
        check("130K/920.66 (PR placeholder): FAIL CLOSED", !s.success());
    }

    // ------------------------------------------------------------------
    // 2. With the GERG table: dense liquid becomes a valid state.
    // ------------------------------------------------------------------
    thermo.registerResidualTable(&table);
    {
        const auto s = thermo.evaluate(names, YO2, 130.0, 920.66);
        check("130K/920.66 (GERG): SUCCESS", s.success());
        if (s.success())
        {
            // CoolProp HEOS reference at 130 K / 4 MPa / rho=920.66:
            //   cp=2023.27  cv=829.24  a=575.76
            check("GERG: cp within 8% of HEOS (2023.3)",
                relativeError(s.cp, 2023.27) < 0.08, s.cp, 2023.27);
            check("GERG: cv within 8% of HEOS (829.2)",
                relativeError(s.cv, 829.24) < 0.08, s.cv, 829.24);
            check("GERG: a within 8% of HEOS (575.8)",
                relativeError(s.soundSpeed, 575.76) < 0.08,
                s.soundSpeed, 575.76);
            check("GERG: mechanical stability dp/drho > 0",
                s.dpdrho_T > 0.0, s.dpdrho_T);
            check("GERG: pressure finite positive",
                std::isfinite(s.p) && s.p > 0.0, s.p, 4.0e6);
        }
    }

    // ------------------------------------------------------------------
    // 3. rho-e recovery in the dense liquid branch.
    // ------------------------------------------------------------------
    {
        const auto fwd = thermo.evaluate(names, YO2, 130.0, 920.66);
        if (fwd.success())
        {
            // 110-120 K at 920.66 kg/m3 lie in the two-phase dome (correctly
            // fail closed); the 130 K compressed-liquid bracket is stable.
            const auto rev = thermo.recoverRhoEY(
                names, YO2, 920.66, fwd.e, 130.0, 500.0);
            check("GERG: dense-liquid rho-e round trip",
                rev.success() && std::abs(rev.T - 130.0)/130.0 < 1e-4,
                rev.T, 130.0);
        }
    }

    // ------------------------------------------------------------------
    // 4. Gas branch still consistent (300 K / 1 atm).
    // ------------------------------------------------------------------
    {
        const auto s = thermo.evaluate(names, YO2, 300.0, 1.301);
        check("300K gas (GERG): SUCCESS", s.success());
        if (s.success())
        {
            check("300K gas: cp within 3% of HEOS (919.9)",
                relativeError(s.cp, 919.89) < 0.03, s.cp, 919.89);
            check("300K gas: a within 3% of HEOS (329.7)",
                relativeError(s.soundSpeed, 329.72) < 0.03,
                s.soundSpeed, 329.72);
        }
    }

    // ------------------------------------------------------------------
    // 4b. O2/C3H8 binary mixture route: nearest fixed-composition table.
    // ------------------------------------------------------------------
    {
        std::vector<double> Ymix(names.size(), 0.0);
        Ymix[registry.index("O2")] = 0.5;
        Ymix[registry.index("C3H8")] = 0.5;
        const auto s = thermo.evaluate(names, Ymix, 300.0, 10.0);
        check("300K O2/C3H8 50/50 mix (GERG): SUCCESS", s.success());
        if (s.success())
            check("mix: mechanical stability", s.dpdrho_T > 0.0, s.dpdrho_T);
        // Dense single-phase mix (300 K compressed): GERG defined.  Note
        // 130 K / 600 lies in the two-phase dome for xO2=0.5 (correctly
        // fail-closed by the mechanical-stability gate).
        const auto sd = thermo.evaluate(names, Ymix, 300.0, 600.0);
        check("300K dense mix (GERG): SUCCESS", sd.success());
        if (sd.success())
            check("300K dense mix: p finite positive",
                std::isfinite(sd.p) && sd.p > 0.0, sd.p);
        // Two-phase dome mix correctly fails closed.
        const auto dome = thermo.evaluate(names, Ymix, 130.0, 600.0);
        check("130K/600 mix in dome: FAIL CLOSED", !dome.success());
    }

    // ------------------------------------------------------------------
    // 5. Bubble pressure unchanged (GERG table does not affect the PR
    //    isobaric phase-boundary path used for validation).
    // ------------------------------------------------------------------
    {
        const double pBub = thermo.bubblePressure(names, YO2, 110.22);
        check("O2 bubble still ~554 kPa (isobaric PR path)",
            std::abs(pBub - 554.2e3)/554.2e3 < 0.05, pBub/1e3, 554.2);
    }

    std::cout << (failures == 0 ? "GERG_RESIDUAL_PASS\n"
                                : "GERG_RESIDUAL_FAIL\n");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main()
{
    return run();
}
