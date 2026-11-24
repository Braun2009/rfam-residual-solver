// HybridHelmholtzThermo single-component bubble-point chain test.
//
// Validates the isobaric root + unified-potential chain (PR placeholder):
//   - single-component bubble pressure of O2 at 110.22 K must match the
//     CoolProp PR reference (554.2 kPa) and HEOS (551.7 kPa);
//   - C3H8 at 200 K: CoolProp PR ~ 20.2 kPa;
//   - at the bubble pressure the isobaric roots are liquid/vapour density
//     magnitudes and g_L == g_V.
//
// Output marker: HYBRID_BUBBLE_PASS / _FAIL.

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

int run()
{
    SpeciesRegistry registry;
    const std::vector<std::string> names = registry.names();
    std::vector<double> YO2(names.size(), 0.0);
    YO2[registry.index("O2")] = 1.0;
    HybridHelmholtzThermo thermo(registry, 600.0, 900.0);

    // ------------------------------------------------------------------
    // O2 bubble at 110.22 K: CoolProp PR 554.2 kPa, HEOS 551.7 kPa.
    // ------------------------------------------------------------------
    const double pBubO2 = thermo.bubblePressure(names, YO2, 110.22);
    check("O2 bubble found", std::isfinite(pBubO2), pBubO2/1e3, 554.2);
    check("O2 bubble ~ CoolProp PR 554.2 kPa (within 5%)",
        std::abs(pBubO2 - 554.2e3)/554.2e3 < 0.05, pBubO2/1e3, 554.2);
    check("O2 bubble ~ HEOS 551.7 kPa (within 5%)",
        std::abs(pBubO2 - 551.7e3)/551.7e3 < 0.05, pBubO2/1e3, 551.7);

    // At the bubble the isobaric roots must be liquid/vapour magnitudes and
    // g_L == g_V.
    {
        const auto roots = thermo.isobaricRoots(names, YO2, 110.22, pBubO2);
        check("two-phase roots at bubble", roots.valid && roots.twoPhase);
        if (roots.valid && roots.twoPhase)
        {
            check("liquid root dense (>800 kg/m3)",
                roots.rhoLiquidKgM3 > 800.0, roots.rhoLiquidKgM3, 1150.5);
            check("vapour root light (<100 kg/m3)",
                roots.rhoVaporKgM3 < 100.0, roots.rhoVaporKgM3, 21.7);
            const auto sL = thermo.evaluate(
                names, YO2, 110.22, roots.rhoLiquidKgM3);
            const auto sV = thermo.evaluate(
                names, YO2, 110.22, roots.rhoVaporKgM3);
            if (sL.success() && sV.success())
            {
                const double gL = (sL.e - 110.22*sL.s)
                    + sL.p/roots.rhoLiquidKgM3;
                const double gV = (sV.e - 110.22*sV.s)
                    + sV.p/roots.rhoVaporKgM3;
                check("g_L == g_V at bubble (within 1e-4)",
                    std::abs(gL - gV)/std::max(std::abs(gL), 1.0) < 1e-4,
                    gL - gV, 0.0);
            }
        }
    }

    // ------------------------------------------------------------------
    // C3H8 bubble at 200 K: CoolProp PR ~ 20.2 kPa.
    // ------------------------------------------------------------------
    std::vector<double> YC3(names.size(), 0.0);
    YC3[registry.index("C3H8")] = 1.0;
    const double pBubC3 = thermo.bubblePressure(names, YC3, 200.0);
    check("C3H8 bubble found", std::isfinite(pBubC3), pBubC3/1e3, 20.2);
    check("C3H8 bubble ~ CoolProp PR 20.2 kPa (within 15%)",
        std::abs(pBubC3 - 20.2e3)/20.2e3 < 0.15, pBubC3/1e3, 20.2);

    std::cout << (failures == 0 ? "HYBRID_BUBBLE_PASS\n"
                                : "HYBRID_BUBBLE_FAIL\n");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main()
{
    return run();
}
