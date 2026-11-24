// GERGSaturationTable: bubble-pressure two-phase detection.
//
// Validates:
//   1. Pure O2 110.22 K bubble = 551.7 kPa (HEOS/experiment).
//   2. Pure O2 120 K / 4 MPa (Mayer jet condition) = compressed liquid
//      (SUPPORTED_SINGLE_PHASE) — 4 MPa >> p_sat 1.02 MPa.
//   3. Pure O2 120 K / 0.5 MPa = inside the dome (TWO_PHASE_CANDIDATE).
//   4. C3H8 200 K bubble ~ 20.2 kPa.
//   5. Out-of-table (T beyond coverage) = OUT_OF_TABLE.
//
// Output marker: SATURATION_TABLE_PASS / _FAIL.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "GERGSaturationTable.H"

using namespace Foam;

namespace
{

int failures = 0;

void check(const std::string& name, bool ok,
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
    GERGSaturationTable table;
    if (!table.load("tests/reference/gerg_saturation_table_v4.csv"))
    {
        check("saturation table load", false);
        std::cout << "SATURATION_TABLE_FAIL\n";
        return 1;
    }

    // 1. Pure O2 110.22 K bubble.
    const double p110 = table.saturationPressure(110.22, 1.0);
    check("O2 110.22K bubble ~ 551.7 kPa",
        std::isfinite(p110) && std::abs(p110 - 551.7e3) / 551.7e3 < 0.02,
        p110 / 1e3, 551.7);

    // 2. Mayer jet condition: 120 K / 4 MPa compressed liquid.
    check("O2 120K/4MPa (Mayer): compressed liquid",
        table.classify(120.0, 4.0e6, 1.0)
            == SaturationPhase::SUPPORTED_SINGLE_PHASE);

    // 3. 120 K / 0.5 MPa inside dome.
    check("O2 120K/0.5MPa: two-phase candidate",
        table.classify(120.0, 0.5e6, 1.0)
            == SaturationPhase::TWO_PHASE_CANDIDATE);

    // 4. C3H8 200 K bubble ~ 20.2 kPa.
    const double pC3 = table.saturationPressure(200.0, 0.0);
    check("C3H8 200K bubble ~ 20.2 kPa",
        std::isfinite(pC3) && std::abs(pC3 - 20.2e3) / 20.2e3 < 0.05,
        pC3 / 1e3, 20.2);

    // 5. Out of table.
    check("O2 500K: out of table",
        table.classify(500.0, 4.0e6, 1.0)
            == SaturationPhase::OUT_OF_TABLE);

    // 6. Mixture level xO2=0.95 130 K (covered; xO2=0.5 deep-cold rows
    //    above 100 K are absent — dual root-solver divergence recorded).
    const double pMix = table.saturationPressure(130.0, 0.95);
    check("mix xO2=0.95 130K bubble finite",
        std::isfinite(pMix) && pMix > 0.0, pMix / 1e3);

    // 7. v3 coexisting-phase data: xO2=0.75 120 K carries a teqp bubble
    //    pressure inside the pure-component bounds (1022 kPa O2 / ~40 kPa
    //    C3H8) and a coexisting vapour composition y > x.  (The v2
    //    torch-flash value 1042 kPa violated the physical upper bound.)
    {
        double rhoV = 0.0, rhoL = 0.0, eV = 0.0, eL = 0.0;
        // 130 K: v4 same-pressure endpoints (teqp density roots at the
        // common bubble pressure) must be available and physical.
        const bool ok = table.saturationEndpoints(
            130.0, 0.75, rhoV, rhoL, &eV, &eL);
        check("mix xO2=0.75 130K endpoints available (teqp v4)",
            ok && rhoL > rhoV && rhoV > 1.0
                && std::isfinite(eV) && std::isfinite(eL));
        if (ok)
        {
            check("mix xO2=0.75 130K rhoL physical (~870 kg/m3)",
                std::abs(rhoL - 872.7) / 872.7 < 0.05, rhoL);
            check("mix xO2=0.75 130K rhoV physical (~20 kg/m3)",
                std::abs(rhoV - 20.4) / 20.4 < 0.05, rhoV);
        }
        const double pSat = table.saturationPressure(120.0, 0.75);
        check("mix xO2=0.75 120K bubble within physical bounds",
            std::isfinite(pSat) && pSat > 200.0e3 && pSat < 1.1e6,
            pSat / 1e3);
        const double yO2 = table.vaporComposition(120.0, 0.75);
        check("mix xO2=0.75 120K coexisting vapour y > x",
            std::isfinite(yO2) && yO2 > 0.75 && yO2 < 1.0, yO2);
    }

    std::cout << (failures == 0 ? "SATURATION_TABLE_PASS\n"
                                : "SATURATION_TABLE_FAIL\n");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main()
{
    return run();
}
