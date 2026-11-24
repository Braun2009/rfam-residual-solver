// Two-phase continuation hint: a seeded (T, quality) local scan must hit
// for the same state, and a wrong hint must fall back to the full scan.
//
// Output marker: HEM_HINT_PASS / _FAIL.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "GERGResidualTable.H"
#include "GERGSaturationTable.H"
#include "RealFluidState.H"
#include "SpeciesRegistry.H"

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
    static GERGResidualTable gerg;
    static GERGSaturationTable sat;
    if (!gerg.load("tests/reference/gerg_residual_table_v3.csv")
        || !sat.load("tests/reference/gerg_saturation_table_v4.csv"))
    {
        check("tables load", false);
        std::cout << "HEM_HINT_FAIL\n";
        return 1;
    }
    SpeciesRegistry registry;
    const auto names = registry.names();
    std::vector<double> YO2(names.size(), 0.0);
    YO2[registry.index("O2")] = 1.0;

    // Reference two-phase state at 100 K, q=0.4.
    double rhoV = 0.0, rhoL = 0.0, eV = 0.0, eL = 0.0;
    if (!sat.saturationEndpoints(100.0, 1.0, rhoV, rhoL, &eV, &eL))
    {
        check("endpoints", false);
        std::cout << "HEM_HINT_FAIL\n";
        return 1;
    }
    const double x = 0.4;
    const double rhoMix = 1.0 / (x / rhoV + (1.0 - x) / rhoL);
    const double eMix = x * eV + (1.0 - x) * eL;

    // 1. Correct hint: seeded (100 K, 0.4) -> local scan hits.
    {
        RealFluidState st;
        st.setSpecies(names);
        st.setGERGTable(&gerg);
        st.setSaturationTable(&sat);
        st.setTwoPhaseHint(100.0, 0.4);
        const auto status = st.updateRhoEY(rhoMix, eMix, YO2);
        check("hint recovery SUCCESS",
            status == RealFluidState::Status::SUCCESS && st.isTwoPhase());
        const auto& p = st.recoveryProfile();
        check("hint: local scan hit (no full scan)",
            p.hemHintHits >= 1 && p.hemFullScans == 0,
            static_cast<double>(p.hemHintHits), 1.0);
        check("hint: quality preserved ~ 0.4",
            std::abs(st.quality() - x) < 0.05, st.quality(), x);
        std::cout << "  [INFO] hint: hits=" << p.hemHintHits
                  << " fullScans=" << p.hemFullScans << "\n";
    }

    // 2. Wrong hint (200 K, q=1 — far): local scan misses, full scan
    //    falls back and still recovers the same state.
    {
        RealFluidState st;
        st.setSpecies(names);
        st.setGERGTable(&gerg);
        st.setSaturationTable(&sat);
        st.setTwoPhaseHint(200.0, 1.0);
        const auto status = st.updateRhoEY(rhoMix, eMix, YO2);
        check("wrong-hint recovery SUCCESS (fallback)",
            status == RealFluidState::Status::SUCCESS && st.isTwoPhase());
        const auto& p = st.recoveryProfile();
        check("wrong-hint: full scan fallback used",
            p.hemFullScans >= 1,
            static_cast<double>(p.hemFullScans), 1.0);
        check("wrong-hint: quality still ~ 0.4 (hint not authoritative)",
            std::abs(st.quality() - x) < 0.08, st.quality(), x);
        std::cout << "  [INFO] wrong-hint: hits=" << p.hemHintHits
                  << " fullScans=" << p.hemFullScans << "\n";
    }

    std::cout << (failures == 0 ? "HEM_HINT_PASS\n" : "HEM_HINT_FAIL\n");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main()
{
    return run();
}
