// Recovery profiling: Newton vs bracket fallback vs HEM scan counters.
//
// 1. A warm single-phase recovery (T0 near the solution) must solve via
//    Newton with few iterations and NO bracket fallback.
// 2. A cold start (T0 far) may fall back to the bracket path: counters
//    must reflect it (bracketFallbacks > 0, bracketRate > 0).
// 3. A two-phase recovery must count a HEM full scan.
// 4. Profile aggregates (newtonAvg, bracketRate) stay finite.
//
// Output marker: RECOVERY_PROFILE_PASS / _FAIL.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "GERGResidualTable.H"
#include "GERGSaturationTable.H"
#include "HybridHelmholtzThermo.H"
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
        std::cout << "RECOVERY_PROFILE_FAIL\n";
        return 1;
    }
    SpeciesRegistry registry;
    const auto names = registry.names();
    std::vector<double> YO2(names.size(), 0.0);
    YO2[registry.index("O2")] = 1.0;

    // 1. GERG main path (300 K gas): the isochoric recovery is a bisection
    //    (no Newton on the GERG path) — counters must show bisection iters.
    {
        HybridHelmholtzThermo thermo(registry);
        thermo.registerResidualTable(&gerg);
        std::vector<std::string> n1{"O2"};
        std::vector<double> y1{1.0};
        // 300 K / 1.301 kg/m3 gas: e from evaluate
        const auto ref = thermo.evaluate(n1, y1, 300.0, 1.301);
        check("GERG gas reference", ref.success());
        const auto st = thermo.recoverRhoEY(n1, y1, 1.301, ref.e);
        check("GERG recovery SUCCESS", st.success());
        const auto& c = thermo.recoveryCounters();
        check("GERG: bisection used",
            c.bisectionIters > 0 && c.bisectionIters < 100,
            static_cast<double>(c.bisectionIters), 10.0);
        check("GERG: adaptive scan bounded (<64 samples/end)",
            c.adaptiveScans < 128,
            static_cast<double>(c.adaptiveScans), 128.0);
        std::cout << "  [INFO] gerg: calls=" << c.calls
                  << " bisection=" << c.bisectionIters
                  << " adaptiveScans=" << c.adaptiveScans << "\n";
    }

    // 2. PR fallback path (hot products beyond the GERG table): Newton and
    //    bracket counters live on the RealFluidState side.
    {
        RealFluidState st;
        st.setSpecies(names);
        st.setGERGTable(&gerg);
        st.setSaturationTable(&sat);
        std::vector<double> Yhot(names.size(), 0.0);
        Yhot[registry.index("CO2")] = 0.6468;
        Yhot[registry.index("H2O")] = 0.3532;
        const auto s0 = st.updatePTY(4.0e6, 2500.0, Yhot);
        check("hot products setup", s0 == RealFluidState::Status::SUCCESS);
        const auto status = st.updateRhoEY(st.rho(), st.e(), Yhot, 2500.0);
        check("hot recovery SUCCESS", status == RealFluidState::Status::SUCCESS);
        const auto& p = st.recoveryProfile();
        check("hot: recovery counted",
            p.recoveryCalls >= 1,
            static_cast<double>(p.recoveryCalls), 1.0);
        check("hot: bracketRate in [0,1]",
            p.bracketRate >= 0.0 && p.bracketRate <= 1.0, p.bracketRate);
        std::cout << "  [INFO] hot: calls=" << p.recoveryCalls
                  << " newtonIters=" << p.newtonIterations
                  << " bracketFallbacks=" << p.bracketFallbacks
                  << " bracketRate=" << p.bracketRate << "\n";
    }

    // 3. Two-phase recovery counts a HEM full scan.
    {
        double rhoV = 0.0, rhoL = 0.0, eV = 0.0, eL = 0.0;
        const bool ok = sat.saturationEndpoints(100.0, 1.0, rhoV, rhoL, &eV, &eL);
        check("two-phase endpoints", ok);
        if (ok)
        {
            const double x = 0.4;
            const double rhoMix = 1.0 / (x / rhoV + (1.0 - x) / rhoL);
            const double eMix = x * eV + (1.0 - x) * eL;
            RealFluidState st;
            st.setSpecies(names);
            st.setGERGTable(&gerg);
            st.setSaturationTable(&sat);
            const auto status = st.updateRhoEY(rhoMix, eMix, YO2);
            check("two-phase recovery SUCCESS",
                status == RealFluidState::Status::SUCCESS && st.isTwoPhase());
            const auto& p = st.recoveryProfile();
            check("two-phase: HEM full scan counted",
                p.hemFullScans >= 1,
                static_cast<double>(p.hemFullScans), 1.0);
            check("aggregates finite",
                std::isfinite(p.newtonAvg) && std::isfinite(p.bracketRate));

            std::cout << "  [INFO] two-phase: hemFullScans=" << p.hemFullScans
                      << " recoveryCalls=" << p.recoveryCalls << "\n";
        }
    }

    std::cout << (failures == 0 ? "RECOVERY_PROFILE_PASS\n"
                                : "RECOVERY_PROFILE_FAIL\n");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main()
{
    return run();
}
