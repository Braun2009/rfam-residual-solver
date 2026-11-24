// RealFluidState minimal HEM two-phase recovery.
//
// A (rho, e) state inside the GERG two-phase dome (e.g. 120 K O2 liquid
// mixed by diffusion toward the vapour) must recover through the minimal
// HEM path instead of failing: T -> Tsat, quality from the energy lever
// rule, with the two-phase flag set.
//
// Test plan:
//   1. Pure O2 two-phase box: pick (T=100 K, x=0.4) from the saturation
//      table endpoints (CoolProp HEOS); ask updateRhoEY with (rho_mix,
//      e_mix): must return SUCCESS, two-phase, quality ~ 0.4, T ~ 100 K.
//   2. Compressed liquid (130 K/4 MPa) still single-phase (regression).
//   3. Hot gas (300 K) still single-phase.
//
// Output marker: HEM_TWO_PHASE_PASS / _FAIL.

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
        std::cout << "HEM_TWO_PHASE_FAIL\n";
        return 1;
    }
    SpeciesRegistry registry;
    const auto names = registry.names();
    std::vector<double> YO2(names.size(), 0.0);
    YO2[registry.index("O2")] = 1.0;

    // 1. Two-phase box at 100 K, quality x=0.4.
    {
        const double T = 100.0;
        const double pSat = sat.saturationPressure(T, 1.0);
        check("O2 100K p_sat finite", std::isfinite(pSat) && pSat > 0.0,
            pSat / 1e3);
        double rhoV = 0.0, rhoL = 0.0, eV = 0.0, eL = 0.0;
        const bool endpoints = sat.saturationEndpoints(
            T, 1.0, rhoV, rhoL, &eV, &eL);
        check("O2 100K saturation endpoints", endpoints);
        if (endpoints)
        {
            const double x = 0.4;
            const double rhoMix = 1.0 / (x / rhoV + (1.0 - x) / rhoL);
            const double eMix = x * eV + (1.0 - x) * eL;
            RealFluidState rec;
            rec.setSpecies(names);
            rec.setGERGTable(&gerg);
            rec.setSaturationTable(&sat);
            const auto status = rec.updateRhoEY(rhoMix, eMix, YO2);
            check("HEM recovery SUCCESS + two-phase",
                status == RealFluidState::Status::SUCCESS && rec.isTwoPhase());
            if (status == RealFluidState::Status::SUCCESS)
            {
                check("HEM quality ~ 0.4",
                    std::abs(rec.quality() - x) < 0.05,
                    rec.quality(), x);
                check("HEM T ~ 100 K",
                    std::abs(rec.T() - T) / T < 0.02, rec.T(), T);
                check("HEM p ~ p_sat",
                    std::abs(rec.p() - pSat) / pSat < 0.02,
                    rec.p(), pSat);
            }
        }
    }

    // 2. Clamped two-phase box: energy beyond the saturated-vapour value
    //    (super-heated vapour inside the dome) -> quality clamps to 1.
    {
        const double T = 100.0;
        double rhoV = 0.0, rhoL = 0.0, eV = 0.0, eL = 0.0;
        if (sat.saturationEndpoints(T, 1.0, rhoV, rhoL, &eV, &eL))
        {
            // v2 sensible e is negative: super-heat via addition.
            const double eHot = eV + 0.05 * std::abs(eV);
            const double rhoMix = rhoV;         // quality clamps to 1 (vapour)
            RealFluidState rec;
            rec.setSpecies(names);
            rec.setGERGTable(&gerg);
            rec.setSaturationTable(&sat);
            const auto status = rec.updateRhoEY(rhoMix, eHot, YO2);
            check("clamped HEM: SUCCESS + two-phase + q=1",
                status == RealFluidState::Status::SUCCESS && rec.isTwoPhase()
                    && std::abs(rec.quality() - 1.0) < 0.05,
                rec.quality(), 1.0);
        }
    }

    // 3. Binary mixture xO2=0.75 at 120 K: v4 same-pressure endpoints are
    //    available, but a NON-PHYSICAL placeholder box (rho=100,
    //    e=-50000 outside the 120 K band) must still FAIL CLOSED (no
    //    silent wrong answer).
    {
        std::vector<double> Ymix(names.size(), 0.0);
        Ymix[registry.index("O2")] = 0.75;
        Ymix[registry.index("C3H8")] = 0.25;
        double rhoV = 0.0, rhoL = 0.0, eV = 0.0, eL = 0.0;
        const bool ok = sat.saturationEndpoints(120.0, 0.75, rhoV, rhoL, &eV, &eL);
        check("mix xO2=0.75 120K endpoints available (teqp v4)",
            ok);
        if (ok)
        {
            const double rhoMix = 100.0;   // dome-range placeholder
            const double eMix = -50000.0;
            RealFluidState rec;
            rec.setSpecies(names);
            rec.setGERGTable(&gerg);
            rec.setSaturationTable(&sat);
            const auto status = rec.updateRhoEY(rhoMix, eMix, Ymix);
            check("mixture HEM: FAIL CLOSED (no silent wrong answer)",
                status != RealFluidState::Status::SUCCESS
                    || !rec.isTwoPhase());
        }
    }

    // 3b. Binary mixture xO2=0.75 at 130 K (teqp v4 SAME-PRESSURE
    //     endpoints): a two-phase box must recover with quality ~0.35 and
    //     T ~ 130 K.  (v3 had NaN endpoints — CoolProp Q=0/1 unreliable;
    //     v4 evaluates both roots at the common bubble pressure.)
    {
        // mole fraction 0.75 O2 -> mass fraction 0.686 (the HEM band
        // [0.72, 0.78] is in MOLE fraction; 0.75 mass would be 0.805)
        const double yO2mass = 0.75 * 0.0319988
            / (0.75 * 0.0319988 + 0.25 * 0.044097);
        std::vector<double> Ymix(names.size(), 0.0);
        Ymix[registry.index("O2")] = yO2mass;
        Ymix[registry.index("C3H8")] = 1.0 - yO2mass;
        double rhoV = 0.0, rhoL = 0.0, eV = 0.0, eL = 0.0;
        const bool endpoints = sat.saturationEndpoints(
            130.0, 0.75, rhoV, rhoL, &eV, &eL);
        check("mix xO2=0.75 130K endpoints available (teqp v4)",
            endpoints && std::isfinite(eV) && std::isfinite(eL));
        if (endpoints && std::isfinite(eV) && std::isfinite(eL))
        {
            const double x = 0.35;
            const double rhoMix = 1.0 / (x / rhoV + (1.0 - x) / rhoL);
            const double eMix = x * eV + (1.0 - x) * eL;
            RealFluidState rec;
            rec.setSpecies(names);
            rec.setGERGTable(&gerg);
            rec.setSaturationTable(&sat);
            const auto status = rec.updateRhoEY(rhoMix, eMix, Ymix);
            check("mix xO2=0.75 130K HEM recovery SUCCESS + two-phase",
                status == RealFluidState::Status::SUCCESS
                    && rec.isTwoPhase(),
                static_cast<int>(status));
            if (status == RealFluidState::Status::SUCCESS)
            {
                check("mix xO2=0.75 130K quality ~ 0.35",
                    std::abs(rec.quality() - x) < 0.05,
                    rec.quality(), x);
                check("mix xO2=0.75 130K T ~ 130",
                    std::abs(rec.T() - 130.0) / 130.0 < 0.02,
                    rec.T(), 130.0);
            }
        }
    }

    // 4. Compressed liquid still single-phase (regression).
    {
        RealFluidState st;
        st.setSpecies(names);
        st.setGERGTable(&gerg);
        st.setSaturationTable(&sat);
        const auto status = st.updatePTY(4.0e6, 130.0, YO2);
        check("130K/4MPa compressed liquid: single-phase",
            status == RealFluidState::Status::SUCCESS && !st.isTwoPhase());
    }

    // 3. Hot gas single-phase (regression).
    {
        RealFluidState st;
        st.setSpecies(names);
        st.setGERGTable(&gerg);
        st.setSaturationTable(&sat);
        const auto status = st.updatePTY(4.0e6, 300.0, YO2);
        check("300K gas: single-phase",
            status == RealFluidState::Status::SUCCESS && !st.isTwoPhase());
    }

    std::cout << (failures == 0 ? "HEM_TWO_PHASE_PASS\n"
                                : "HEM_TWO_PHASE_FAIL\n");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main()
{
    return run();
}
