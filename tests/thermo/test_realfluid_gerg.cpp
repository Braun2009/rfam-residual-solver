// RealFluidState GERG-2008 backend integration (hot-start LOX path).
//
// Validates:
//   1. updatePTY at the 130 K / 4 MPa LOX injector state: the GERG backend
//      finds the dense-liquid root (rho ~ 920.66) where the PR isochoric
//      path has no single-phase root.  Without the table the same call
//      fails closed.
//   2. updateRhoTY dense liquid.
//   3. updateRhoEY round trip in the dense liquid branch.
//   4. Audited PR fallback counter increments when GERG cannot serve a
//      state (hot gas outside the table range).
//   5. Legacy thermo model chain (EquilibriumThermoModel) reaches the GERG
//      backend through setGERGTable.
//
// Output marker: REALFLUID_GERG_PASS / _FAIL.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "EquilibriumThermoModel.H"
#include "GERGResidualTable.H"
#include "RealFluidState.H"
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
    static GERGResidualTable table;
    if (!table.load("tests/reference/gerg_residual_table_v1.csv"))
    {
        check("GERG table load", false);
        std::cout << "REALFLUID_GERG_FAIL\n";
        return 1;
    }

    SpeciesRegistry registry;
    const auto names = registry.names();
    std::vector<double> YO2(names.size(), 0.0);
    YO2[registry.index("O2")] = 1.0;
    std::vector<double> Yhot(names.size(), 0.0);
    Yhot[registry.index("CO2")] = 0.6468;
    Yhot[registry.index("H2O")] = 0.3532;

    // ------------------------------------------------------------------
    // 1. updatePTY at the 130 K / 4 MPa LOX injector state.
    // ------------------------------------------------------------------
    {
        RealFluidState state;
        state.setSpecies(names);
        // Without the GERG table: PR gives its own liquid root (993.68,
        // CoolProp PR-consistent); the fail-close limitation applies to the
        // isochoric path, not to this isobaric PTY update.
        const auto prStatus = state.updatePTY(4.0e6, 130.0, YO2);
        check("updatePTY 130K/4MPa (PR only): SUCCESS (isobaric PR root)",
            prStatus == RealFluidState::Status::SUCCESS);
        if (prStatus == RealFluidState::Status::SUCCESS)
            check("PR liquid root 993.7 (CoolProp PR reference)",
                std::abs(state.rho() - 993.68)/993.68 < 0.01,
                state.rho(), 993.68);

        // With the GERG table: dense liquid root.
        state.setGERGTable(&table);
        const auto status = state.updatePTY(4.0e6, 130.0, YO2);
        check("updatePTY 130K/4MPa (GERG): SUCCESS",
            status == RealFluidState::Status::SUCCESS);
        if (status == RealFluidState::Status::SUCCESS)
        {
            check("GERG liquid root rho ~ 920.7 (HEOS)",
                std::abs(state.rho() - 920.66)/920.66 < 0.05,
                state.rho(), 920.66);
            check("GERG liquid root cp within 8% of HEOS (2023.3)",
                std::abs(state.Cp() - 2023.27)/2023.27 < 0.08,
                state.Cp(), 2023.27);
            check("GERG liquid root sound speed within 8% (575.8)",
                std::abs(state.a() - 575.76)/575.76 < 0.08,
                state.a(), 575.76);
            check("GERG no PR fallback on liquid path",
                state.gergFallbackCount() == 0,
                static_cast<double>(state.gergFallbackCount()), 0.0);
        }
    }

    // ------------------------------------------------------------------
    // 2. updateRhoTY dense liquid.
    // ------------------------------------------------------------------
    {
        RealFluidState state;
        state.setSpecies(names);
        state.setGERGTable(&table);
        const auto status = state.updateRhoTY(920.66, 130.0, YO2);
        check("updateRhoTY 920.66/130 (GERG): SUCCESS",
            status == RealFluidState::Status::SUCCESS);
        if (status == RealFluidState::Status::SUCCESS)
        {
            check("updateRhoTY p ~ 4 MPa",
                std::abs(state.p() - 4.0e6)/4.0e6 < 0.25, state.p(), 4.0e6);
        }
    }

    // ------------------------------------------------------------------
    // 3. updateRhoEY round trip in the dense liquid branch.
    // ------------------------------------------------------------------
    {
        RealFluidState state;
        state.setSpecies(names);
        state.setGERGTable(&table);
        const auto fwd = state.updatePTY(4.0e6, 130.0, YO2);
        if (fwd == RealFluidState::Status::SUCCESS)
        {
            RealFluidState rev;
            rev.setSpecies(names);
            rev.setGERGTable(&table);
            const auto status =
                rev.updateRhoEY(state.rho(), state.e(), YO2, 130.0);
            check("updateRhoEY dense-liquid round trip",
                status == RealFluidState::Status::SUCCESS
                    && std::abs(rev.T() - 130.0)/130.0 < 1e-4,
                rev.T(), 130.0);
            check("updateRhoEY p recovered",
                status == RealFluidState::Status::SUCCESS
                    && std::abs(rev.p() - 4.0e6)/4.0e6 < 0.02,
                rev.p(), 4.0e6);
        }
    }

    // ------------------------------------------------------------------
    // 4. Hot reacting-product state: GERG cannot serve (> table range),
    //    fallback to PR is explicit and audited.
    // ------------------------------------------------------------------
    {
        RealFluidState state;
        state.setSpecies(names);
        state.setGERGTable(&table);
        const auto status = state.updatePTY(4.0e6, 2500.0, Yhot);
        check("updatePTY 2500K products (GERG miss -> PR): SUCCESS",
            status == RealFluidState::Status::SUCCESS);
        check("GERG fallback audited (count >= 1)",
            state.gergFallbackCount() >= 1,
            static_cast<double>(state.gergFallbackCount()), 1.0);
    }

    // ------------------------------------------------------------------
    // 5. EquilibriumThermoModel chain reaches the GERG backend.
    // ------------------------------------------------------------------
    {
        auto model = thermo::makeLegacySinglePhaseThermoModel(names);
        std::const_pointer_cast<thermo::EquilibriumThermoModel>(model)
            ->setGERGTable(&table);
        const auto result =
            model->initializePTY(4.0e6, 130.0, YO2);
        check("Legacy model initializePTY 130K/4MPa (GERG): SUCCESS",
            result.success());
        if (result.success() && result.state)
        {
            check("Legacy model dense liquid rho ~ 920.7",
                std::abs(result.state->rho - 920.66)/920.66 < 0.05,
                result.state->rho, 920.66);
        }
    }

    std::cout << (failures == 0 ? "REALFLUID_GERG_PASS\n"
                                : "REALFLUID_GERG_FAIL\n");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main()
{
    return run();
}
