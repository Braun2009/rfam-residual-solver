#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "RealFluidState.H"

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
    return std::abs(got - expected)
        / std::max({std::abs(got), std::abs(expected), floor});
}

void checkRoundTrip(
    const std::string& label,
    const std::vector<std::string>& species,
    const std::vector<double>& Y,
    double p, double T)
{
    RealFluidState original;
    original.setSpecies(species);
    const auto initialStatus = original.updatePTY(p, T, Y);
    check(label + " PTY succeeds",
        initialStatus == RealFluidState::Status::SUCCESS);
    if (initialStatus != RealFluidState::Status::SUCCESS) return;

    check(label + " Cv positive", original.Cv() > 0.0, original.Cv());
    check(label + " Cp greater than Cv",
        original.Cp() > original.Cv(), original.Cp(), original.Cv());
    check(label + " sound speed finite",
        std::isfinite(original.a()) && original.a() > 0.0, original.a());

    RealFluidState recovered;
    recovered.setSpecies(species);
    const auto recoveredStatus = recovered.updateRhoEY(
        original.rho(), original.e(), Y, 1.05*T);
    check(label + " rhoE succeeds",
        recoveredStatus == RealFluidState::Status::SUCCESS,
        static_cast<double>(recoveredStatus),
        static_cast<double>(RealFluidState::Status::SUCCESS));
    if (recoveredStatus != RealFluidState::Status::SUCCESS) return;

    check(label + " temperature round trip",
        relativeError(recovered.T(), T) < 2.0e-8, recovered.T(), T);
    check(label + " pressure round trip",
        relativeError(recovered.p(), p) < 1.0e-7, recovered.p(), p);
    check(label + " density round trip",
        relativeError(recovered.rho(), original.rho()) < 2.0e-8,
        recovered.rho(), original.rho());
    check(label + " energy round trip",
        relativeError(recovered.e(), original.e()) < 1.0e-8,
        recovered.e(), original.e());
}

} // namespace

int main()
{
    std::cout << "=== RealFluidState validation and inversion contract ===\n";

    using Status = RealFluidState::Status;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    check("valid scalar set succeeds",
        RealFluidState::validateThermodynamicScalars(700.0, 1000.0, 1.0e5)
            == Status::SUCCESS);
    check("negative Cv is preserved as failure",
        RealFluidState::validateThermodynamicScalars(-1.0, 1000.0, 1.0e5)
            == Status::NEGATIVE_CV);
    check("NaN Cv is preserved as failure",
        RealFluidState::validateThermodynamicScalars(nan, 1000.0, 1.0e5)
            == Status::NEGATIVE_CV);
    check("Cp <= Cv is preserved as failure",
        RealFluidState::validateThermodynamicScalars(700.0, 700.0, 1.0e5)
            == Status::NEGATIVE_CP);
    check("invalid sound speed is preserved as failure",
        RealFluidState::validateThermodynamicScalars(700.0, 1000.0, nan)
            == Status::NEGATIVE_SOUNDSPEED);

    {
        RealFluidState state;
        check("NaN pressure fails",
            state.updatePTY(nan, 300.0, {1.0, 0.0}) != Status::SUCCESS);
        check("failed state is not valid", !state.isValid());
        check("negative temperature fails",
            state.updatePTY(1.0e5, -1.0, {1.0, 0.0}) != Status::SUCCESS);
        check("NaN composition fails",
            state.updatePTY(1.0e5, 300.0, {nan, 0.0}) != Status::SUCCESS);
        check("wrong composition size fails",
            state.updatePTY(1.0e5, 300.0, {1.0}) != Status::SUCCESS);
    }

    {
        const double p = 2.0e5;
        const double T = 300.0;
        const std::vector<std::string> species = {"O2"};
        const std::vector<double> Y = {1.0};
        RealFluidState center;
        center.setSpecies(species);
        const auto status = center.updatePTY(p, T, Y);
        check("derivative reference state succeeds", status == Status::SUCCESS);
        if (status == Status::SUCCESS)
        {
            const double identity =
                (center.p() - center.T()*center.dpdT_rho())
                /(center.rho()*center.rho());
            check("dedrho_T identity",
                relativeError(center.dedrho_T(), identity) < 1.0e-12,
                center.dedrho_T(), identity);

            const double drho = 1.0e-4*center.rho();
            RealFluidState lower;
            RealFluidState upper;
            lower.setSpecies(species);
            upper.setSpecies(species);
            const auto lowerStatus = lower.updateRhoTY(
                center.rho() - drho, T, Y);
            const auto upperStatus = upper.updateRhoTY(
                center.rho() + drho, T, Y);
            check("lower density derivative sample succeeds",
                lowerStatus == Status::SUCCESS);
            check("upper density derivative sample succeeds",
                upperStatus == Status::SUCCESS);
            if (lowerStatus == Status::SUCCESS && upperStatus == Status::SUCCESS)
            {
                const double numeric = (upper.e() - lower.e())/(2.0*drho);
                check("dedrho_T finite difference",
                    relativeError(center.dedrho_T(), numeric, 1.0e-7) < 2.0e-4,
                    center.dedrho_T(), numeric);
            }
        }
    }

    checkRoundTrip("O2 liquid", {"O2"}, {1.0}, 4.0e6, 130.0);
    checkRoundTrip("C3H8 dense gas", {"C3H8"}, {1.0}, 5.0e6, 450.0);
    checkRoundTrip(
        "O2/C3H8 mixture", {"O2", "C3H8"},
        {0.784, 0.216}, 4.0e6, 400.0);
    checkRoundTrip(
        "burned mixture", {"O2", "CO2", "H2O", "N2"},
        {0.10, 0.35, 0.20, 0.35}, 4.0e6, 2000.0);
    checkRoundTrip(
        "nine-species hot reactants",
        {"O2", "N2", "H2", "CH4", "C2H6", "C3H8", "H2O", "CO", "CO2"},
        {0.784, 0.0, 0.0, 0.0, 0.0, 0.216, 0.0, 0.0, 0.0},
        4.0e6, 1200.0);
    checkRoundTrip(
        "nine-species cold mixture",
        {"O2", "N2", "H2", "CH4", "C2H6", "C3H8", "H2O", "CO", "CO2"},
        {0.8, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0},
        1.0e6, 500.0);

    {
        const double pressure = 6.0e6;
        const double temperature = 100.0;
        const std::vector<double> composition = {1.0, 0.0};
        RealFluidState original;
        const auto initialStatus =
            original.updatePTY(pressure,temperature,composition);
        check("legacy dense PTY succeeds",initialStatus == Status::SUCCESS);
        if (initialStatus == Status::SUCCESS)
        {
            RealFluidState recovered;
            const auto recoveredStatus = recovered.updateRhoEY(
                original.rho(),original.e(),composition,1.05*temperature);
            check("legacy dense rhoE succeeds from a nearby branch seed",
                recoveredStatus == Status::SUCCESS,
                static_cast<double>(recoveredStatus),
                static_cast<double>(Status::SUCCESS));
            if (recoveredStatus == Status::SUCCESS)
            {
                check("legacy dense temperature round trip",
                    relativeError(recovered.T(),temperature) < 2.0e-8,
                    recovered.T(),temperature);
                check("legacy dense pressure round trip",
                    relativeError(recovered.p(),pressure) < 1.0e-7,
                    recovered.p(),pressure);
            }
        }
    }

    {
        RealFluidState oxygen;
        oxygen.setSpecies({"O2"});
        check("O2 pseudo-critical state succeeds",
            oxygen.updatePTY(5.043e6, 154.58, {1.0}) == Status::SUCCESS);
        check("O2 pseudo-critical Tr follows composition",
            std::abs(oxygen.Tr() - 1.0) < 1.0e-12, oxygen.Tr(), 1.0);
        check("O2 pseudo-critical Pr follows composition",
            std::abs(oxygen.Pr() - 1.0) < 1.0e-12, oxygen.Pr(), 1.0);

        RealFluidState propane;
        propane.setSpecies({"C3H8"});
        check("C3H8 pseudo-critical state succeeds",
            propane.updatePTY(4.251e6, 369.89, {1.0}) == Status::SUCCESS);
        check("C3H8 pseudo-critical Tr follows composition",
            std::abs(propane.Tr() - 1.0) < 1.0e-12, propane.Tr(), 1.0);
        check("C3H8 pseudo-critical Pr follows composition",
            std::abs(propane.Pr() - 1.0) < 1.0e-12, propane.Pr(), 1.0);
    }

    {
        RealFluidState state;
        state.setSpecies({"O2"});
        const auto status = state.updateRhoEY(1000.0, 1.0e12, {1.0}, 130.0);
        check("unbracketed rhoE target fails closed", status != Status::SUCCESS);
        check("unbracketed rhoE state is invalid", !state.isValid());
    }

    std::cout << (failures == 0
        ? "REAL_FLUID_STATE_CONTRACT_PASS\n"
        : "REAL_FLUID_STATE_CONTRACT_FAIL\n");
    return failures == 0 ? 0 : 1;
}
