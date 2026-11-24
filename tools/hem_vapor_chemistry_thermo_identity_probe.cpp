#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "MixingRules/BinaryInteractionDB.H"
#include "MixingRules/SpeciesCriticalDB.H"
#include "PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "RealFluidState.H"
#include "SpeciesRegistry.H"

using namespace Foam;
using namespace Foam::phaseEquilibrium;
using namespace Foam::thermo;

namespace
{

RhoEFlashOptions diagnosticOptions()
{
    RhoEFlashOptions options;
    options.minimumPressure = 1.0e5;
    options.maximumPressure = 1.0e7;
    options.minimumTemperature = 250.0;
    options.maximumTemperature = 2000.0;
    options.initialPressure = 4.0e6;
    options.initialTemperature = 800.0;
    options.pressureScanIntervals = 64;
    options.maximumIterations = 120;
    return options;
}

bool finitePositive(double value)
{
    return std::isfinite(value) && value > 0.0;
}

} // namespace

int main()
{
    const SpeciesRegistry registry;
    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();

    BinaryInteractionDB configuredKij = BinaryInteractionDB::createDefault();
    constexpr double configuredO2C3H8Kij = 0.20;
    configuredKij.add("O2", "C3H8", configuredO2C3H8Kij);

    PREquilibriumThermoModel configuredModel(
        critical, configuredKij, registry.names(), diagnosticOptions());

    std::vector<double> massFractions(
        static_cast<std::size_t>(registry.size()), 0.0);
    massFractions[registry.index("O2")] = 0.80;
    massFractions[registry.index("C3H8")] = 0.20;

    constexpr double targetPressure = 4.0e6;
    constexpr double targetTemperature = 800.0;
    const EquilibriumResult configured = configuredModel.initializePTY(
        targetPressure, targetTemperature, massFractions);
    if (!configured.success() || !configured.state.has_value())
    {
        std::cerr << "HEM_VAPOR_CHEMISTRY_THERMO_IDENTITY_PROBE_FAIL "
                  << "configured_model_initialization\n";
        return 1;
    }

    const EquilibriumState& authoritative = *configured.state;
    if (authoritative.region == thermo::EquilibriumRegion::TWO_PHASE
        || !finitePositive(authoritative.rho)
        || !std::isfinite(authoritative.e))
    {
        std::cerr << "HEM_VAPOR_CHEMISTRY_THERMO_IDENTITY_PROBE_FAIL "
                  << "diagnostic_state_not_single_phase\n";
        return 1;
    }

    RealFluidState cvodeThermo;
    cvodeThermo.setSpecies(registry.names());
    const RealFluidState::Status cvodeStatus = cvodeThermo.updateRhoEY(
        authoritative.rho, authoritative.e, massFractions,
        authoritative.T);
    if (cvodeStatus != RealFluidState::Status::SUCCESS)
    {
        std::cerr << "HEM_VAPOR_CHEMISTRY_THERMO_IDENTITY_PROBE_FAIL "
                  << "legacy_cvode_thermo_recovery status="
                  << static_cast<int>(cvodeStatus) << '\n';
        return 1;
    }

    const double pressureDifference = cvodeThermo.p() - authoritative.p;
    const double temperatureDifference = cvodeThermo.T() - authoritative.T;
    const double pressureRelativeDifference = std::abs(pressureDifference)
        /std::max(std::abs(authoritative.p), 1.0);
    const double temperatureRelativeDifference = std::abs(temperatureDifference)
        /std::max(std::abs(authoritative.T), 1.0);

    std::cout << std::setprecision(17)
              << "scope=offline-nonproduction-thermo-identity-diagnostic\n"
              << "configured_model_id=" << configuredModel.modelId() << '\n'
              << "configured_o2_c3h8_kij=" << configuredO2C3H8Kij << '\n'
              << "legacy_static_o2_c3h8_kij="
              << BinaryInteractionDB::createDefault().get("O2", "C3H8")
              << '\n'
              << "rho=" << authoritative.rho << '\n'
              << "specific_energy=" << authoritative.e << '\n'
              << "configured_pressure=" << authoritative.p << '\n'
              << "configured_temperature=" << authoritative.T << '\n'
              << "cvode_pressure=" << cvodeThermo.p() << '\n'
              << "cvode_temperature=" << cvodeThermo.T() << '\n'
              << "pressure_relative_difference="
              << pressureRelativeDifference << '\n'
              << "temperature_relative_difference="
              << temperatureRelativeDifference << '\n';

    constexpr double identityTolerance = 1.0e-8;
    if (pressureRelativeDifference <= identityTolerance
        && temperatureRelativeDifference <= identityTolerance)
    {
        std::cerr << "HEM_VAPOR_CHEMISTRY_THERMO_IDENTITY_PROBE_FAIL "
                  << "diagnostic_did_not_separate_backends\n";
        return 1;
    }

    std::cout << "HEM_VAPOR_CHEMISTRY_THERMO_IDENTITY_GAP_CONFIRMED\n";
    return 0;
}
