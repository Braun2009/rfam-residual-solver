#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "PhaseEquilibrium/PRPHFlash.H"

using namespace Foam;
using namespace Foam::phaseEquilibrium;

namespace
{

const char* regionName(EquilibriumRegion region)
{
    switch (region)
    {
        case EquilibriumRegion::UNKNOWN: return "UNKNOWN";
        case EquilibriumRegion::SINGLE_LIQUID: return "SINGLE_LIQUID";
        case EquilibriumRegion::SINGLE_VAPOR: return "SINGLE_VAPOR";
        case EquilibriumRegion::SINGLE_SUPERCRITICAL:
            return "SINGLE_SUPERCRITICAL";
        case EquilibriumRegion::TWO_PHASE: return "TWO_PHASE";
        case EquilibriumRegion::PURE_COMPONENT_COEXISTENCE:
            return "PURE_COMPONENT_COEXISTENCE";
    }
    return "UNKNOWN";
}

ConservedPhaseProperties conserved(
    const PRPhaseCaloricProperties& phase)
{
    ConservedPhaseProperties result;
    result.molecularWeight = phase.molecularWeight;
    result.rho = phase.rho;
    result.h = phase.h;
    result.e = phase.e;
    result.s = phase.s;
    result.massFractions = phase.massFractions;
    return result;
}

std::optional<double> equilibriumEnthalpy(
    const PRTPFlash& flash,
    const PRPhaseCaloric& caloric,
    double pressure,
    double temperature,
    const std::vector<double>& massFractions)
{
    const TPFlashResult tp = flash.flashMassFractions(
        pressure, temperature, massFractions);
    if (!tp.valid()) return std::nullopt;

    if (tp.region == EquilibriumRegion::TWO_PHASE)
    {
        const PRPhaseCaloricResult liquid = caloric.evaluate(
            pressure, temperature, tp.liquidMoleFractions,
            PhaseRoot::LIQUID);
        const PRPhaseCaloricResult vapor = caloric.evaluate(
            pressure, temperature, tp.vaporMoleFractions,
            PhaseRoot::VAPOR);
        if (!liquid.success() || !vapor.success()) return std::nullopt;
        const ConservedMixtureResult mixture = assembleConservedPhaseMixture(
            pressure, tp.vaporMoleFraction,
            conserved(*liquid.state), conserved(*vapor.state),
            tp.feedMassFractions);
        return mixture.success()
            ? std::optional<double>(mixture.state->h) : std::nullopt;
    }

    if (!tp.singlePhaseState.valid()
        || tp.singlePhaseMoleFractions.empty())
        return std::nullopt;
    const PRPhaseCaloricResult single = caloric.evaluate(
        pressure, temperature, tp.singlePhaseMoleFractions,
        tp.singlePhaseState.requestedRoot);
    return single.success()
        ? std::optional<double>(single.state->h) : std::nullopt;
}

double vaporMassFraction(const PHFlashResult& result)
{
    if (!result.success()) return std::numeric_limits<double>::quiet_NaN();
    if (result.conservedMixture)
        return result.conservedMixture->vaporMassFraction;
    if (result.tpFlash.region == EquilibriumRegion::SINGLE_LIQUID)
        return 0.0;
    if (result.tpFlash.region == EquilibriumRegion::SINGLE_VAPOR
        || result.tpFlash.region == EquilibriumRegion::SINGLE_SUPERCRITICAL)
        return 1.0;
    return std::numeric_limits<double>::quiet_NaN();
}

} // namespace

int main(int argc, char** argv)
{
    std::size_t scanIntervals = 512;
    if (argc == 3 && std::string(argv[1]) == "--scan-intervals")
    {
        try
        {
            const unsigned long parsed = std::stoul(argv[2]);
            if (parsed == 0 || parsed > 8192) throw std::out_of_range("scan");
            scanIntervals = static_cast<std::size_t>(parsed);
        }
        catch (const std::exception&)
        {
            std::cerr << "HEM_INLET_MIXING_PATH_FAIL reason="
                      << "invalid-scan-intervals\n";
            return 2;
        }
    }
    else if (argc != 1)
    {
        std::cerr << "HEM_INLET_MIXING_PATH_FAIL reason=invalid-arguments\n";
        return 2;
    }

    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const std::vector<std::string> species = {"O2", "C3H8"};
    const std::vector<double> lox = {1.0, 0.0};
    const std::vector<double> propane = {0.0, 1.0};
    const std::vector<double> pressures = {4.0e6, 5.0e6, 6.0e6};
    const std::vector<double> kijValues = {-0.05, 0.0, 0.055, 0.10, 0.15};
    const std::vector<double> fuelMassFractions = {
        0.0, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25,
        1.05/3.85, 0.30, 0.40, 0.50, 0.70, 0.90, 1.0};

    std::size_t rows = 0;
    std::size_t successes = 0;
    std::size_t failures = 0;
    std::size_t twoPhase = 0;
    std::size_t unbracketed = 0;

    std::cout << std::setprecision(17);
    std::cout
        << "# rfam.hem-inlet-mixing-path-diagnostic/v1\n"
        << "pressure_pa,kij,fuel_mass_fraction,o2_mole_fraction,"
        << "target_h_j_kg,status,phase_status,temperature_k,region,"
        << "vapor_mass_fraction,enthalpy_relative_residual,"
        << "valid_temperature_samples,invalid_temperature_samples,"
        << "phase_topology_validated,production_ready\n";

    for (double pressure : pressures)
    {
        for (double kijValue : kijValues)
        {
            BinaryInteractionDB kij;
            kij.add("O2", "C3H8", kijValue);
            PRTPFlash tp(critical, kij, species);
            PRPhaseCaloric caloric(critical, kij, species);
            const std::optional<double> loxH = equilibriumEnthalpy(
                tp, caloric, pressure, 130.0, lox);
            const std::optional<double> propaneH = equilibriumEnthalpy(
                tp, caloric, pressure, 450.0, propane);
            if (!loxH || !propaneH)
            {
                std::cerr << "HEM_INLET_MIXING_PATH_FAIL reason="
                          << "endpoint-enthalpy pressure=" << pressure
                          << " kij=" << kijValue << '\n';
                return 1;
            }

            PHFlashOptions options;
            options.minimumTemperature = 80.0;
            options.maximumTemperature = 700.0;
            options.scanIntervals = scanIntervals;
            options.maximumIterations = 240;
            PRPHFlash ph(critical, kij, species, options);

            for (double fuelMassFraction : fuelMassFractions)
            {
                const std::vector<double> overall = {
                    1.0 - fuelMassFraction, fuelMassFraction};
                std::vector<double> moleFractions;
                const PhaseEquilibriumStatus conversion =
                    tp.kernel().massToMoleFractions(overall, moleFractions);
                if (conversion != PhaseEquilibriumStatus::SUCCESS)
                {
                    std::cerr << "HEM_INLET_MIXING_PATH_FAIL reason="
                              << "mass-to-mole status="
                              << statusName(conversion) << '\n';
                    return 1;
                }

                const double targetH = (1.0 - fuelMassFraction)*(*loxH)
                    + fuelMassFraction*(*propaneH);
                const PHFlashResult result = ph.flashMassFractions(
                    pressure, targetH, overall);
                ++rows;
                if (result.success())
                {
                    ++successes;
                    if (result.tpFlash.region == EquilibriumRegion::TWO_PHASE)
                        ++twoPhase;
                }
                else
                {
                    ++failures;
                    if (result.status == PHFlashStatus::NO_ENTHALPY_BRACKET)
                        ++unbracketed;
                }

                std::cout
                    << pressure << ',' << kijValue << ','
                    << fuelMassFraction << ',' << moleFractions[0] << ','
                    << targetH << ',' << phFlashStatusName(result.status)
                    << ',' << statusName(result.phaseStatus) << ','
                    << result.temperature << ','
                    << regionName(result.tpFlash.region) << ','
                    << vaporMassFraction(result) << ','
                    << result.enthalpyRelativeResidual << ','
                    << result.validTemperatureSamples << ','
                    << result.invalidTemperatureSamples
                    << ",0,0\n";
            }
        }
    }

    if (rows != pressures.size()*kijValues.size()*fuelMassFractions.size()
        || successes + failures != rows || twoPhase == 0 || failures == 0)
    {
        std::cerr << "HEM_INLET_MIXING_PATH_FAIL rows=" << rows
                  << " successes=" << successes
                  << " failures=" << failures
                  << " two_phase=" << twoPhase << '\n';
        return 1;
    }

    std::cout << "# HEM_INLET_MIXING_PATH_DIAGNOSTIC_BLOCKED rows="
              << rows << " successes=" << successes
              << " failures=" << failures
              << " unbracketed=" << unbracketed
              << " two_phase=" << twoPhase
              << " scan_intervals=" << scanIntervals
              << " phase_topology_validated=0 production_ready=0\n";
    return 0;
}
