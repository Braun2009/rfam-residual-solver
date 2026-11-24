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

const char* conservedStatusName(ConservedMixtureStatus status)
{
    switch (status)
    {
        case ConservedMixtureStatus::SUCCESS: return "SUCCESS";
        case ConservedMixtureStatus::INVALID_INPUT: return "INVALID_INPUT";
        case ConservedMixtureStatus::INVALID_PHASE_COMPOSITION:
            return "INVALID_PHASE_COMPOSITION";
        case ConservedMixtureStatus::PHASE_ENERGY_IDENTITY_FAILED:
            return "PHASE_ENERGY_IDENTITY_FAILED";
        case ConservedMixtureStatus::MATERIAL_BALANCE_FAILED:
            return "MATERIAL_BALANCE_FAILED";
        case ConservedMixtureStatus::MIXTURE_ENERGY_IDENTITY_FAILED:
            return "MIXTURE_ENERGY_IDENTITY_FAILED";
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

struct EnthalpyEvaluation
{
    TPFlashResult tp;
    PRPhaseCaloricStatus caloricStatus =
        PRPhaseCaloricStatus::INVALID_CONFIGURATION;
    ConservedMixtureStatus conservedStatus =
        ConservedMixtureStatus::INVALID_INPUT;
    double enthalpy = std::numeric_limits<double>::quiet_NaN();

    bool valid() const
    {
        return tp.valid()
            && caloricStatus == PRPhaseCaloricStatus::SUCCESS
            && conservedStatus == ConservedMixtureStatus::SUCCESS
            && std::isfinite(enthalpy);
    }
};

EnthalpyEvaluation evaluate(
    const PRTPFlash& flash,
    const PRPhaseCaloric& caloric,
    double pressure,
    double temperature,
    const std::vector<double>& massFractions)
{
    EnthalpyEvaluation result;
    result.tp = flash.flashMassFractions(
        pressure, temperature, massFractions);
    if (!result.tp.valid()) return result;

    if (result.tp.region == EquilibriumRegion::TWO_PHASE)
    {
        const PRPhaseCaloricResult liquid = caloric.evaluate(
            pressure,
            temperature,
            result.tp.liquidMoleFractions,
            PhaseRoot::LIQUID);
        const PRPhaseCaloricResult vapor = caloric.evaluate(
            pressure,
            temperature,
            result.tp.vaporMoleFractions,
            PhaseRoot::VAPOR);
        result.caloricStatus = !liquid.success()
            ? liquid.status : vapor.status;
        if (!liquid.success() || !vapor.success()) return result;

        const ConservedMixtureResult mixture =
            assembleConservedPhaseMixture(
                pressure,
                result.tp.vaporMoleFraction,
                conserved(*liquid.state),
                conserved(*vapor.state),
                result.tp.feedMassFractions);
        result.conservedStatus = mixture.status;
        if (!mixture.success()) return result;
        result.enthalpy = mixture.state->h;
        result.caloricStatus = PRPhaseCaloricStatus::SUCCESS;
        return result;
    }

    if (!result.tp.singlePhaseState.valid()
        || result.tp.singlePhaseMoleFractions.empty())
    {
        return result;
    }
    const PRPhaseCaloricResult single = caloric.evaluate(
        pressure,
        temperature,
        result.tp.singlePhaseMoleFractions,
        result.tp.singlePhaseState.requestedRoot);
    result.caloricStatus = single.status;
    if (!single.success()) return result;
    result.enthalpy = single.state->h;
    result.conservedStatus = ConservedMixtureStatus::SUCCESS;
    return result;
}

double parseFinite(const char* text, bool positive)
{
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != std::string(text).size()
        || !std::isfinite(value) || (positive && value <= 0.0))
    {
        throw std::invalid_argument("finite-number");
    }
    return value;
}

std::size_t parseIntervals(const char* text)
{
    std::size_t consumed = 0;
    const unsigned long value = std::stoul(text, &consumed);
    if (consumed != std::string(text).size()
        || value < 16 || value > 20000)
    {
        throw std::invalid_argument("intervals");
    }
    return static_cast<std::size_t>(value);
}

} // namespace

int main(int argc, char** argv)
{
    double pressure = 6.0e6;
    double kijValue = 0.0;
    double fuelMassFraction = 0.10;
    std::size_t intervals = 2048;
    if (argc == 6 && std::string(argv[1]) == "--state")
    {
        try
        {
            pressure = parseFinite(argv[2], true);
            kijValue = parseFinite(argv[3], false);
            fuelMassFraction = parseFinite(argv[4], false);
            intervals = parseIntervals(argv[5]);
            if (fuelMassFraction < 0.0 || fuelMassFraction > 1.0)
                throw std::invalid_argument("fuel-fraction");
        }
        catch (const std::exception&)
        {
            std::cerr << "HEM_PH_TEMPERATURE_SCAN_FAIL reason=invalid-state\n";
            return 2;
        }
    }
    else if (argc != 1)
    {
        std::cerr << "HEM_PH_TEMPERATURE_SCAN_FAIL reason=invalid-arguments\n";
        return 2;
    }

    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    BinaryInteractionDB kij;
    kij.add("O2", "C3H8", kijValue);
    const std::vector<std::string> species = {"O2", "C3H8"};
    PRTPFlash flash(critical, kij, species);
    PRPhaseCaloric caloric(critical, kij, species);
    const std::vector<double> oxygen = {1.0, 0.0};
    const std::vector<double> propane = {0.0, 1.0};
    const EnthalpyEvaluation oxygenInlet = evaluate(
        flash, caloric, pressure, 130.0, oxygen);
    const EnthalpyEvaluation propaneInlet = evaluate(
        flash, caloric, pressure, 450.0, propane);
    if (!oxygenInlet.valid() || !propaneInlet.valid())
    {
        std::cerr << "HEM_PH_TEMPERATURE_SCAN_FAIL reason=inlet-enthalpy\n";
        return 1;
    }

    const std::vector<double> overall = {
        1.0 - fuelMassFraction, fuelMassFraction};
    const double targetEnthalpy =
        (1.0 - fuelMassFraction)*oxygenInlet.enthalpy
        + fuelMassFraction*propaneInlet.enthalpy;

    PHFlashOptions phOptions;
    phOptions.minimumTemperature = 80.0;
    phOptions.maximumTemperature = 700.0;
    phOptions.scanIntervals = intervals;
    phOptions.maximumIterations = 240;
    const PRPHFlash ph(critical, kij, species, phOptions);
    const PHFlashResult phResult = ph.flashMassFractions(
        pressure, targetEnthalpy, overall);

    std::cout << std::setprecision(17);
    std::cout << "# rfam.hem-ph-temperature-scan/v1\n"
              << "temperature_k,phase_status,caloric_status,"
              << "conserved_status,region,enthalpy_j_kg,residual_j_kg,"
              << "flash_iterations,fugacity_residual,material_residual\n";
    std::size_t valid = 0;
    std::size_t invalid = 0;
    for (std::size_t sample = 0; sample <= intervals; ++sample)
    {
        const double fraction = static_cast<double>(sample)
            /static_cast<double>(intervals);
        const double temperature = sample == intervals
            ? 700.0 : 80.0 + fraction*(700.0 - 80.0);
        const EnthalpyEvaluation state = evaluate(
            flash, caloric, pressure, temperature, overall);
        if (state.valid()) ++valid;
        else ++invalid;
        std::cout << temperature << ',' << statusName(state.tp.status)
                  << ',' << phaseCaloricStatusName(state.caloricStatus)
                  << ',' << conservedStatusName(state.conservedStatus)
                  << ',' << regionName(state.tp.region) << ','
                  << state.enthalpy << ','
                  << state.enthalpy - targetEnthalpy << ','
                  << state.tp.flashIterations << ','
                  << state.tp.maximumLogFugacityResidual << ','
                  << state.tp.maximumMoleMaterialResidual << '\n';
    }
    std::cout << "# HEM_PH_TEMPERATURE_SCAN_BLOCKED pressure_pa="
              << pressure << " kij=" << kijValue
              << " fuel_mass_fraction=" << fuelMassFraction
              << " target_h_j_kg=" << targetEnthalpy
              << " valid=" << valid << " invalid=" << invalid
              << " ph_status=" << phFlashStatusName(phResult.status)
              << " ph_temperature_k=" << phResult.temperature
              << " production_ready=0\n";
    return 0;
}
