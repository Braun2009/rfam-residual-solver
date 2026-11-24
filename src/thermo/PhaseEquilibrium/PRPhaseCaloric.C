#include "PRPhaseCaloric.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include "NASA7Thermo.H"
#include "TabulatedIdealCaloric.H"

namespace Foam
{
namespace phaseEquilibrium
{

namespace
{

bool finitePositive(double value)
{
    return std::isfinite(value) && value > 0.0;
}

double relativeResidual(double lhs, double rhs)
{
    return std::abs(lhs - rhs)
        / std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

bool c1HandoffWindow(
    const std::string& species,
    double& minimum,
    double& maximum)
{
    if (species == "O2")
    {
        minimum = 195.0;
        maximum = 200.0;
        return true;
    }
    if (species == "C3H8")
    {
        minimum = 295.0;
        maximum = NASA7Thermo::sensibleReferenceTemperature;
        return true;
    }
    return false;
}

double smoothstep(double u)
{
    const double bounded = std::clamp(u, 0.0, 1.0);
    return bounded*bounded*(3.0 - 2.0*bounded);
}

double smoothstepIntegral(double u)
{
    const double bounded = std::clamp(u, 0.0, 1.0);
    const double squared = bounded*bounded;
    return squared*bounded - 0.5*squared*squared;
}

double smoothstepCpOverTemperatureIntegral(
    double minimum,
    double maximum,
    double u)
{
    const double bounded = std::clamp(u, 0.0, 1.0);
    const double ratio = (maximum - minimum)/minimum;
    double ratioPower = 1.0;
    double uPower = bounded*bounded*bounded;
    double sum = 0.0;
    for (int order = 0; order < 16; ++order)
    {
        const double moment =
            3.0*uPower/static_cast<double>(order + 3)
          - 2.0*uPower*bounded/static_cast<double>(order + 4);
        sum += (order % 2 == 0 ? 1.0 : -1.0)*ratioPower*moment;
        ratioPower *= ratio;
        uPower *= bounded;
    }
    return ratio*sum;
}

TabulatedIdealCaloricResult correctedLowTemperatureCaloric(
    const TabulatedIdealCaloric& table,
    const std::string& species,
    double temperature,
    double handoffMinimum,
    double handoffMaximum,
    double handoffCpCorrection,
    bool hasC1Handoff)
{
    TabulatedIdealCaloricResult result = table.evaluate(species, temperature);
    if (!result.success() || !hasC1Handoff) return result;

    const double width = handoffMaximum - handoffMinimum;
    const double u = (temperature - handoffMinimum)/width;
    result.state->cp += handoffCpCorrection*smoothstep(u);
    result.state->sensibleEnthalpy += handoffCpCorrection*width
        *smoothstepIntegral(u);
    result.state->sensibleEntropy += handoffCpCorrection
        *smoothstepCpOverTemperatureIntegral(
            handoffMinimum, handoffMaximum, u);
    return result;
}

double standardEntropyMass(
    const NASA7Thermo& nasa,
    const std::string& species,
    double temperature,
    double nasaMinimum,
    double molecularWeight)
{
    if (nasaMinimum > NASA7Thermo::sensibleReferenceTemperature
        && temperature >= NASA7Thermo::sensibleReferenceTemperature
        && temperature < nasaMinimum)
    {
        const double cpAtMinimum =
            nasa.cp_mole(species, nasaMinimum)/molecularWeight;
        return nasa.s_mole(species, nasaMinimum)/molecularWeight
            + cpAtMinimum*std::log(temperature/nasaMinimum);
    }
    return nasa.s_mole(species, temperature)/molecularWeight;
}

struct IdealCaloricResult
{
    bool valid = false;
    double cp = std::numeric_limits<double>::quiet_NaN();
    double sensibleEnthalpy = std::numeric_limits<double>::quiet_NaN();
    double standardEntropy = std::numeric_limits<double>::quiet_NaN();
};

IdealCaloricResult evaluateIdealCaloric(
    const NASA7Thermo& nasa,
    const TabulatedIdealCaloric& lowTemperature,
    const std::string& species,
    double temperature,
    double nasaMinimum,
    double nasaMaximum,
    double lowTemperatureOffset,
    double lowTemperatureEntropyOffset,
    bool hasLowTemperatureTable,
    double handoffMinimum,
    double handoffMaximum,
    double handoffCpCorrection,
    bool hasC1Handoff)
{
    IdealCaloricResult result;
    const double molecularWeight = nasa.MW(species);
    if (!finitePositive(molecularWeight)) return result;

    double massCp = std::numeric_limits<double>::quiet_NaN();
    double molarEnthalpy = std::numeric_limits<double>::quiet_NaN();
    double massEntropy = std::numeric_limits<double>::quiet_NaN();
    if (temperature >= nasaMinimum && temperature <= nasaMaximum)
    {
        massCp = nasa.cp_mole(species, temperature)/molecularWeight;
        molarEnthalpy = nasa.hsens_mole(species, temperature);
        massEntropy = standardEntropyMass(
            nasa, species, temperature, nasaMinimum, molecularWeight);
    }
    else if (temperature < nasaMinimum)
    {
        if (!hasLowTemperatureTable) return result;

        // For species whose NASA range starts above 298.15 K (C3H8), retain
        // the existing NASA linear bridge from the sensible reference to the
        // published lower bound. Below the reference, use the frozen ideal
        // table. For species whose NASA range starts below 298.15 K (O2), the
        // table is shifted once at the NASA lower bound so h remains continuous.
        if (nasaMinimum > NASA7Thermo::sensibleReferenceTemperature
            && temperature >= NASA7Thermo::sensibleReferenceTemperature)
        {
            massCp = nasa.cp_mole(species, temperature)/molecularWeight;
            molarEnthalpy = nasa.hsens_mole(species, temperature);
            massEntropy = standardEntropyMass(
                nasa, species, temperature, nasaMinimum, molecularWeight);
        }
        else
        {
            const TabulatedIdealCaloricResult table =
                correctedLowTemperatureCaloric(
                    lowTemperature,
                    species,
                    temperature,
                    handoffMinimum,
                    handoffMaximum,
                    handoffCpCorrection,
                    hasC1Handoff);
            if (!table.success()) return result;
            massCp = table.state->cp;
            molarEnthalpy = molecularWeight*table.state->sensibleEnthalpy
                + molecularWeight*lowTemperatureOffset;
            massEntropy = table.state->sensibleEntropy
                + lowTemperatureEntropyOffset;
        }
    }

    if (!finitePositive(massCp)
        || !std::isfinite(molarEnthalpy) || !std::isfinite(massEntropy))
        return result;
    result.cp = massCp;
    result.sensibleEnthalpy = molarEnthalpy/molecularWeight;
    result.standardEntropy = massEntropy;
    result.valid = finitePositive(result.cp)
        && std::isfinite(result.sensibleEnthalpy)
        && std::isfinite(result.standardEntropy);
    return result;
}

} // End anonymous namespace

const char* phaseCaloricStatusName(PRPhaseCaloricStatus status)
{
    switch (status)
    {
        case PRPhaseCaloricStatus::SUCCESS: return "SUCCESS";
        case PRPhaseCaloricStatus::INVALID_CONFIGURATION:
            return "INVALID_CONFIGURATION";
        case PRPhaseCaloricStatus::MISSING_CALORIC_DATA:
            return "MISSING_CALORIC_DATA";
        case PRPhaseCaloricStatus::MOLECULAR_WEIGHT_MISMATCH:
            return "MOLECULAR_WEIGHT_MISMATCH";
        case PRPhaseCaloricStatus::MISSING_BINARY_INTERACTION:
            return "MISSING_BINARY_INTERACTION";
        case PRPhaseCaloricStatus::INVALID_PRESSURE:
            return "INVALID_PRESSURE";
        case PRPhaseCaloricStatus::INVALID_TEMPERATURE:
            return "INVALID_TEMPERATURE";
        case PRPhaseCaloricStatus::INVALID_COMPOSITION:
            return "INVALID_COMPOSITION";
        case PRPhaseCaloricStatus::INVALID_PHASE_ROOT:
            return "INVALID_PHASE_ROOT";
        case PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE:
            return "UNSUPPORTED_CALORIC_RANGE";
        case PRPhaseCaloricStatus::PHASE_ROOT_FAILURE:
            return "PHASE_ROOT_FAILURE";
        case PRPhaseCaloricStatus::NONFINITE_PROPERTY:
            return "NONFINITE_PROPERTY";
        case PRPhaseCaloricStatus::EOS_CLOSURE_FAILURE:
            return "EOS_CLOSURE_FAILURE";
    }
    return "UNKNOWN_STATUS";
}

PRPhaseCaloric::PRPhaseCaloric(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    double compositionTolerance)
:
    criticalData_(criticalData),
    kij_(kij),
    speciesNames_(speciesNames),
    kernel_(criticalData_, kij_, speciesNames, compositionTolerance),
    eos_(criticalData_, kij_),
    configurationStatus_(PRPhaseCaloricStatus::INVALID_CONFIGURATION),
    configurationPhaseStatus_(kernel_.configurationStatus())
{
    if (!kernel_.configured()) return;

    caloricRanges_.reserve(speciesNames_.size());
    const NASA7Thermo& caloric = nasa7Thermo();
    const TabulatedIdealCaloric& lowTemperature =
        tabulatedIdealCaloric();
    for (std::size_t i = 0; i < speciesNames_.size(); ++i)
    {
        const std::string& species = speciesNames_[i];
        TemperatureRange range;
        if (!caloric.has(species) ||
            !caloric.temperatureRange(
                species, range.nasaMinimum, range.nasaMaximum))
        {
            configurationStatus_ =
                PRPhaseCaloricStatus::MISSING_CALORIC_DATA;
            configurationSpecies_ = species;
            return;
        }

        range.minimum = range.nasaMinimum;
        range.maximum = range.nasaMaximum;
        double tableMinimum = 0.0;
        double tableMaximum = 0.0;
        if (lowTemperature.temperatureRange(
                species, tableMinimum, tableMaximum))
        {
            range.hasLowTemperatureTable = true;
            range.minimum = std::min(range.minimum, tableMinimum);
            const double anchorTemperature = range.nasaMinimum
                <= NASA7Thermo::sensibleReferenceTemperature
                ? range.nasaMinimum
                : NASA7Thermo::sensibleReferenceTemperature;
            double handoffMinimum = 0.0;
            double handoffMaximum = 0.0;
            if (c1HandoffWindow(
                    species, handoffMinimum, handoffMaximum))
            {
                const TabulatedIdealCaloricResult rawBoundary =
                    lowTemperature.evaluate(species, handoffMaximum);
                const double targetCp = caloric.cp_mole(
                    species, range.nasaMinimum)/caloric.MW(species);
                if (!rawBoundary.success()
                    || anchorTemperature != handoffMaximum
                    || !finitePositive(targetCp))
                {
                    configurationStatus_ =
                        PRPhaseCaloricStatus::MISSING_CALORIC_DATA;
                    configurationSpecies_ = species;
                    return;
                }
                range.handoffMinimum = handoffMinimum;
                range.handoffMaximum = handoffMaximum;
                range.handoffCpCorrection =
                    targetCp - rawBoundary.state->cp;
                range.hasC1Handoff =
                    std::isfinite(range.handoffCpCorrection);
                if (!range.hasC1Handoff)
                {
                    configurationStatus_ =
                        PRPhaseCaloricStatus::MISSING_CALORIC_DATA;
                    configurationSpecies_ = species;
                    return;
                }
            }
            const TabulatedIdealCaloricResult boundary =
                correctedLowTemperatureCaloric(
                    lowTemperature,
                    species,
                    anchorTemperature,
                    range.handoffMinimum,
                    range.handoffMaximum,
                    range.handoffCpCorrection,
                    range.hasC1Handoff);
            if (!boundary.success())
            {
                configurationStatus_ =
                    PRPhaseCaloricStatus::MISSING_CALORIC_DATA;
                configurationSpecies_ = species;
                return;
            }
            const double nasaBoundary = caloric.hsens_mole(
                species, anchorTemperature)/caloric.MW(species);
            const double nasaEntropyBoundary = standardEntropyMass(
                caloric,
                species,
                anchorTemperature,
                range.nasaMinimum,
                caloric.MW(species));
            range.lowTemperatureOffset =
                nasaBoundary - boundary.state->sensibleEnthalpy;
            range.lowTemperatureEntropyOffset =
                nasaEntropyBoundary - boundary.state->sensibleEntropy;
        }

        const double criticalMW = kernel_.properties()[i].MW;
        const double caloricMW = caloric.MW(species);
        if (!finitePositive(criticalMW) || !finitePositive(caloricMW) ||
            relativeResidual(criticalMW, caloricMW) > 1.0e-12)
        {
            configurationStatus_ =
                PRPhaseCaloricStatus::MOLECULAR_WEIGHT_MISMATCH;
            configurationSpecies_ = species;
            return;
        }
        caloricRanges_.push_back(range);
    }

    try
    {
        eos_.setSpeciesNames(speciesNames_);
    }
    catch (...)
    {
        configurationStatus_ =
            PRPhaseCaloricStatus::INVALID_CONFIGURATION;
        return;
    }

    configurationStatus_ = PRPhaseCaloricStatus::SUCCESS;
    configurationPhaseStatus_ = PhaseEquilibriumStatus::SUCCESS;
}

PRPhaseCaloricTemperatureRange PRPhaseCaloric::supportedTemperatureRange(
    const std::vector<double>& fractions) const
{
    PRPhaseCaloricTemperatureRange result;
    result.status = configurationStatus_;
    result.phaseStatus = configurationPhaseStatus_;
    result.limitingSpecies = configurationSpecies_;
    if (!configured()) return result;

    result.phaseStatus = kernel_.validateMoleFractions(fractions);
    if (result.phaseStatus != PhaseEquilibriumStatus::SUCCESS)
    {
        result.status = PRPhaseCaloricStatus::INVALID_COMPOSITION;
        return result;
    }

    result.minimum = 0.0;
    result.maximum = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < fractions.size(); ++i)
    {
        if (fractions[i] == 0.0) continue;
        const TemperatureRange& range = caloricRanges_[i];
        result.minimum = std::max(result.minimum, range.minimum);
        result.maximum = std::min(result.maximum, range.maximum);
        if (result.minimum > result.maximum)
        {
            result.status =
                PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE;
            result.limitingSpecies = speciesNames_[i];
            return result;
        }
    }

    if (!finitePositive(result.minimum)
        || !finitePositive(result.maximum))
    {
        result.status = PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE;
        return result;
    }

    result.status = PRPhaseCaloricStatus::SUCCESS;
    result.phaseStatus = PhaseEquilibriumStatus::SUCCESS;
    result.limitingSpecies.clear();
    return result;
}

std::optional<double> PRPhaseCaloric::idealCpMass(
    double temperature,
    const std::vector<double>& massFractions) const
{
    if (!configured() || !finitePositive(temperature)
        || kernel_.validateMoleFractions(massFractions)
            != PhaseEquilibriumStatus::SUCCESS)
        return std::nullopt;

    const NASA7Thermo& caloric = nasa7Thermo();
    const TabulatedIdealCaloric& lowTemperature = tabulatedIdealCaloric();
    double mixtureCp = 0.0;
    for (std::size_t i = 0; i < massFractions.size(); ++i)
    {
        if (massFractions[i] == 0.0) continue;
        const TemperatureRange& range = caloricRanges_[i];
        if (temperature < range.minimum || temperature > range.maximum)
            return std::nullopt;
        const IdealCaloricResult speciesCaloric = evaluateIdealCaloric(
            caloric,
            lowTemperature,
            speciesNames_[i],
            temperature,
            range.nasaMinimum,
            range.nasaMaximum,
            range.lowTemperatureOffset,
            range.lowTemperatureEntropyOffset,
            range.hasLowTemperatureTable,
            range.handoffMinimum,
            range.handoffMaximum,
            range.handoffCpCorrection,
            range.hasC1Handoff);
        if (!speciesCaloric.valid || !finitePositive(speciesCaloric.cp))
            return std::nullopt;
        mixtureCp += massFractions[i]*speciesCaloric.cp;
    }
    return finitePositive(mixtureCp)
        ? std::optional<double>(mixtureCp) : std::nullopt;
}

PRPhaseCaloricResult PRPhaseCaloric::evaluate(
    double pressure,
    double temperature,
    const std::vector<double>& moleFractions,
    PhaseRoot root) const
{
    PRPhaseCaloricResult result;
    result.status = configurationStatus_;
    result.phaseStatus = configurationPhaseStatus_;
    result.limitingSpecies = configurationSpecies_;
    if (!configured()) return result;

    if (!finitePositive(pressure))
    {
        result.status = PRPhaseCaloricStatus::INVALID_PRESSURE;
        result.phaseStatus = PhaseEquilibriumStatus::INVALID_PRESSURE;
        return result;
    }
    if (!finitePositive(temperature))
    {
        result.status = PRPhaseCaloricStatus::INVALID_TEMPERATURE;
        result.phaseStatus = PhaseEquilibriumStatus::INVALID_TEMPERATURE;
        return result;
    }
    if (root != PhaseRoot::LIQUID && root != PhaseRoot::VAPOR)
    {
        result.status = PRPhaseCaloricStatus::INVALID_PHASE_ROOT;
        return result;
    }

    result.phaseStatus = kernel_.validateMoleFractions(moleFractions);
    if (result.phaseStatus != PhaseEquilibriumStatus::SUCCESS)
    {
        result.status = PRPhaseCaloricStatus::INVALID_COMPOSITION;
        return result;
    }

    for (std::size_t i = 0; i < moleFractions.size(); ++i)
    {
        if (moleFractions[i] == 0.0) continue;
        const TemperatureRange& range = caloricRanges_[i];
        if (temperature < range.minimum || temperature > range.maximum)
        {
            result.status =
                PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE;
            result.limitingSpecies = speciesNames_[i];
            result.supportedTemperatureMinimum = range.minimum;
            result.supportedTemperatureMaximum = range.maximum;
            return result;
        }
    }

    const PRPhaseState rootState =
        kernel_.phaseState(pressure, temperature, moleFractions, root);
    result.phaseStatus = rootState.status;
    if (!rootState.valid())
    {
        result.status = rootState.status ==
            PhaseEquilibriumStatus::MISSING_BINARY_INTERACTION
            ? PRPhaseCaloricStatus::MISSING_BINARY_INTERACTION
            : PRPhaseCaloricStatus::PHASE_ROOT_FAILURE;
        return result;
    }

    std::vector<double> massFractions;
    result.phaseStatus =
        kernel_.moleToMassFractions(moleFractions, massFractions);
    if (result.phaseStatus != PhaseEquilibriumStatus::SUCCESS)
    {
        result.status = PRPhaseCaloricStatus::INVALID_COMPOSITION;
        return result;
    }

    PRPhaseCaloricProperties candidate;
    candidate.root = root;
    candidate.moleFractions = moleFractions;
    candidate.massFractions = massFractions;
    candidate.pressure = pressure;
    candidate.temperature = temperature;
    candidate.molecularWeight =
        kernel_.mixtureMolecularWeight(moleFractions);
    candidate.compressibility = rootState.compressibility;
    candidate.rho = pressure*candidate.molecularWeight
        /(candidate.compressibility
            *PRMixtureKernel::universalGasConstant*temperature);

    const NASA7Thermo& caloric = nasa7Thermo();
    const TabulatedIdealCaloric& lowTemperature = tabulatedIdealCaloric();
    double idealSensibleEnthalpy = 0.0;
    double idealEntropy = 0.0;
    static constexpr double standardPressure = 1.0e5;
    for (std::size_t i = 0; i < massFractions.size(); ++i)
    {
        if (massFractions[i] == 0.0) continue;
        const TemperatureRange& range = caloricRanges_[i];
        const IdealCaloricResult speciesCaloric = evaluateIdealCaloric(
            caloric,
            lowTemperature,
            speciesNames_[i],
            temperature,
            range.nasaMinimum,
            range.nasaMaximum,
            range.lowTemperatureOffset,
            range.lowTemperatureEntropyOffset,
            range.hasLowTemperatureTable,
            range.handoffMinimum,
            range.handoffMaximum,
            range.handoffCpCorrection,
            range.hasC1Handoff);
        if (!speciesCaloric.valid)
        {
            result.status = PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE;
            result.limitingSpecies = speciesNames_[i];
            result.supportedTemperatureMinimum = range.minimum;
            result.supportedTemperatureMaximum = range.maximum;
            return result;
        }
        idealSensibleEnthalpy += massFractions[i]
            *speciesCaloric.sensibleEnthalpy;
        const double partialPressureRatio =
            moleFractions[i]*pressure/standardPressure;
        if (!finitePositive(partialPressureRatio))
        {
            result.status = PRPhaseCaloricStatus::NONFINITE_PROPERTY;
            return result;
        }
        idealEntropy += massFractions[i]
            *(speciesCaloric.standardEntropy
                - PRMixtureKernel::universalGasConstant
                    /caloric.MW(speciesNames_[i])
                    *std::log(partialPressureRatio));
    }
    const double departureEnthalpy = eos_.enthalpyDeparture(
        pressure, temperature, candidate.compressibility, massFractions);
    const double departureEntropy = eos_.entropyDeparture(
        pressure, temperature, candidate.compressibility, massFractions);
    candidate.h = idealSensibleEnthalpy + departureEnthalpy;
    candidate.e = candidate.h - pressure/candidate.rho;
    candidate.s = idealEntropy + departureEntropy;

    if (!finitePositive(candidate.molecularWeight) ||
        !finitePositive(candidate.compressibility) ||
        !finitePositive(candidate.rho) ||
        !std::isfinite(idealSensibleEnthalpy) ||
        !std::isfinite(idealEntropy) ||
        !std::isfinite(departureEnthalpy) ||
        !std::isfinite(departureEntropy) ||
        !std::isfinite(candidate.h) || !std::isfinite(candidate.e) ||
        !std::isfinite(candidate.s))
    {
        result.status = PRPhaseCaloricStatus::NONFINITE_PROPERTY;
        return result;
    }

    const double recoveredPressure = eos_.pressure(
        candidate.rho, temperature, massFractions);
    result.pressureRelativeResidual =
        relativeResidual(recoveredPressure, pressure);
    result.energyIdentityRelativeResidual = relativeResidual(
        candidate.e, candidate.h - pressure/candidate.rho);
    if (!std::isfinite(recoveredPressure) ||
        !std::isfinite(result.pressureRelativeResidual) ||
        !std::isfinite(result.energyIdentityRelativeResidual) ||
        result.pressureRelativeResidual > 1.0e-9 ||
        result.energyIdentityRelativeResidual > 1.0e-13)
    {
        result.status = PRPhaseCaloricStatus::EOS_CLOSURE_FAILURE;
        return result;
    }

    result.status = PRPhaseCaloricStatus::SUCCESS;
    result.phaseStatus = PhaseEquilibriumStatus::SUCCESS;
    result.limitingSpecies.clear();
    result.state = std::move(candidate);
    return result;
}

} // End namespace phaseEquilibrium
} // End namespace Foam
