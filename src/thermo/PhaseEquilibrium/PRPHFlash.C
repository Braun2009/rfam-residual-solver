#include "PRPHFlash.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

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

bool finiteNonNegative(double value)
{
    return std::isfinite(value) && value >= 0.0;
}

double relativeResidual(double lhs, double rhs)
{
    return std::abs(lhs - rhs)
        / std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

bool enthalpyConverged(
    double computed,
    double target,
    const PHFlashOptions& options)
{
    if (!std::isfinite(computed) || !std::isfinite(target)) return false;
    const double tolerance = options.enthalpyAbsoluteTolerance
        + options.enthalpyRelativeTolerance
            * std::max({std::abs(computed), std::abs(target), 1.0});
    return std::abs(computed - target) <= tolerance;
}

double temperatureTolerance(
    double lower,
    double upper,
    const PHFlashOptions& options)
{
    return options.temperatureAbsoluteTolerance
        + options.temperatureRelativeTolerance
            * std::max({std::abs(lower), std::abs(upper), 1.0});
}

ConservedPhaseProperties conservedProperties(
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

bool verifyTPClosure(
    const TPFlashResult& tp,
    const TPFlashOptions& options,
    double& fugacityResidual,
    double& materialResidual)
{
    fugacityResidual = 0.0;
    materialResidual = 0.0;
    if (!tp.valid() || tp.feedMoleFractions.empty()) return false;

    const std::vector<double>& feed = tp.feedMoleFractions;
    if (tp.region == EquilibriumRegion::TWO_PHASE)
    {
        if (!tp.phaseFractionDetermined
            || !std::isfinite(tp.vaporMoleFraction)
            || tp.vaporMoleFraction < 0.0
            || tp.vaporMoleFraction > 1.0
            || tp.liquidMoleFractions.size() != feed.size()
            || tp.vaporMoleFractions.size() != feed.size()
            || !tp.liquidState.valid() || !tp.vaporState.valid()
            || tp.liquidState.lnFugacityCoefficients.size() != feed.size()
            || tp.vaporState.lnFugacityCoefficients.size() != feed.size())
            return false;

        for (std::size_t species = 0; species < feed.size(); ++species)
        {
            const double liquid = tp.liquidMoleFractions[species];
            const double vapor = tp.vaporMoleFractions[species];
            if (!std::isfinite(liquid) || !std::isfinite(vapor)
                || liquid < 0.0 || vapor < 0.0)
                return false;
            const double reconstructed =
                (1.0 - tp.vaporMoleFraction)*liquid
                + tp.vaporMoleFraction*vapor;
            materialResidual = std::max(
                materialResidual, std::abs(feed[species] - reconstructed));

            if (feed[species] > 0.0)
            {
                if (!finitePositive(liquid) || !finitePositive(vapor))
                    return false;
                const double residual = std::log(liquid)
                    + tp.liquidState.lnFugacityCoefficients[species]
                    - std::log(vapor)
                    - tp.vaporState.lnFugacityCoefficients[species];
                if (!std::isfinite(residual)) return false;
                fugacityResidual = std::max(
                    fugacityResidual, std::abs(residual));
            }
        }
        if (!std::isfinite(tp.maximumLogFugacityResidual)
            || !std::isfinite(tp.maximumMoleMaterialResidual)
            || !std::isfinite(tp.rachfordRiceResidual)
            || tp.maximumLogFugacityResidual > options.fugacityTolerance
            || tp.maximumMoleMaterialResidual > options.materialTolerance
            || std::abs(tp.rachfordRiceResidual) > options.materialTolerance)
            return false;
    }
    else if (tp.region == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE)
    {
        if (feed.size() != 1 || tp.phaseFractionDetermined
            || !std::isnan(tp.vaporMoleFraction)
            || !tp.liquidState.valid() || !tp.vaporState.valid()
            || tp.liquidState.lnFugacityCoefficients.size() != 1
            || tp.vaporState.lnFugacityCoefficients.size() != 1)
            return false;
        fugacityResidual = std::abs(
            tp.liquidState.lnFugacityCoefficients[0]
            - tp.vaporState.lnFugacityCoefficients[0]);
        if (!std::isfinite(tp.maximumLogFugacityResidual)
            || !std::isfinite(tp.maximumMoleMaterialResidual)
            || tp.maximumLogFugacityResidual > options.fugacityTolerance
            || tp.maximumMoleMaterialResidual > options.materialTolerance)
            return false;
    }
    else
    {
        if (!tp.singlePhaseState.valid()
            || tp.singlePhaseMoleFractions.size() != feed.size())
            return false;
        for (std::size_t species = 0; species < feed.size(); ++species)
        {
            const double value = tp.singlePhaseMoleFractions[species];
            if (!std::isfinite(value) || value < 0.0) return false;
            materialResidual = std::max(
                materialResidual, std::abs(value - feed[species]));
        }
        if (!std::isfinite(tp.maximumMoleMaterialResidual)
            || tp.maximumMoleMaterialResidual > options.materialTolerance)
            return false;
    }

    return std::isfinite(fugacityResidual)
        && std::isfinite(materialResidual)
        && fugacityResidual <= options.fugacityTolerance
        && materialResidual <= options.materialTolerance;
}

} // End anonymous namespace

const char* phFlashStatusName(PHFlashStatus status)
{
    switch (status)
    {
        case PHFlashStatus::SUCCESS: return "SUCCESS";
        case PHFlashStatus::INVALID_CONFIGURATION:
            return "INVALID_CONFIGURATION";
        case PHFlashStatus::INVALID_OPTIONS: return "INVALID_OPTIONS";
        case PHFlashStatus::INVALID_PRESSURE: return "INVALID_PRESSURE";
        case PHFlashStatus::NONFINITE_TARGET_ENTHALPY:
            return "NONFINITE_TARGET_ENTHALPY";
        case PHFlashStatus::INVALID_COMPOSITION:
            return "INVALID_COMPOSITION";
        case PHFlashStatus::INVALID_TEMPERATURE_INTERVAL:
            return "INVALID_TEMPERATURE_INTERVAL";
        case PHFlashStatus::TP_FLASH_FAILURE: return "TP_FLASH_FAILURE";
        case PHFlashStatus::TP_CLOSURE_FAILURE:
            return "TP_CLOSURE_FAILURE";
        case PHFlashStatus::UNSUPPORTED_CALORIC_RANGE:
            return "UNSUPPORTED_CALORIC_RANGE";
        case PHFlashStatus::CALORIC_FAILURE: return "CALORIC_FAILURE";
        case PHFlashStatus::NONPOSITIVE_LATENT_ENTHALPY:
            return "NONPOSITIVE_LATENT_ENTHALPY";
        case PHFlashStatus::CONSERVED_MIXTURE_FAILURE:
            return "CONSERVED_MIXTURE_FAILURE";
        case PHFlashStatus::NO_VALID_CALORIC_INTERVAL:
            return "NO_VALID_CALORIC_INTERVAL";
        case PHFlashStatus::NO_ENTHALPY_BRACKET:
            return "NO_ENTHALPY_BRACKET";
        case PHFlashStatus::NOT_CONVERGED: return "NOT_CONVERGED";
        case PHFlashStatus::ENTHALPY_RESIDUAL_FAILURE:
            return "ENTHALPY_RESIDUAL_FAILURE";
    }
    return "UNKNOWN_STATUS";
}

PRPHFlash::PRPHFlash(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    const PHFlashOptions& options,
    const TPFlashOptions& tpOptions)
:
    speciesNames_(speciesNames),
    tpFlash_(criticalData, kij, speciesNames, tpOptions),
    caloric_(criticalData, kij, speciesNames),
    options_(options),
    configurationStatus_(PHFlashStatus::SUCCESS)
{
    if (!finitePositive(options_.minimumTemperature)
        || !finitePositive(options_.maximumTemperature)
        || options_.minimumTemperature >= options_.maximumTemperature)
    {
        configurationStatus_ =
            PHFlashStatus::INVALID_TEMPERATURE_INTERVAL;
        return;
    }
    if (options_.scanIntervals == 0 || options_.maximumIterations == 0
        || !finiteNonNegative(options_.temperatureAbsoluteTolerance)
        || !finiteNonNegative(options_.temperatureRelativeTolerance)
        || !finiteNonNegative(options_.enthalpyAbsoluteTolerance)
        || !finiteNonNegative(options_.enthalpyRelativeTolerance)
        || !finiteNonNegative(options_.rootRelativeTolerance)
        || !finiteNonNegative(
            options_.conservedTolerances.compositionAbsolute)
        || !finiteNonNegative(
            options_.conservedTolerances.energyRelative))
    {
        configurationStatus_ = PHFlashStatus::INVALID_OPTIONS;
        return;
    }
    if (tpFlash_.configurationStatus()
            != PhaseEquilibriumStatus::SUCCESS
        || !caloric_.configured())
    {
        configurationStatus_ = PHFlashStatus::INVALID_CONFIGURATION;
    }
}

PRPHFlash::TemperatureEvaluation PRPHFlash::evaluateTemperature(
    double pressure,
    double targetEnthalpy,
    const std::vector<double>& massFractions,
    double temperature) const
{
    TemperatureEvaluation result;
    result.temperature = temperature;
    result.tp = tpFlash_.flashMassFractions(
        pressure, temperature, massFractions);
    result.phaseStatus = result.tp.status;
    if (!result.tp.valid())
    {
        result.status = PHFlashStatus::TP_FLASH_FAILURE;
        return result;
    }

    if (!verifyTPClosure(
            result.tp,
            tpFlash_.options(),
            result.tpFugacityResidual,
            result.tpMaterialResidual))
    {
        result.status = PHFlashStatus::TP_CLOSURE_FAILURE;
        return result;
    }

    if (result.tp.region == EquilibriumRegion::TWO_PHASE
        || result.tp.region
            == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE)
    {
        const PRPhaseCaloricResult liquid = caloric_.evaluate(
            pressure,
            temperature,
            result.tp.liquidMoleFractions,
            PhaseRoot::LIQUID);
        const PRPhaseCaloricResult vapor = caloric_.evaluate(
            pressure,
            temperature,
            result.tp.vaporMoleFractions,
            PhaseRoot::VAPOR);
        result.caloricStatus = !liquid.success()
            ? liquid.status : vapor.status;
        result.phaseStatus = !liquid.success()
            ? liquid.phaseStatus : vapor.phaseStatus;
        if (!liquid.success() || !vapor.success())
        {
            result.status =
                result.caloricStatus
                    == PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE
                ? PHFlashStatus::UNSUPPORTED_CALORIC_RANGE
                : PHFlashStatus::CALORIC_FAILURE;
            return result;
        }

        result.liquid = liquid.state;
        result.vapor = vapor.state;
        result.phaseEnergyIdentityResidual = std::max(
            liquid.energyIdentityRelativeResidual,
            vapor.energyIdentityRelativeResidual);

        double vaporMoleFraction = result.tp.vaporMoleFraction;
        if (result.tp.region
            == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE)
        {
            const double liquidH = result.liquid->h;
            const double vaporH = result.vapor->h;
            const double latentH = vaporH - liquidH;
            const double latentTolerance =
                options_.enthalpyAbsoluteTolerance
                + options_.enthalpyRelativeTolerance
                    *std::max({std::abs(liquidH), std::abs(vaporH), 1.0});
            if (!std::isfinite(latentH) || latentH <= latentTolerance)
            {
                result.status =
                    PHFlashStatus::NONPOSITIVE_LATENT_ENTHALPY;
                return result;
            }
            if (targetEnthalpy < liquidH || targetEnthalpy > vaporH)
            {
                // A pure TP coexistence point is a real discontinuity in the
                // single-phase h(T) graph. An out-of-lever target breaks the
                // outer bracket; it is never bridged by interpolation.
                result.status = PHFlashStatus::NO_ENTHALPY_BRACKET;
                return result;
            }
            const double vaporMassFraction =
                (targetEnthalpy - liquidH)/latentH;
            if (!std::isfinite(vaporMassFraction)
                || vaporMassFraction < 0.0 || vaporMassFraction > 1.0
                || relativeResidual(
                    result.liquid->molecularWeight,
                    result.vapor->molecularWeight)
                    > options_.rootRelativeTolerance)
            {
                result.status =
                    PHFlashStatus::CONSERVED_MIXTURE_FAILURE;
                return result;
            }
            // For one species both phases have the same molecular weight, so
            // mass and mole vapor fractions are identical.
            vaporMoleFraction = vaporMassFraction;
            result.pureLeverRule = true;
        }

        const ConservedMixtureResult mixture =
            assembleConservedPhaseMixture(
                pressure,
                vaporMoleFraction,
                conservedProperties(*result.liquid),
                conservedProperties(*result.vapor),
                massFractions,
                options_.conservedTolerances);
        result.conservedStatus = mixture.status;
        result.mixtureMaterialResidual = mixture.materialResidual;
        result.mixtureEnergyIdentityResidual =
            mixture.mixtureEnergyIdentityResidual;
        result.phaseEnergyIdentityResidual = std::max(
            result.phaseEnergyIdentityResidual,
            mixture.phaseEnergyIdentityResidual);
        if (!mixture.success())
        {
            result.status = PHFlashStatus::CONSERVED_MIXTURE_FAILURE;
            return result;
        }

        result.vaporMoleFraction = vaporMoleFraction;
        result.mixture = mixture.state;
        result.enthalpy = result.mixture->h;
        result.entropy = result.mixture->s;
        result.residual = result.enthalpy - targetEnthalpy;
        result.caloricStatus = PRPhaseCaloricStatus::SUCCESS;
        result.phaseStatus = PhaseEquilibriumStatus::SUCCESS;
        result.status = PHFlashStatus::SUCCESS;
        result.valid = true;
        return result;
    }

    PhaseRoot root = PhaseRoot::VAPOR;
    if (result.tp.region == EquilibriumRegion::SINGLE_LIQUID)
        root = PhaseRoot::LIQUID;
    else if (result.tp.region != EquilibriumRegion::SINGLE_VAPOR
        && result.tp.region != EquilibriumRegion::SINGLE_SUPERCRITICAL)
    {
        result.status = PHFlashStatus::TP_CLOSURE_FAILURE;
        return result;
    }

    const PRPhaseCaloricResult single = caloric_.evaluate(
        pressure,
        temperature,
        result.tp.singlePhaseMoleFractions,
        root);
    result.caloricStatus = single.status;
    result.phaseStatus = single.phaseStatus;
    if (!single.success())
    {
        result.status =
            single.status == PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE
            ? PHFlashStatus::UNSUPPORTED_CALORIC_RANGE
            : PHFlashStatus::CALORIC_FAILURE;
        return result;
    }

    if (relativeResidual(
            single.state->compressibility,
            result.tp.singlePhaseState.compressibility)
            > options_.rootRelativeTolerance)
    {
        result.status = PHFlashStatus::TP_CLOSURE_FAILURE;
        return result;
    }

    result.mixtureMaterialResidual = 0.0;
    if (single.state->massFractions.size() != massFractions.size())
    {
        result.status = PHFlashStatus::CONSERVED_MIXTURE_FAILURE;
        return result;
    }
    for (std::size_t species = 0; species < massFractions.size(); ++species)
    {
        result.mixtureMaterialResidual = std::max(
            result.mixtureMaterialResidual,
            std::abs(single.state->massFractions[species]
                - massFractions[species]));
    }
    if (result.mixtureMaterialResidual
        > options_.conservedTolerances.compositionAbsolute)
    {
        result.status = PHFlashStatus::CONSERVED_MIXTURE_FAILURE;
        return result;
    }

    result.single = single.state;
    result.conservedStatus = ConservedMixtureStatus::SUCCESS;
    result.phaseEnergyIdentityResidual =
        single.energyIdentityRelativeResidual;
    result.mixtureEnergyIdentityResidual =
        single.energyIdentityRelativeResidual;
    result.vaporMoleFraction = result.tp.vaporMoleFraction;
    result.enthalpy = result.single->h;
    result.entropy = result.single->s;
    result.residual = result.enthalpy - targetEnthalpy;
    result.caloricStatus = PRPhaseCaloricStatus::SUCCESS;
    result.phaseStatus = PhaseEquilibriumStatus::SUCCESS;
    result.status = PHFlashStatus::SUCCESS;
    result.valid = true;
    return result;
}

PRPHFlash::TemperatureEvaluation PRPHFlash::evaluatePureCoexistence(
    double pressure,
    double targetEnthalpy,
    const std::vector<double>& massFractions,
    const PureSaturationResult& saturation) const
{
    TemperatureEvaluation result;
    if (!saturation.valid())
    {
        result.status = PHFlashStatus::TP_FLASH_FAILURE;
        result.phaseStatus = saturation.status;
        return result;
    }
    result = evaluateTemperature(
        pressure,
        targetEnthalpy,
        massFractions,
        saturation.temperature);
    if (result.tp.valid()
        && result.tp.region
            != EquilibriumRegion::PURE_COMPONENT_COEXISTENCE)
    {
        result.valid = false;
        result.status = PHFlashStatus::TP_CLOSURE_FAILURE;
    }
    return result;
}

PHFlashResult PRPHFlash::commit(
    double pressure,
    double targetEnthalpy,
    const TemperatureEvaluation& evaluation,
    std::size_t iterations,
    std::size_t validSamples,
    std::size_t invalidSamples) const
{
    PHFlashResult result;
    result.status = evaluation.status;
    result.phaseStatus = evaluation.phaseStatus;
    result.caloricStatus = evaluation.caloricStatus;
    result.conservedStatus = evaluation.conservedStatus;
    result.pressure = pressure;
    result.targetEnthalpy = targetEnthalpy;
    result.temperature = evaluation.temperature;
    result.vaporMoleFraction = evaluation.vaporMoleFraction;
    result.computedEnthalpy = evaluation.enthalpy;
    result.enthalpyResidual = evaluation.residual;
    result.enthalpyRelativeResidual = relativeResidual(
        evaluation.enthalpy, targetEnthalpy);
    result.entropy = evaluation.entropy;
    result.tpFugacityResidual = evaluation.tpFugacityResidual;
    result.tpMaterialResidual = evaluation.tpMaterialResidual;
    result.phaseEnergyIdentityResidual =
        evaluation.phaseEnergyIdentityResidual;
    result.mixtureMaterialResidual = evaluation.mixtureMaterialResidual;
    result.mixtureEnergyIdentityResidual =
        evaluation.mixtureEnergyIdentityResidual;
    result.pureComponentEnthalpyLeverRule = evaluation.pureLeverRule;
    result.iterations = iterations;
    result.validTemperatureSamples = validSamples;
    result.invalidTemperatureSamples = invalidSamples;
    result.tpFlash = evaluation.tp;

    if (!evaluation.valid) return result;
    if (!enthalpyConverged(
            evaluation.enthalpy, targetEnthalpy, options_))
    {
        result.status = PHFlashStatus::ENTHALPY_RESIDUAL_FAILURE;
        return result;
    }

    result.status = PHFlashStatus::SUCCESS;
    result.liquidPhase = evaluation.liquid;
    result.vaporPhase = evaluation.vapor;
    result.singlePhase = evaluation.single;
    result.conservedMixture = evaluation.mixture;
    return result;
}

PHFlashResult PRPHFlash::flashMassFractions(
    double pressure,
    double targetEnthalpy,
    const std::vector<double>& overallMassFractions) const
{
    PHFlashResult result;
    result.status = configurationStatus_;
    result.pressure = pressure;
    result.targetEnthalpy = targetEnthalpy;
    result.phaseStatus = tpFlash_.configurationStatus();
    result.caloricStatus = caloric_.configurationStatus();
    if (configurationStatus_ != PHFlashStatus::SUCCESS) return result;
    if (!finitePositive(pressure))
    {
        result.status = PHFlashStatus::INVALID_PRESSURE;
        result.phaseStatus = PhaseEquilibriumStatus::INVALID_PRESSURE;
        return result;
    }
    if (!std::isfinite(targetEnthalpy))
    {
        result.status = PHFlashStatus::NONFINITE_TARGET_ENTHALPY;
        return result;
    }

    std::vector<double> overallMoleFractions;
    result.phaseStatus = tpFlash_.kernel().massToMoleFractions(
        overallMassFractions, overallMoleFractions);
    if (result.phaseStatus != PhaseEquilibriumStatus::SUCCESS)
    {
        result.status = PHFlashStatus::INVALID_COMPOSITION;
        return result;
    }

    const PRPhaseCaloricTemperatureRange caloricRange =
        caloric_.supportedTemperatureRange(overallMassFractions);
    result.phaseStatus = caloricRange.phaseStatus;
    result.caloricStatus = caloricRange.status;
    if (!caloricRange.success())
    {
        result.status = caloricRange.status
                == PRPhaseCaloricStatus::INVALID_COMPOSITION
            ? PHFlashStatus::INVALID_COMPOSITION
            : PHFlashStatus::NO_VALID_CALORIC_INTERVAL;
        return result;
    }

    const double scanMinimum = std::max(
        options_.minimumTemperature, caloricRange.minimum);
    const double scanMaximum = std::min(
        options_.maximumTemperature, caloricRange.maximum);
    if (!finitePositive(scanMinimum) || !finitePositive(scanMaximum)
        || scanMinimum > scanMaximum)
    {
        result.status = PHFlashStatus::NO_VALID_CALORIC_INTERVAL;
        result.caloricStatus =
            PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE;
        return result;
    }

    std::size_t validSamples = 0;
    std::size_t invalidSamples = 0;
    TemperatureEvaluation lastEvaluation;
    bool haveLastEvaluation = false;

    const auto oppositeResiduals = []
    (
        const TemperatureEvaluation& lhs,
        const TemperatureEvaluation& rhs
    )
    {
        return (lhs.residual < 0.0 && rhs.residual > 0.0)
            || (lhs.residual > 0.0 && rhs.residual < 0.0);
    };

    const auto locateOppositeBeforeInvalid = [&]
    (
        const TemperatureEvaluation& validEndpoint,
        const TemperatureEvaluation& invalidEndpoint
    ) -> std::optional<TemperatureEvaluation>
    {
        if (!validEndpoint.valid || invalidEndpoint.valid) return std::nullopt;

        TemperatureEvaluation valid = validEndpoint;
        TemperatureEvaluation invalid = invalidEndpoint;
        for (std::size_t iteration = 0;
             iteration < options_.maximumIterations;
             ++iteration)
        {
            const double middleTemperature =
                0.5*(valid.temperature + invalid.temperature);
            if (!(middleTemperature > std::min(
                        valid.temperature, invalid.temperature)
                    && middleTemperature < std::max(
                        valid.temperature, invalid.temperature)))
            {
                break;
            }

            TemperatureEvaluation middle = evaluateTemperature(
                pressure,
                targetEnthalpy,
                overallMassFractions,
                middleTemperature);
            haveLastEvaluation = true;
            lastEvaluation = middle;
            if (middle.valid)
            {
                ++validSamples;
                if (enthalpyConverged(
                        middle.enthalpy, targetEnthalpy, options_)
                    || oppositeResiduals(validEndpoint, middle))
                {
                    return middle;
                }
                valid = std::move(middle);
            }
            else
            {
                ++invalidSamples;
                invalid = std::move(middle);
            }

            if (std::abs(invalid.temperature - valid.temperature)
                <= temperatureTolerance(
                    valid.temperature, invalid.temperature, options_))
            {
                break;
            }
        }
        return std::nullopt;
    };

    // TP alone cannot determine beta for a pure coexistence state. Resolve
    // that point first with the PH enthalpy lever before scanning h(T).
    if (speciesNames_.size() == 1)
    {
        const PureSaturationResult saturation =
            tpFlash_.solvePureSaturationTemperature(
                pressure,
                scanMinimum,
                scanMaximum);
        if (saturation.valid())
        {
            TemperatureEvaluation coexistence = evaluatePureCoexistence(
                pressure,
                targetEnthalpy,
                overallMassFractions,
                saturation);
            haveLastEvaluation = true;
            lastEvaluation = coexistence;
            if (coexistence.valid)
            {
                ++validSamples;
                if (enthalpyConverged(
                        coexistence.enthalpy,
                        targetEnthalpy,
                        options_))
                    return commit(
                        pressure,
                        targetEnthalpy,
                        coexistence,
                        0,
                        validSamples,
                        invalidSamples);
            }
            else
            {
                ++invalidSamples;
                if (coexistence.status == PHFlashStatus::TP_CLOSURE_FAILURE
                    || coexistence.status
                        == PHFlashStatus::CALORIC_FAILURE
                    || coexistence.status
                        == PHFlashStatus::NONPOSITIVE_LATENT_ENTHALPY
                    || coexistence.status
                        == PHFlashStatus::CONSERVED_MIXTURE_FAILURE)
                {
                    return commit(
                        pressure,
                        targetEnthalpy,
                        coexistence,
                        0,
                        validSamples,
                        invalidSamples);
                }
            }
        }
        else if (saturation.status
                != PhaseEquilibriumStatus::PURE_SATURATION_NO_BRACKET)
        {
            TemperatureEvaluation failure;
            failure.status = PHFlashStatus::TP_FLASH_FAILURE;
            failure.phaseStatus = saturation.status;
            return commit(
                pressure,
                targetEnthalpy,
                failure,
                0,
                validSamples,
                invalidSamples);
        }
    }

    std::optional<TemperatureEvaluation> previous;
    const double temperatureSpan = scanMaximum - scanMinimum;

    if (temperatureSpan == 0.0)
    {
        TemperatureEvaluation only = evaluateTemperature(
            pressure,
            targetEnthalpy,
            overallMassFractions,
            scanMinimum);
        if (only.valid)
        {
            ++validSamples;
            if (enthalpyConverged(
                    only.enthalpy, targetEnthalpy, options_))
                return commit(
                    pressure,
                    targetEnthalpy,
                    only,
                    0,
                    validSamples,
                    invalidSamples);
            only.status = PHFlashStatus::NO_ENTHALPY_BRACKET;
        }
        else
        {
            ++invalidSamples;
            only.status = PHFlashStatus::NO_VALID_CALORIC_INTERVAL;
        }
        only.valid = false;
        return commit(
            pressure,
            targetEnthalpy,
            only,
            0,
            validSamples,
            invalidSamples);
    }

    for (std::size_t sample = 0;
         sample <= options_.scanIntervals;
         ++sample)
    {
        const double fraction = static_cast<double>(sample)
            / static_cast<double>(options_.scanIntervals);
        const double temperature = sample == options_.scanIntervals
            ? scanMaximum
            : scanMinimum + fraction*temperatureSpan;
        TemperatureEvaluation current = evaluateTemperature(
            pressure,
            targetEnthalpy,
            overallMassFractions,
            temperature);
        haveLastEvaluation = true;
        lastEvaluation = current;
        bool currentAlreadyCounted = false;

        if (!current.valid)
        {
            ++invalidSamples;
            if (previous.has_value())
            {
                const std::optional<TemperatureEvaluation> boundary =
                    locateOppositeBeforeInvalid(*previous, current);
                if (boundary.has_value())
                {
                    current = *boundary;
                    currentAlreadyCounted = true;
                    lastEvaluation = current;
                }
            }
            if (!current.valid)
            {
                previous.reset();
                continue;
            }
        }
        if (!currentAlreadyCounted) ++validSamples;

        if (enthalpyConverged(
                current.enthalpy, targetEnthalpy, options_))
            return commit(
                pressure,
                targetEnthalpy,
                current,
                0,
                validSamples,
                invalidSamples);

        if (previous.has_value()
            && oppositeResiduals(*previous, current))
        {
            TemperatureEvaluation lower = *previous;
            TemperatureEvaluation upper = current;
            TemperatureEvaluation best =
                std::abs(lower.residual) <= std::abs(upper.residual)
                ? lower : upper;
            const auto retainBest = [&best]
            (
                const TemperatureEvaluation& candidate
            )
            {
                if (candidate.valid
                    && std::abs(candidate.residual) < std::abs(best.residual))
                {
                    best = candidate;
                }
            };
            bool bracketBroken = false;

            for (std::size_t iteration = 1;
                 iteration <= options_.maximumIterations;
                 ++iteration)
            {
                const double middleTemperature =
                    0.5*(lower.temperature + upper.temperature);
                if (!(middleTemperature > lower.temperature
                    && middleTemperature < upper.temperature))
                {
                    lastEvaluation = best;
                    lastEvaluation.status = PHFlashStatus::NOT_CONVERGED;
                    return commit(
                        pressure,
                        targetEnthalpy,
                        lastEvaluation,
                        iteration,
                        validSamples,
                        invalidSamples);
                }

                TemperatureEvaluation middle = evaluateTemperature(
                    pressure,
                    targetEnthalpy,
                    overallMassFractions,
                    middleTemperature);
                haveLastEvaluation = true;
                lastEvaluation = middle;
                if (!middle.valid)
                {
                    ++invalidSamples;
                    const std::optional<TemperatureEvaluation> leftBoundary =
                        locateOppositeBeforeInvalid(lower, middle);
                    if (leftBoundary.has_value())
                    {
                        if (enthalpyConverged(
                                leftBoundary->enthalpy,
                                targetEnthalpy,
                                options_))
                        {
                            return commit(
                                pressure,
                                targetEnthalpy,
                                *leftBoundary,
                                iteration,
                                validSamples,
                                invalidSamples);
                        }
                        upper = *leftBoundary;
                        retainBest(*leftBoundary);
                        continue;
                    }

                    const std::optional<TemperatureEvaluation> rightBoundary =
                        locateOppositeBeforeInvalid(upper, middle);
                    if (rightBoundary.has_value())
                    {
                        if (enthalpyConverged(
                                rightBoundary->enthalpy,
                                targetEnthalpy,
                                options_))
                        {
                            return commit(
                                pressure,
                                targetEnthalpy,
                                *rightBoundary,
                                iteration,
                                validSamples,
                                invalidSamples);
                        }
                        lower = *rightBoundary;
                        retainBest(*rightBoundary);
                        continue;
                    }

                    // Neither continuous side of the invalid state contains
                    // the sign change. Do not bridge the topology hole.
                    bracketBroken = true;
                    break;
                }
                ++validSamples;
                retainBest(middle);

                if (enthalpyConverged(
                        middle.enthalpy, targetEnthalpy, options_))
                    return commit(
                        pressure,
                        targetEnthalpy,
                        middle,
                        iteration,
                        validSamples,
                        invalidSamples);

                if (oppositeResiduals(lower, middle))
                    upper = std::move(middle);
                else
                    lower = std::move(middle);

                if (iteration == options_.maximumIterations)
                {
                    lastEvaluation = best;
                    lastEvaluation.status = PHFlashStatus::NOT_CONVERGED;
                    return commit(
                        pressure,
                        targetEnthalpy,
                        lastEvaluation,
                        iteration,
                        validSamples,
                        invalidSamples);
                }
            }

            if (bracketBroken)
            {
                // The failed bracket says nothing about the interval to the
                // right of the current scan point. Retain that point as the
                // first endpoint of a new, independently checked segment.
                previous = std::move(current);
                continue;
            }
        }

        previous = std::move(current);
    }

    if (!haveLastEvaluation)
    {
        lastEvaluation.status = PHFlashStatus::NO_VALID_CALORIC_INTERVAL;
    }
    else if (validSamples == 0)
    {
        lastEvaluation.status =
            PHFlashStatus::NO_VALID_CALORIC_INTERVAL;
    }
    else
    {
        lastEvaluation.status = PHFlashStatus::NO_ENTHALPY_BRACKET;
    }
    lastEvaluation.valid = false;
    return commit(
        pressure,
        targetEnthalpy,
        lastEvaluation,
        0,
        validSamples,
        invalidSamples);
}

} // End namespace phaseEquilibrium
} // End namespace Foam
