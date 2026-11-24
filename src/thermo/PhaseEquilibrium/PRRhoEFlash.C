#include "PRRhoEFlash.H"

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

bool densityConverged(
    double computed,
    double target,
    const RhoEFlashOptions& options)
{
    if (!finitePositive(computed) || !finitePositive(target)) return false;
    return std::abs(computed - target)
        <= options.densityAbsoluteTolerance
            + options.densityRelativeTolerance
                *std::max({std::abs(computed), std::abs(target), 1.0});
}

bool energyConverged(
    double computed,
    double target,
    const RhoEFlashOptions& options)
{
    if (!std::isfinite(computed) || !std::isfinite(target)) return false;
    return std::abs(computed - target)
        <= options.energyAbsoluteTolerance
            + options.energyRelativeTolerance
                *std::max({std::abs(computed), std::abs(target), 1.0});
}

double scaledDensityResidual(double computed, double target)
{
    return (computed - target)/std::max({std::abs(target), 1.0});
}

double scaledEnergyResidual(double computed, double target)
{
    return (computed - target)/std::max({std::abs(target), 1.0});
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
    if (tp.region == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE)
        return false;

    if (tp.region == EquilibriumRegion::TWO_PHASE)
    {
        if (!tp.phaseFractionDetermined
            || !std::isfinite(tp.vaporMoleFraction)
            || tp.vaporMoleFraction < 0.0
            || tp.vaporMoleFraction > 1.0
            || tp.liquidMoleFractions.size() != feed.size()
            || tp.vaporMoleFractions.size() != feed.size()
            || !tp.liquidState.valid() || !tp.vaporState.valid())
            return false;

        for (std::size_t species = 0; species < feed.size(); ++species)
        {
            const double liquid = tp.liquidMoleFractions[species];
            const double vapor = tp.vaporMoleFractions[species];
            if (!finiteNonNegative(liquid) || !finiteNonNegative(vapor))
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
        return std::isfinite(fugacityResidual)
            && std::isfinite(materialResidual)
            && std::isfinite(tp.rachfordRiceResidual)
            && fugacityResidual <= options.fugacityTolerance
            && materialResidual <= options.materialTolerance
            && std::abs(tp.rachfordRiceResidual)
                <= options.materialTolerance;
    }

    if (!tp.singlePhaseState.valid()
        || tp.singlePhaseMoleFractions.size() != feed.size())
        return false;
    for (std::size_t species = 0; species < feed.size(); ++species)
    {
        const double value = tp.singlePhaseMoleFractions[species];
        if (!finiteNonNegative(value)) return false;
        materialResidual = std::max(
            materialResidual, std::abs(value - feed[species]));
    }
    return std::isfinite(materialResidual)
        && materialResidual <= options.materialTolerance;
}

PhaseRoot rootForRegion(EquilibriumRegion region)
{
    return region == EquilibriumRegion::SINGLE_LIQUID
        ? PhaseRoot::LIQUID : PhaseRoot::VAPOR;
}

} // End anonymous namespace

const char* rhoEFlashStatusName(RhoEFlashStatus status)
{
    switch (status)
    {
        case RhoEFlashStatus::SUCCESS: return "SUCCESS";
        case RhoEFlashStatus::INVALID_CONFIGURATION:
            return "INVALID_CONFIGURATION";
        case RhoEFlashStatus::INVALID_OPTIONS: return "INVALID_OPTIONS";
        case RhoEFlashStatus::INVALID_DENSITY: return "INVALID_DENSITY";
        case RhoEFlashStatus::NONFINITE_TARGET_ENERGY:
            return "NONFINITE_TARGET_ENERGY";
        case RhoEFlashStatus::INVALID_COMPOSITION:
            return "INVALID_COMPOSITION";
        case RhoEFlashStatus::OUT_OF_CONFIGURED_RANGE:
            return "OUT_OF_CONFIGURED_RANGE";
        case RhoEFlashStatus::REQUIRED_REGION_MISMATCH:
            return "REQUIRED_REGION_MISMATCH";
        case RhoEFlashStatus::TP_FLASH_FAILURE: return "TP_FLASH_FAILURE";
        case RhoEFlashStatus::TP_CLOSURE_FAILURE:
            return "TP_CLOSURE_FAILURE";
        case RhoEFlashStatus::CALORIC_FAILURE: return "CALORIC_FAILURE";
        case RhoEFlashStatus::UNSUPPORTED_CALORIC_RANGE:
            return "UNSUPPORTED_CALORIC_RANGE";
        case RhoEFlashStatus::CONSERVED_MIXTURE_FAILURE:
            return "CONSERVED_MIXTURE_FAILURE";
        case RhoEFlashStatus::PURE_SATURATION_NO_BRACKET:
            return "PURE_SATURATION_NO_BRACKET";
        case RhoEFlashStatus::NO_VALID_SATURATION_INTERVAL:
            return "NO_VALID_SATURATION_INTERVAL";
        case RhoEFlashStatus::NO_RHOE_BRACKET:
            return "NO_RHOE_BRACKET";
        case RhoEFlashStatus::NOT_CONVERGED: return "NOT_CONVERGED";
        case RhoEFlashStatus::RESIDUAL_FAILURE:
            return "RESIDUAL_FAILURE";
    }
    return "UNKNOWN_STATUS";
}

PRRhoEFlash::PRRhoEFlash(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    const RhoEFlashOptions& options)
:
    speciesNames_(speciesNames),
    tpFlash_(criticalData, kij, speciesNames, options.tpOptions),
    caloric_(criticalData, kij, speciesNames),
    options_(options),
    configurationStatus_(RhoEFlashStatus::SUCCESS)
{
    if (!finitePositive(options_.minimumPressure)
        || !finitePositive(options_.maximumPressure)
        || options_.minimumPressure >= options_.maximumPressure
        || !finitePositive(options_.minimumTemperature)
        || !finitePositive(options_.maximumTemperature)
        || options_.minimumTemperature >= options_.maximumTemperature
        || options_.pressureScanIntervals == 0
        || options_.maximumIterations == 0
        || !finitePositive(options_.finiteDifferenceStep)
        || options_.finiteDifferenceStep >= 0.5
        || !finitePositive(options_.saturationLogPressureTolerance)
        || !finitePositive(options_.logStepLimit)
        || !finiteNonNegative(options_.densityAbsoluteTolerance)
        || !finiteNonNegative(options_.densityRelativeTolerance)
        || !finiteNonNegative(options_.energyAbsoluteTolerance)
        || !finiteNonNegative(options_.energyRelativeTolerance)
        || !finiteNonNegative(options_.tpRootRelativeTolerance)
        || !finiteNonNegative(
            options_.conservedTolerances.compositionAbsolute)
        || !finiteNonNegative(options_.conservedTolerances.energyRelative))
    {
        configurationStatus_ = RhoEFlashStatus::INVALID_OPTIONS;
        return;
    }
    if (tpFlash_.configurationStatus() != PhaseEquilibriumStatus::SUCCESS
        || !caloric_.configured())
    {
        configurationStatus_ = RhoEFlashStatus::INVALID_CONFIGURATION;
    }
}

PRRhoEFlash::Evaluation PRRhoEFlash::assemblePT(
    double pressure,
    double temperature,
    const std::vector<double>& massFractions) const
{
    Evaluation result;
    result.pressure = pressure;
    result.temperature = temperature;
    result.tp = tpFlash_.flashMassFractions(
        pressure, temperature, massFractions);
    result.phaseStatus = result.tp.status;
    if (!result.tp.valid())
    {
        result.status = RhoEFlashStatus::TP_FLASH_FAILURE;
        return result;
    }
    if (!verifyTPClosure(
            result.tp,
            options_.tpOptions,
            result.tpFugacityResidual,
            result.tpMaterialResidual))
    {
        result.status = RhoEFlashStatus::TP_CLOSURE_FAILURE;
        return result;
    }

    if (result.tp.region == EquilibriumRegion::TWO_PHASE)
    {
        const PRPhaseCaloricResult liquid = caloric_.evaluate(
            pressure, temperature, result.tp.liquidMoleFractions,
            PhaseRoot::LIQUID);
        const PRPhaseCaloricResult vapor = caloric_.evaluate(
            pressure, temperature, result.tp.vaporMoleFractions,
            PhaseRoot::VAPOR);
        result.caloricStatus = !liquid.success()
            ? liquid.status : vapor.status;
        result.phaseStatus = !liquid.success()
            ? liquid.phaseStatus : vapor.phaseStatus;
        if (!liquid.success() || !vapor.success())
        {
            result.status = result.caloricStatus ==
                PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE
                ? RhoEFlashStatus::UNSUPPORTED_CALORIC_RANGE
                : RhoEFlashStatus::CALORIC_FAILURE;
            return result;
        }
        result.liquid = liquid.state;
        result.vapor = vapor.state;
        result.phaseEnergyIdentityResidual = std::max(
            liquid.energyIdentityRelativeResidual,
            vapor.energyIdentityRelativeResidual);
        const ConservedMixtureResult mixture =
            assembleConservedPhaseMixture(
                pressure,
                result.tp.vaporMoleFraction,
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
            result.status = RhoEFlashStatus::CONSERVED_MIXTURE_FAILURE;
            return result;
        }
        result.mixture = mixture.state;
        result.vaporMoleFraction = result.tp.vaporMoleFraction;
        result.density = result.mixture->rho;
        result.energy = result.mixture->e;
        result.entropy = result.mixture->s;
    }
    else
    {
        if (result.tp.region == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE)
        {
            result.status = RhoEFlashStatus::TP_CLOSURE_FAILURE;
            return result;
        }
        const PhaseRoot root = rootForRegion(result.tp.region);
        const PRPhaseCaloricResult single = caloric_.evaluate(
            pressure, temperature, result.tp.singlePhaseMoleFractions, root);
        result.caloricStatus = single.status;
        result.phaseStatus = single.phaseStatus;
        if (!single.success())
        {
            result.status = single.status ==
                PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE
                ? RhoEFlashStatus::UNSUPPORTED_CALORIC_RANGE
                : RhoEFlashStatus::CALORIC_FAILURE;
            return result;
        }
        if (relativeResidual(
                single.state->compressibility,
                result.tp.singlePhaseState.compressibility)
                > options_.tpRootRelativeTolerance)
        {
            result.status = RhoEFlashStatus::TP_CLOSURE_FAILURE;
            return result;
        }
        result.single = single.state;
        result.mixtureMaterialResidual = 0.0;
        for (std::size_t species = 0;
             species < massFractions.size(); ++species)
        {
            result.mixtureMaterialResidual = std::max(
                result.mixtureMaterialResidual,
                std::abs(single.state->massFractions[species]
                    - massFractions[species]));
        }
        if (result.mixtureMaterialResidual
            > options_.conservedTolerances.compositionAbsolute)
        {
            result.status = RhoEFlashStatus::CONSERVED_MIXTURE_FAILURE;
            return result;
        }
        result.conservedStatus = ConservedMixtureStatus::SUCCESS;
        result.phaseEnergyIdentityResidual =
            single.energyIdentityRelativeResidual;
        result.mixtureEnergyIdentityResidual =
            single.energyIdentityRelativeResidual;
        result.density = single.state->rho;
        result.energy = single.state->e;
        result.entropy = single.state->s;
    }

    if (!finitePositive(result.density) || !std::isfinite(result.energy))
    {
        result.status = RhoEFlashStatus::RESIDUAL_FAILURE;
        return result;
    }
    result.status = RhoEFlashStatus::SUCCESS;
    result.valid = true;
    return result;
}

PRRhoEFlash::Evaluation PRRhoEFlash::evaluatePT(
    double pressure,
    double temperature,
    double targetDensity,
    double targetEnergy,
    const std::vector<double>& massFractions) const
{
    Evaluation result = assemblePT(pressure, temperature, massFractions);
    if (!result.valid) return result;

    result.densityResidual = result.density - targetDensity;
    result.densityRelativeResidual =
        relativeResidual(result.density, targetDensity);
    result.energyResidual = result.energy - targetEnergy;
    result.energyRelativeResidual = relativeResidual(
        result.energy, targetEnergy);
    if (!std::isfinite(result.densityResidual)
        || !std::isfinite(result.energyResidual))
    {
        result.status = RhoEFlashStatus::RESIDUAL_FAILURE;
        result.valid = false;
    }
    return result;
}

PRRhoEFlash::Evaluation PRRhoEFlash::evaluatePureSaturationPressure(
    double pressure,
    double targetDensity,
    double targetEnergy,
    const std::vector<double>& massFractions) const
{
    Evaluation result;
    result.pressure = pressure;
    const PureSaturationResult saturation =
        tpFlash_.solvePureSaturationTemperature(
            pressure,
            options_.minimumTemperature,
            options_.maximumTemperature);
    if (!saturation.valid())
    {
        result.status = saturation.status ==
            PhaseEquilibriumStatus::PURE_SATURATION_NO_BRACKET
            ? RhoEFlashStatus::PURE_SATURATION_NO_BRACKET
            : RhoEFlashStatus::TP_FLASH_FAILURE;
        result.phaseStatus = saturation.status;
        return result;
    }
    result = evaluatePT(
        pressure,
        saturation.temperature,
        targetDensity,
        targetEnergy,
        massFractions);
    if (result.status == RhoEFlashStatus::TP_CLOSURE_FAILURE)
    {
        // TP at a pure saturation point intentionally has an underdetermined
        // beta. Rebuild the phase payload explicitly and determine beta from
        // the target specific volume below.
        result = Evaluation();
        result.pressure = pressure;
        result.temperature = saturation.temperature;
        result.tp = tpFlash_.flashMassFractions(
            pressure, saturation.temperature, massFractions);
        result.phaseStatus = result.tp.status;
        if (!result.tp.valid()
            || result.tp.region != EquilibriumRegion::PURE_COMPONENT_COEXISTENCE)
        {
            result.status = RhoEFlashStatus::TP_CLOSURE_FAILURE;
            return result;
        }
        const PRPhaseCaloricResult liquid = caloric_.evaluate(
            pressure, saturation.temperature, massFractions, PhaseRoot::LIQUID);
        const PRPhaseCaloricResult vapor = caloric_.evaluate(
            pressure, saturation.temperature, massFractions, PhaseRoot::VAPOR);
        result.caloricStatus = !liquid.success()
            ? liquid.status : vapor.status;
        if (!liquid.success() || !vapor.success())
        {
            result.status = result.caloricStatus ==
                PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE
                ? RhoEFlashStatus::UNSUPPORTED_CALORIC_RANGE
                : RhoEFlashStatus::CALORIC_FAILURE;
            return result;
        }
        result.liquid = liquid.state;
        result.vapor = vapor.state;
        result.phaseEnergyIdentityResidual = std::max(
            liquid.energyIdentityRelativeResidual,
            vapor.energyIdentityRelativeResidual);
        const double liquidVolume = 1.0/result.liquid->rho;
        const double vaporVolume = 1.0/result.vapor->rho;
        const double targetVolume = 1.0/targetDensity;
        const double volumeDenominator = vaporVolume - liquidVolume;
        if (!finitePositive(volumeDenominator))
        {
            result.status = RhoEFlashStatus::NO_RHOE_BRACKET;
            return result;
        }
        double betaMass = (targetVolume - liquidVolume)/volumeDenominator;
        const double fractionTolerance = 1.0e-10;
        if (!std::isfinite(betaMass)
            || betaMass < -fractionTolerance
            || betaMass > 1.0 + fractionTolerance)
        {
            result.status = RhoEFlashStatus::NO_RHOE_BRACKET;
            return result;
        }
        if (betaMass < 0.0) betaMass = 0.0;
        if (betaMass > 1.0) betaMass = 1.0;
        const ConservedMixtureResult mixture =
            assembleConservedPhaseMixture(
                pressure,
                betaMass,
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
            result.status = RhoEFlashStatus::CONSERVED_MIXTURE_FAILURE;
            return result;
        }
        result.mixture = mixture.state;
        result.vaporMoleFraction = betaMass;
        result.density = result.mixture->rho;
        result.energy = result.mixture->e;
        result.entropy = result.mixture->s;
        result.densityResidual = result.density - targetDensity;
        result.densityRelativeResidual =
            relativeResidual(result.density, targetDensity);
        result.energyResidual = result.energy - targetEnergy;
        result.energyRelativeResidual =
            relativeResidual(result.energy, targetEnergy);
        result.pureVolumeEnergyLeverRule = true;
        result.tpFugacityResidual = result.tp.maximumLogFugacityResidual;
        result.tpMaterialResidual = result.tp.maximumMoleMaterialResidual;
        result.status = RhoEFlashStatus::SUCCESS;
        result.valid = true;
    }
    return result;
}

RhoEFlashResult PRRhoEFlash::commit(
    double targetDensity,
    double targetEnergy,
    const Evaluation& evaluation,
    std::size_t iterations,
    std::size_t validSamples,
    std::size_t invalidSamples) const
{
    RhoEFlashResult result;
    result.status = evaluation.status;
    result.phaseStatus = evaluation.phaseStatus;
    result.caloricStatus = evaluation.caloricStatus;
    result.conservedStatus = evaluation.conservedStatus;
    result.targetDensity = targetDensity;
    result.targetEnergy = targetEnergy;
    result.pressure = evaluation.pressure;
    result.temperature = evaluation.temperature;
    result.density = evaluation.density;
    result.energy = evaluation.energy;
    result.entropy = evaluation.entropy;
    result.densityResidual = evaluation.densityResidual;
    result.densityRelativeResidual = evaluation.densityRelativeResidual;
    result.energyResidual = evaluation.energyResidual;
    result.energyRelativeResidual = evaluation.energyRelativeResidual;
    result.tpFugacityResidual = evaluation.tpFugacityResidual;
    result.tpMaterialResidual = evaluation.tpMaterialResidual;
    result.phaseEnergyIdentityResidual =
        evaluation.phaseEnergyIdentityResidual;
    result.mixtureMaterialResidual = evaluation.mixtureMaterialResidual;
    result.mixtureEnergyIdentityResidual =
        evaluation.mixtureEnergyIdentityResidual;
    result.vaporMoleFraction = evaluation.vaporMoleFraction;
    result.pureVolumeEnergyLeverRule = evaluation.pureVolumeEnergyLeverRule;
    result.iterations = iterations;
    result.validSamples = validSamples;
    result.invalidSamples = invalidSamples;
    result.tpFlash = evaluation.tp;

    if (!evaluation.valid) return result;
    if (!densityConverged(evaluation.density, targetDensity, options_)
        || !energyConverged(evaluation.energy, targetEnergy, options_))
    {
        result.status = evaluation.status == RhoEFlashStatus::NOT_CONVERGED
            ? RhoEFlashStatus::NOT_CONVERGED
            : RhoEFlashStatus::RESIDUAL_FAILURE;
        return result;
    }
    result.status = RhoEFlashStatus::SUCCESS;
    result.liquidPhase = evaluation.liquid;
    result.vaporPhase = evaluation.vapor;
    result.singlePhase = evaluation.single;
    result.conservedMixture = evaluation.mixture;
    return result;
}

std::optional<RhoEFlashResult> PRRhoEFlash::refinePureSaturation(
    double targetDensity,
    double targetEnergy,
    const std::vector<double>& massFractions,
    Evaluation& best,
    std::size_t& iterationsPerformed,
    std::size_t& validSamples,
    std::size_t& invalidSamples) const
{
    if (!best.valid || !best.pureVolumeEnergyLeverRule
        || best.tp.region != EquilibriumRegion::PURE_COMPONENT_COEXISTENCE)
        return std::nullopt;
    Evaluation current = best;
    double currentLogPressure = std::log(current.pressure);
    const double logPressureMinimum = std::log(options_.minimumPressure);
    const double logPressureMaximum = std::log(options_.maximumPressure);
    for (std::size_t iteration = 1;
         iteration <= options_.maximumIterations; ++iteration)
    {
        iterationsPerformed = iteration;
        if (densityConverged(current.density, targetDensity, options_)
            && energyConverged(current.energy, targetEnergy, options_))
        {
            return commit(
                targetDensity, targetEnergy, current, iteration,
                validSamples, invalidSamples);
        }

        double h = options_.finiteDifferenceStep;
        Evaluation plus;
        Evaluation minus;
        double plusLogPressure = currentLogPressure;
        double minusLogPressure = currentLogPressure;
        bool plusValid = false;
        bool minusValid = false;
        for (int refinement = 0; refinement < 20; ++refinement)
        {
            plusLogPressure = std::min(
                logPressureMaximum, currentLogPressure + h);
            minusLogPressure = std::max(
                logPressureMinimum, currentLogPressure - h);
            plus = evaluatePureSaturationPressure(
                std::exp(plusLogPressure), targetDensity,
                targetEnergy, massFractions);
            minus = evaluatePureSaturationPressure(
                std::exp(minusLogPressure), targetDensity,
                targetEnergy, massFractions);
            plusValid = plus.valid
                && plusLogPressure > currentLogPressure;
            minusValid = minus.valid
                && minusLogPressure < currentLogPressure;
            if (plusValid) ++validSamples;
            else ++invalidSamples;
            if (minusValid) ++validSamples;
            else ++invalidSamples;
            if (plusValid || minusValid) break;
            h *= 0.5;
        }
        if (!plusValid && !minusValid) break;

        double derivative = std::numeric_limits<double>::quiet_NaN();
        if (plusValid && minusValid)
        {
            derivative = (plus.energyResidual - minus.energyResidual)
                /(plusLogPressure - minusLogPressure);
        }
        else if (plusValid)
        {
            derivative = (plus.energyResidual - current.energyResidual)
                /(plusLogPressure - currentLogPressure);
        }
        else
        {
            derivative = (current.energyResidual - minus.energyResidual)
                /(currentLogPressure - minusLogPressure);
        }
        if (!std::isfinite(derivative)
            || std::abs(derivative) < 1.0e-12)
            break;

        double delta = -current.energyResidual/derivative;
        delta = std::clamp(
            delta, -options_.logStepLimit, options_.logStepLimit);
        bool accepted = false;
        for (int line = 0; line < 24; ++line)
        {
            const double trialLogPressure = std::clamp(
                currentLogPressure + delta,
                logPressureMinimum, logPressureMaximum);
            if (trialLogPressure == currentLogPressure) break;
            const double continuationCount = std::ceil(
                std::abs(trialLogPressure - currentLogPressure)
                /options_.finiteDifferenceStep);
            if (!std::isfinite(continuationCount)
                || continuationCount > 64.0)
            {
                delta *= 0.5;
                continue;
            }
            if (continuationCount < 1.0) break;
            const std::size_t continuationSamples =
                static_cast<std::size_t>(continuationCount);

            Evaluation trial;
            bool connected = true;
            for (std::size_t continuation = 1;
                 continuation <= continuationSamples; ++continuation)
            {
                const double fraction = static_cast<double>(continuation)
                    /static_cast<double>(continuationSamples);
                trial = evaluatePureSaturationPressure(
                    std::exp(currentLogPressure
                        + fraction*(trialLogPressure - currentLogPressure)),
                    targetDensity, targetEnergy, massFractions);
                if (!trial.valid)
                {
                    ++invalidSamples;
                    connected = false;
                    break;
                }
                ++validSamples;
            }
            if (connected
                && std::abs(trial.energyResidual)
                    < std::abs(current.energyResidual))
            {
                current = std::move(trial);
                currentLogPressure = trialLogPressure;
                accepted = true;
                if (std::abs(current.energyResidual)
                    < std::abs(best.energyResidual))
                    best = current;
                break;
            }
            delta *= 0.5;
        }
        if (!accepted) break;
    }
    return std::nullopt;
}

std::optional<RhoEFlashResult> PRRhoEFlash::solveGeneralSeed(
    double logPressure,
    double logTemperature,
    bool connectedOnly,
    EquilibriumRegion requiredRegion,
    double targetDensity,
    double targetEnergy,
    const std::vector<double>& massFractions,
    Evaluation& best,
    bool& haveBest,
    Evaluation& regionMismatch,
    bool& haveRegionMismatch,
    std::size_t& iterationsPerformed,
    std::size_t& validSamples,
    std::size_t& invalidSamples) const
{
    const double logPMin = std::log(options_.minimumPressure);
    const double logPMax = std::log(options_.maximumPressure);
    const double logTMin = std::log(options_.minimumTemperature);
    const double logTMax = std::log(options_.maximumTemperature);
    if (logPressure < logPMin || logPressure > logPMax
        || logTemperature < logTMin || logTemperature > logTMax)
        return std::nullopt;

    double logP = logPressure;
    double logT = logTemperature;
    auto onRequiredRegion = [&](const Evaluation& evaluation)
    {
        if (!evaluation.valid) return false;
        if (requiredRegion != EquilibriumRegion::UNKNOWN
            && evaluation.tp.region != requiredRegion)
        {
            regionMismatch = evaluation;
            haveRegionMismatch = true;
            return false;
        }
        return true;
    };
    Evaluation current = evaluatePT(
        std::exp(logP), std::exp(logT), targetDensity, targetEnergy,
        massFractions);
    if (!onRequiredRegion(current))
    {
        ++invalidSamples;
        return std::nullopt;
    }
    ++validSamples;
    if (!haveBest || std::hypot(
            scaledDensityResidual(current.density, targetDensity),
            scaledEnergyResidual(current.energy, targetEnergy))
            < std::hypot(
                scaledDensityResidual(best.density, targetDensity),
                scaledEnergyResidual(best.energy, targetEnergy)))
    {
        best = current;
        haveBest = true;
    }

    for (std::size_t iteration = 1;
         iteration <= options_.maximumIterations; ++iteration)
    {
        iterationsPerformed = iteration;
        if (densityConverged(current.density, targetDensity, options_)
            && energyConverged(current.energy, targetEnergy, options_))
        {
            return commit(
                targetDensity, targetEnergy, current, iteration,
                validSamples, invalidSamples);
        }
        const double h = options_.finiteDifferenceStep;
        Evaluation pPlus = evaluatePT(
            std::exp(std::min(logPMax, logP + h)),
            std::exp(logT), targetDensity, targetEnergy, massFractions);
        Evaluation pMinus = evaluatePT(
            std::exp(std::max(logPMin, logP - h)),
            std::exp(logT), targetDensity, targetEnergy, massFractions);
        Evaluation tPlus = evaluatePT(
            std::exp(logP),
            std::exp(std::min(logTMax, logT + h)),
            targetDensity, targetEnergy, massFractions);
        Evaluation tMinus = evaluatePT(
            std::exp(logP),
            std::exp(std::max(logTMin, logT - h)),
            targetDensity, targetEnergy, massFractions);
        const bool pPlusValid = onRequiredRegion(pPlus);
        const bool pMinusValid = onRequiredRegion(pMinus);
        const bool tPlusValid = onRequiredRegion(tPlus);
        const bool tMinusValid = onRequiredRegion(tMinus);
        validSamples += static_cast<std::size_t>(pPlusValid)
            + static_cast<std::size_t>(pMinusValid)
            + static_cast<std::size_t>(tPlusValid)
            + static_cast<std::size_t>(tMinusValid);
        invalidSamples += static_cast<std::size_t>(!pPlusValid)
            + static_cast<std::size_t>(!pMinusValid)
            + static_cast<std::size_t>(!tPlusValid)
            + static_cast<std::size_t>(!tMinusValid);
        if (!pPlusValid || !pMinusValid || !tPlusValid || !tMinusValid)
        {
            break;
        }
        const double pDenominator =
            std::min(logPMax, logP + h)
                - std::max(logPMin, logP - h);
        const double tDenominator =
            std::min(logTMax, logT + h)
                - std::max(logTMin, logT - h);
        if (!(pDenominator > 0.0 && tDenominator > 0.0)) break;
        const double rRho = scaledDensityResidual(
            current.density, targetDensity);
        const double rE = scaledEnergyResidual(
            current.energy, targetEnergy);
        const double j11 = (
            scaledDensityResidual(pPlus.density, targetDensity)
            - scaledDensityResidual(pMinus.density, targetDensity))
            /pDenominator;
        const double j21 = (
            scaledEnergyResidual(pPlus.energy, targetEnergy)
            - scaledEnergyResidual(pMinus.energy, targetEnergy))
            /pDenominator;
        const double j12 = (
            scaledDensityResidual(tPlus.density, targetDensity)
            - scaledDensityResidual(tMinus.density, targetDensity))
            /tDenominator;
        const double j22 = (
            scaledEnergyResidual(tPlus.energy, targetEnergy)
            - scaledEnergyResidual(tMinus.energy, targetEnergy))
            /tDenominator;
        const double determinant = j11*j22 - j12*j21;
        if (!std::isfinite(determinant)
            || std::abs(determinant) < 1.0e-14)
            break;
        double deltaP = (-rRho*j22 + j12*rE)/determinant;
        double deltaT = (-j11*rE + rRho*j21)/determinant;
        deltaP = std::clamp(
            deltaP, -options_.logStepLimit, options_.logStepLimit);
        deltaT = std::clamp(
            deltaT, -options_.logStepLimit, options_.logStepLimit);
        const double oldNorm = std::hypot(rRho, rE);
        bool accepted = false;
        for (int line = 0; line < 12; ++line)
        {
            const double trialP = std::clamp(
                logP + deltaP, logPMin, logPMax);
            const double trialT = std::clamp(
                logT + deltaT, logTMin, logTMax);
            Evaluation trial;
            bool connected = true;
            if (connectedOnly)
            {
                const double continuationCount = std::ceil(
                    std::max(
                        std::abs(trialP - logP),
                        std::abs(trialT - logT))
                    /options_.finiteDifferenceStep);
                if (!std::isfinite(continuationCount)
                    || continuationCount > 64.0)
                {
                    deltaP *= 0.5;
                    deltaT *= 0.5;
                    continue;
                }
                if (continuationCount < 1.0) break;
                const std::size_t continuationSamples =
                    static_cast<std::size_t>(continuationCount);
                for (std::size_t continuation = 1;
                     continuation <= continuationSamples; ++continuation)
                {
                    const double fraction =
                        static_cast<double>(continuation)
                        /static_cast<double>(continuationSamples);
                    trial = evaluatePT(
                        std::exp(logP + fraction*(trialP - logP)),
                        std::exp(logT + fraction*(trialT - logT)),
                        targetDensity, targetEnergy, massFractions);
                    if (!onRequiredRegion(trial))
                    {
                        ++invalidSamples;
                        connected = false;
                        break;
                    }
                    ++validSamples;
                }
            }
            else
            {
                trial = evaluatePT(
                    std::exp(trialP), std::exp(trialT), targetDensity,
                    targetEnergy, massFractions);
                connected = onRequiredRegion(trial);
                if (connected) ++validSamples;
                else ++invalidSamples;
            }
            if (connected)
            {
                const double trialNorm = std::hypot(
                    scaledDensityResidual(trial.density, targetDensity),
                    scaledEnergyResidual(trial.energy, targetEnergy));
                if (std::isfinite(trialNorm) && trialNorm < oldNorm)
                {
                    logP = trialP;
                    logT = trialT;
                    current = std::move(trial);
                    accepted = true;
                    if (!haveBest || trialNorm < std::hypot(
                            scaledDensityResidual(
                                best.density, targetDensity),
                            scaledEnergyResidual(
                                best.energy, targetEnergy)))
                    {
                        best = current;
                        haveBest = true;
                    }
                    break;
                }
            }
            deltaP *= 0.5;
            deltaT *= 0.5;
        }
        if (!accepted) break;
    }
    return std::nullopt;
}

RhoEFlashResult PRRhoEFlash::flashMassFractions(
    double targetDensity,
    double targetEnergy,
    const std::vector<double>& overallMassFractions) const
{
    return flashMassFractions(
        targetDensity, targetEnergy, overallMassFractions, RhoEFlashSeed{});
}

RhoEFlashResult PRRhoEFlash::initializePTYMassFractions(
    double pressure,
    double temperature,
    const std::vector<double>& overallMassFractions) const
{
    RhoEFlashResult invalid;
    invalid.status = configurationStatus_;
    invalid.phaseStatus = tpFlash_.configurationStatus();
    invalid.caloricStatus = caloric_.configurationStatus();
    if (configurationStatus_ != RhoEFlashStatus::SUCCESS) return invalid;

    if (!finitePositive(pressure))
    {
        invalid.status = RhoEFlashStatus::TP_FLASH_FAILURE;
        invalid.phaseStatus = PhaseEquilibriumStatus::INVALID_PRESSURE;
        return invalid;
    }
    if (pressure < options_.minimumPressure
        || pressure > options_.maximumPressure)
    {
        invalid.status = RhoEFlashStatus::OUT_OF_CONFIGURED_RANGE;
        invalid.phaseStatus = PhaseEquilibriumStatus::INVALID_PRESSURE;
        return invalid;
    }
    if (!finitePositive(temperature))
    {
        invalid.status = RhoEFlashStatus::TP_FLASH_FAILURE;
        invalid.phaseStatus = PhaseEquilibriumStatus::INVALID_TEMPERATURE;
        return invalid;
    }
    if (temperature < options_.minimumTemperature
        || temperature > options_.maximumTemperature)
    {
        invalid.status = RhoEFlashStatus::OUT_OF_CONFIGURED_RANGE;
        invalid.phaseStatus = PhaseEquilibriumStatus::INVALID_TEMPERATURE;
        return invalid;
    }

    std::vector<double> feedMoles;
    const PhaseEquilibriumStatus compositionStatus =
        tpFlash_.kernel().massToMoleFractions(
            overallMassFractions, feedMoles);
    if (compositionStatus != PhaseEquilibriumStatus::SUCCESS)
    {
        invalid.status = RhoEFlashStatus::INVALID_COMPOSITION;
        invalid.phaseStatus = compositionStatus;
        return invalid;
    }

    Evaluation evaluation = assemblePT(
        pressure, temperature, overallMassFractions);
    if (!evaluation.valid)
    {
        return commit(
            evaluation.density,
            evaluation.energy,
            evaluation,
            0,
            0,
            1);
    }

    evaluation.densityResidual = 0.0;
    evaluation.densityRelativeResidual = 0.0;
    evaluation.energyResidual = 0.0;
    evaluation.energyRelativeResidual = 0.0;
    return commit(
        evaluation.density,
        evaluation.energy,
        evaluation,
        0,
        1,
        0);
}

RhoEFlashResult PRRhoEFlash::flashMassFractions(
    double targetDensity,
    double targetEnergy,
    const std::vector<double>& overallMassFractions,
    const RhoEFlashSeed& seed) const
{
    RhoEFlashResult invalid;
    invalid.targetDensity = targetDensity;
    invalid.targetEnergy = targetEnergy;
    invalid.status = configurationStatus_;
    invalid.phaseStatus = tpFlash_.configurationStatus();
    invalid.caloricStatus = caloric_.configurationStatus();
    if (configurationStatus_ != RhoEFlashStatus::SUCCESS) return invalid;
    if (!finitePositive(targetDensity))
    {
        invalid.status = RhoEFlashStatus::INVALID_DENSITY;
        return invalid;
    }
    if (!std::isfinite(targetEnergy))
    {
        invalid.status = RhoEFlashStatus::NONFINITE_TARGET_ENERGY;
        return invalid;
    }
    std::vector<double> feedMoles;
    if (tpFlash_.kernel().massToMoleFractions(
            overallMassFractions, feedMoles)
            != PhaseEquilibriumStatus::SUCCESS)
    {
        invalid.status = RhoEFlashStatus::INVALID_COMPOSITION;
        return invalid;
    }

    std::size_t validSamples = 0;
    std::size_t invalidSamples = 0;
    Evaluation best;
    bool haveBest = false;
    Evaluation regionMismatch;
    bool haveRegionMismatch = false;
    std::size_t localIterations = 0;
    Evaluation bestSaturation;
    bool haveBestSaturation = false;
    auto considerBest = [&](const Evaluation& candidate)
    {
        const double candidateNorm = std::hypot(
            scaledDensityResidual(candidate.density, targetDensity),
            scaledEnergyResidual(candidate.energy, targetEnergy));
        const double bestNorm = haveBest
            ? std::hypot(
                scaledDensityResidual(best.density, targetDensity),
                scaledEnergyResidual(best.energy, targetEnergy))
            : std::numeric_limits<double>::infinity();
        if (std::isfinite(candidateNorm) && candidateNorm < bestNorm)
        {
            best = candidate;
            haveBest = true;
        }
    };

    const bool haveCallSeed = finitePositive(seed.pressure)
        && finitePositive(seed.temperature)
        && seed.pressure >= options_.minimumPressure
        && seed.pressure <= options_.maximumPressure
        && seed.temperature >= options_.minimumTemperature
        && seed.temperature <= options_.maximumTemperature;
    if (haveCallSeed && seed.localOnly && speciesNames_.size() == 1
        && (seed.requiredRegion == EquilibriumRegion::UNKNOWN
            || seed.requiredRegion
                == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE))
    {
        Evaluation seededSaturation = evaluatePureSaturationPressure(
            seed.pressure, targetDensity, targetEnergy,
            overallMassFractions);
        if (seededSaturation.valid)
        {
            ++validSamples;
            bestSaturation = std::move(seededSaturation);
            haveBestSaturation = true;
            const auto solved = refinePureSaturation(
                    targetDensity, targetEnergy, overallMassFractions,
                    bestSaturation, localIterations,
                    validSamples, invalidSamples);
            if (solved && seed.localOnly) return *solved;
            considerBest(bestSaturation);
        }
        else
        {
            ++invalidSamples;
        }
    }
    if (haveCallSeed && seed.localOnly
        && seed.requiredRegion
            != EquilibriumRegion::PURE_COMPONENT_COEXISTENCE)
    {
        const auto solved = solveGeneralSeed(
                std::log(seed.pressure), std::log(seed.temperature),
                seed.localOnly,
                seed.requiredRegion,
                targetDensity, targetEnergy, overallMassFractions,
                best, haveBest, regionMismatch, haveRegionMismatch,
                localIterations,
                validSamples, invalidSamples);
        if (solved && seed.localOnly) return *solved;
    }
    if (seed.localOnly)
    {
        if (haveRegionMismatch)
        {
            regionMismatch.status =
                RhoEFlashStatus::REQUIRED_REGION_MISMATCH;
            regionMismatch.valid = false;
            return commit(
                targetDensity, targetEnergy, regionMismatch,
                localIterations, validSamples, invalidSamples);
        }
        if (haveBest)
        {
            best.status = RhoEFlashStatus::NOT_CONVERGED;
            return commit(
                targetDensity, targetEnergy, best,
                localIterations, validSamples, invalidSamples);
        }
        invalid.status = RhoEFlashStatus::NO_RHOE_BRACKET;
        invalid.validSamples = validSamples;
        invalid.invalidSamples = invalidSamples;
        return invalid;
    }

    // A pure two-phase state has a one-dimensional saturation manifold. Scan
    // pressure first and use both the volume and energy lever rules; this
    // avoids asking TP to invent a phase fraction.
    if (speciesNames_.size() == 1)
    {
        bool havePrevious = false;
        Evaluation previous;
        double previousLogPressure = 0.0;
        for (std::size_t sample = 0;
             sample <= options_.pressureScanIntervals; ++sample)
        {
            const double fraction = static_cast<double>(sample)
                /static_cast<double>(options_.pressureScanIntervals);
            const double logPressure = std::log(options_.minimumPressure)
                + fraction*(std::log(options_.maximumPressure)
                    - std::log(options_.minimumPressure));
            const double pressure = std::exp(logPressure);
            Evaluation current = evaluatePureSaturationPressure(
                pressure, targetDensity, targetEnergy, overallMassFractions);
            if (!current.valid)
            {
                ++invalidSamples;
                havePrevious = false;
                continue;
            }
            ++validSamples;
            if (!haveBestSaturation || std::abs(current.energyResidual)
                    < std::abs(bestSaturation.energyResidual))
            {
                bestSaturation = current;
                haveBestSaturation = true;
            }
            if (densityConverged(current.density, targetDensity, options_)
                && energyConverged(current.energy, targetEnergy, options_))
            {
                return commit(
                    targetDensity, targetEnergy, current, 0,
                    validSamples, invalidSamples);
            }
            if (havePrevious
                && ((previous.energyResidual < 0.0
                        && current.energyResidual > 0.0)
                    || (previous.energyResidual > 0.0
                        && current.energyResidual < 0.0)))
            {
                Evaluation lower = previous;
                Evaluation upper = current;
                double lowerLogPressure = previousLogPressure;
                double upperLogPressure = logPressure;
                for (std::size_t iteration = 1;
                     iteration <= options_.maximumIterations; ++iteration)
                {
                    const double middleLogPressure =
                        0.5*(lowerLogPressure + upperLogPressure);
                    Evaluation middle = evaluatePureSaturationPressure(
                        std::exp(middleLogPressure),
                        targetDensity,
                        targetEnergy,
                        overallMassFractions);
                    if (!middle.valid)
                    {
                        ++invalidSamples;
                        lower = Evaluation();
                        break;
                    }
                    ++validSamples;
                    if (!haveBestSaturation
                        || std::abs(middle.energyResidual)
                            < std::abs(bestSaturation.energyResidual))
                    {
                        bestSaturation = middle;
                        haveBestSaturation = true;
                    }
                    if (densityConverged(
                            middle.density, targetDensity, options_)
                        && energyConverged(
                            middle.energy, targetEnergy, options_))
                    {
                        return commit(
                            targetDensity, targetEnergy, middle, iteration,
                            validSamples, invalidSamples);
                    }
                    if (upperLogPressure - lowerLogPressure
                        <= options_.saturationLogPressureTolerance)
                    {
                        break;
                    }
                    if ((lower.energyResidual < 0.0
                            && middle.energyResidual > 0.0)
                        || (lower.energyResidual > 0.0
                            && middle.energyResidual < 0.0))
                    {
                        upper = std::move(middle);
                        upperLogPressure = middleLogPressure;
                    }
                    else
                    {
                        lower = std::move(middle);
                        lowerLogPressure = middleLogPressure;
                    }
                }
            }
            havePrevious = true;
            previous = std::move(current);
            previousLogPressure = logPressure;
        }

        // Near beta=0 or beta=1 the valid saturation-pressure interval can be
        // narrower than one scan cell. Refine the best valid point without
        // requiring a sampled sign change.
        if (haveBestSaturation)
        {
            if (const auto solved = refinePureSaturation(
                    targetDensity, targetEnergy, overallMassFractions,
                    bestSaturation, localIterations,
                    validSamples, invalidSamples))
                return *solved;
            considerBest(bestSaturation);
        }
    }

    // General bounded Newton solve in log(p), log(T). Invalid thermodynamic
    // evaluations are never replaced by endpoint values or an ideal fallback.
    const double logPMin = std::log(options_.minimumPressure);
    const double logPMax = std::log(options_.maximumPressure);
    const double logTMin = std::log(options_.minimumTemperature);
    const double logTMax = std::log(options_.maximumTemperature);
    std::vector<std::pair<double, double>> seeds;
    const double seedP = std::isfinite(options_.initialPressure)
        ? options_.initialPressure
        : std::sqrt(options_.minimumPressure*options_.maximumPressure);
    const double seedT = std::isfinite(options_.initialTemperature)
        ? options_.initialTemperature
        : std::sqrt(options_.minimumTemperature*options_.maximumTemperature);
    if (finitePositive(seedP) && finitePositive(seedT))
        seeds.emplace_back(std::log(seedP), std::log(seedT));
    seeds.emplace_back(0.5*(logPMin + logPMax), 0.5*(logTMin + logTMax));
    if (speciesNames_.size() == 1 && haveBestSaturation
        && finitePositive(bestSaturation.pressure)
        && finitePositive(bestSaturation.temperature))
    {
        const double offset = std::max(
            4.0*options_.finiteDifferenceStep, 1.0e-4);
        const double bestLogPressure = std::log(bestSaturation.pressure);
        const double bestLogTemperature = std::log(bestSaturation.temperature);
        seeds.emplace_back(
            bestLogPressure,
            std::max(logTMin, bestLogTemperature - offset));
        seeds.emplace_back(
            bestLogPressure,
            std::min(logTMax, bestLogTemperature + offset));
    }

    for (const auto& candidateSeed : seeds)
    {
        if (const auto solved = solveGeneralSeed(
                candidateSeed.first, candidateSeed.second,
                false,
                EquilibriumRegion::UNKNOWN,
                targetDensity, targetEnergy, overallMassFractions,
                best, haveBest, regionMismatch, haveRegionMismatch,
                localIterations,
                validSamples, invalidSamples))
            return *solved;
    }

    if (haveBest)
    {
        best.status = RhoEFlashStatus::NOT_CONVERGED;
        return commit(
            targetDensity, targetEnergy, best, options_.maximumIterations,
            validSamples, invalidSamples);
    }
    invalid.status = RhoEFlashStatus::NO_RHOE_BRACKET;
    invalid.validSamples = validSamples;
    invalid.invalidSamples = invalidSamples;
    return invalid;
}

} // End namespace phaseEquilibrium
} // End namespace Foam
