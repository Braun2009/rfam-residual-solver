#include "PRPhaseSoundSpeed.H"

#include <algorithm>
#include <cmath>
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

bool closeScalar(double lhs, double rhs, double relativeTolerance)
{
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) return false;
    return std::abs(lhs - rhs) <= relativeTolerance
        *std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

bool sameComposition(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs,
    double absoluteTolerance)
{
    if (lhs.size() != rhs.size() || lhs.empty()) return false;
    for (std::size_t species = 0; species < lhs.size(); ++species)
    {
        if (!std::isfinite(lhs[species]) || !std::isfinite(rhs[species])
            || std::abs(lhs[species] - rhs[species]) > absoluteTolerance)
            return false;
    }
    return true;
}

} // End anonymous namespace

const char* phaseSoundSpeedStatusName(PhaseSoundSpeedStatus status)
{
    switch (status)
    {
        case PhaseSoundSpeedStatus::SUCCESS: return "SUCCESS";
        case PhaseSoundSpeedStatus::INVALID_CONFIGURATION:
            return "INVALID_CONFIGURATION";
        case PhaseSoundSpeedStatus::INVALID_INPUT: return "INVALID_INPUT";
        case PhaseSoundSpeedStatus::BASE_STATE_INVALID:
            return "BASE_STATE_INVALID";
        case PhaseSoundSpeedStatus::BASE_STATE_MISMATCH:
            return "BASE_STATE_MISMATCH";
        case PhaseSoundSpeedStatus::NONFINITE_ENTROPY:
            return "NONFINITE_ENTROPY";
        case PhaseSoundSpeedStatus::NO_ENTROPY_BRACKET:
            return "NO_ENTROPY_BRACKET";
        case PhaseSoundSpeedStatus::FLASH_FAILURE:
            return "FLASH_FAILURE";
        case PhaseSoundSpeedStatus::NOT_CONVERGED:
            return "NOT_CONVERGED";
        case PhaseSoundSpeedStatus::PHASE_TOPOLOGY_CHANGE:
            return "PHASE_TOPOLOGY_CHANGE";
        case PhaseSoundSpeedStatus::INSUFFICIENT_FLASH_ACCURACY:
            return "INSUFFICIENT_FLASH_ACCURACY";
        case PhaseSoundSpeedStatus::PERTURBATION_NOT_CONVERGED:
            return "PERTURBATION_NOT_CONVERGED";
        case PhaseSoundSpeedStatus::NONPOSITIVE_DERIVATIVE:
            return "NONPOSITIVE_DERIVATIVE";
    }
    return "UNKNOWN_STATUS";
}

PRPhaseSoundSpeed::PRPhaseSoundSpeed(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    const RhoEFlashOptions& rhoEOptions,
    const PhaseSoundSpeedOptions& options)
:
    rhoEFlash_(criticalData, kij, speciesNames, rhoEOptions),
    options_(options),
    configurationStatus_(PhaseSoundSpeedStatus::SUCCESS)
{
    if (rhoEFlash_.configurationStatus() != RhoEFlashStatus::SUCCESS
        || !finitePositive(options_.relativeDensityPerturbation)
        || options_.relativeDensityPerturbation >= 0.25
        || !finitePositive(options_.minimumRelativeDensityPerturbation)
        || options_.minimumRelativeDensityPerturbation
            > options_.relativeDensityPerturbation
        || options_.entropyScanIntervals == 0
        || options_.maximumIterations == 0
        || options_.minimumDerivativeConvergenceComparisons == 0
        || options_.minimumDerivativeConvergenceComparisons
            > options_.maximumTopologyRefinements
        || !finitePositive(options_.energyCorrectionSpanFactor)
        || !std::isfinite(options_.entropyAbsoluteTolerance)
        || options_.entropyAbsoluteTolerance < 0.0
        || !finitePositive(options_.energyIntervalRelativeTolerance)
        || !finitePositive(options_.maximumFlashDensityErrorFraction)
        || options_.maximumFlashDensityErrorFraction >= 0.25
        || !finitePositive(options_.derivativeRelativeTolerance)
        || options_.derivativeRelativeTolerance > 1.0e-2)
        configurationStatus_ = PhaseSoundSpeedStatus::INVALID_CONFIGURATION;
}

bool PRPhaseSoundSpeed::entropyConverged(double value, double target) const
{
    if (!std::isfinite(value) || !std::isfinite(target)) return false;
    return std::abs(value - target) <= options_.entropyAbsoluteTolerance;
}

PRPhaseSoundSpeed::DensitySolve PRPhaseSoundSpeed::solveAtDensity(
    double targetDensity,
    double targetEntropy,
    double energySeed,
    double pressureSeed,
    double temperatureSeed,
    EquilibriumRegion requiredRegion,
    const std::vector<double>& massFractions) const
{
    DensitySolve result;
    if (!finitePositive(targetDensity) || !std::isfinite(targetEntropy)
        || !std::isfinite(energySeed))
    {
        result.status = PhaseSoundSpeedStatus::INVALID_INPUT;
        return result;
    }

    auto oppositeSigns = [](double lhs, double rhs)
    {
        return (lhs < 0.0 && rhs > 0.0)
            || (lhs > 0.0 && rhs < 0.0);
    };

    auto failureStatus = [](const RhoEFlashResult& state)
    {
        return state.status == RhoEFlashStatus::REQUIRED_REGION_MISMATCH
            ? PhaseSoundSpeedStatus::PHASE_TOPOLOGY_CHANGE
            : PhaseSoundSpeedStatus::FLASH_FAILURE;
    };

    auto solveBracket = [&](double firstEnergy,
                            const RhoEFlashResult& firstState,
                            double firstResidual,
                            double secondEnergy,
                            const RhoEFlashResult& secondState,
                            double secondResidual)
    {
        DensitySolve bracket = result;
        double lowerEnergy = firstEnergy;
        double upperEnergy = secondEnergy;
        double lowerResidual = firstResidual;
        double upperResidual = secondResidual;
        RhoEFlashResult lowerState = firstState;
        RhoEFlashResult upperState = secondState;
        if (lowerEnergy > upperEnergy)
        {
            std::swap(lowerEnergy, upperEnergy);
            std::swap(lowerResidual, upperResidual);
            std::swap(lowerState, upperState);
        }

        if (!((lowerResidual < 0.0 && upperResidual > 0.0)
              || (lowerResidual > 0.0 && upperResidual < 0.0)))
        {
            bracket.status = PhaseSoundSpeedStatus::NO_ENTROPY_BRACKET;
            return bracket;
        }

        for (std::size_t iteration = 1;
             iteration <= options_.maximumIterations; ++iteration)
        {
            const double middleEnergy = 0.5*(lowerEnergy + upperEnergy);
            const RhoEFlashResult& warmState =
                std::abs(middleEnergy - lowerEnergy)
                    <= std::abs(upperEnergy - middleEnergy)
                ? lowerState : upperState;
            const RhoEFlashResult middle = rhoEFlash_.flashMassFractions(
                targetDensity, middleEnergy, massFractions,
                {warmState.pressure, warmState.temperature,
                 true, requiredRegion});
            if (!middle.success() || !std::isfinite(middle.entropy))
            {
                ++bracket.invalidSamples;
                bracket.state = middle;
                bracket.iterations = iteration;
                bracket.status = failureStatus(middle);
                return bracket;
            }
            ++bracket.validSamples;
            const double middleResidual = middle.entropy - targetEntropy;
            if (entropyConverged(middle.entropy, targetEntropy))
            {
                bracket.status = PhaseSoundSpeedStatus::SUCCESS;
                bracket.state = middle;
                bracket.entropyResidual = middleResidual;
                bracket.iterations = iteration;
                return bracket;
            }
            if (std::abs(upperEnergy - lowerEnergy)
                <= options_.energyIntervalRelativeTolerance
                    *std::max(std::abs(secondEnergy - firstEnergy), 1.0))
            {
                bracket.state = middle;
                bracket.entropyResidual = middleResidual;
                bracket.iterations = iteration;
                bracket.status = PhaseSoundSpeedStatus::NOT_CONVERGED;
                return bracket;
            }
            if (oppositeSigns(lowerResidual, middleResidual))
            {
                upperEnergy = middleEnergy;
                upperResidual = middleResidual;
                upperState = middle;
            }
            else
            {
                lowerEnergy = middleEnergy;
                lowerResidual = middleResidual;
                lowerState = middle;
            }
        }

        bracket.state = std::abs(lowerResidual) <= std::abs(upperResidual)
            ? lowerState : upperState;
        bracket.entropyResidual = std::abs(lowerResidual)
            <= std::abs(upperResidual) ? lowerResidual : upperResidual;
        bracket.iterations = options_.maximumIterations;
        bracket.status = PhaseSoundSpeedStatus::NOT_CONVERGED;
        return bracket;
    };

    // For an isentropic density perturbation, the first law gives
    // de = -p dv. The caller supplies that first-order prediction as the seed,
    // so search only the connected valid interval around that seed.
    RhoEFlashResult seedState = rhoEFlash_.flashMassFractions(
        targetDensity, energySeed, massFractions,
        {pressureSeed, temperatureSeed, true, requiredRegion});
    if (!seedState.success() || !std::isfinite(seedState.entropy))
    {
        result.state = seedState;
        result.invalidSamples = 1;
        result.status = failureStatus(seedState);
        return result;
    }

    ++result.validSamples;
    double previousEnergy = energySeed;
    double previousResidual = seedState.entropy - targetEntropy;
    RhoEFlashResult previousState = seedState;
    if (entropyConverged(seedState.entropy, targetEntropy))
    {
        result.status = PhaseSoundSpeedStatus::SUCCESS;
        result.state = seedState;
        result.entropyResidual = previousResidual;
        return result;
    }

    const double temperatureScale = std::max(std::abs(seedState.temperature), 1.0);
    const double entropyCorrection =
        temperatureScale*std::abs(previousResidual);
    const double toleranceCorrection = temperatureScale
        *std::max(options_.entropyAbsoluteTolerance, 1.0e-12);
    const double correctionScale = std::max(
        {entropyCorrection, toleranceCorrection, 1.0e-12});
    const double roundoffStep = 32.0*std::numeric_limits<double>::epsilon()
        *std::max(std::abs(energySeed), 1.0);
    const double minimumStep = std::max(
        roundoffStep, 1.0e-6*correctionScale);
    double step = std::max(1.25*entropyCorrection, minimumStep);
    const double maximumSpan =
        options_.energyCorrectionSpanFactor*correctionScale;
    const double direction = previousResidual > 0.0 ? -1.0 : 1.0;

    bool haveBarrier = false;
    double barrierEnergy = std::numeric_limits<double>::quiet_NaN();
    RhoEFlashResult barrierState;
    std::size_t probesUsed = 0;
    for (std::size_t probe = 0;
         probe < options_.entropyScanIntervals; ++probe)
    {
        const double travelled = std::abs(previousEnergy - energySeed);
        double trialStep = step;
        if (haveBarrier)
        {
            const double barrierDistance =
                std::abs(barrierEnergy - previousEnergy);
            if (barrierDistance <= minimumStep)
            {
                result.state = barrierState;
                result.entropyResidual = previousResidual;
                result.iterations = probe;
                result.status = PhaseSoundSpeedStatus::FLASH_FAILURE;
                return result;
            }
            trialStep = 0.5*barrierDistance;
        }
        else
        {
            const double remaining = maximumSpan - travelled;
            if (!(remaining > 0.0)) break;
            trialStep = std::min(trialStep, remaining);
        }
        if (!(trialStep > 0.0) || !std::isfinite(trialStep)) break;

        double candidateEnergy = previousEnergy + direction*trialStep;
        if (!std::isfinite(candidateEnergy)) break;
        probesUsed = probe + 1;
        RhoEFlashResult candidate = rhoEFlash_.flashMassFractions(
            targetDensity, candidateEnergy, massFractions,
            {previousState.pressure, previousState.temperature,
             true, requiredRegion});
        bool candidateValid =
            candidate.success() && std::isfinite(candidate.entropy);
        if (!candidateValid)
        {
            ++result.invalidSamples;
            if (candidate.status
                == RhoEFlashStatus::REQUIRED_REGION_MISMATCH)
            {
                result.state = candidate;
                result.entropyResidual = previousResidual;
                result.iterations = probe + 1;
                result.status =
                    PhaseSoundSpeedStatus::PHASE_TOPOLOGY_CHANGE;
                return result;
            }
            result.state = candidate;
            haveBarrier = true;
            barrierEnergy = candidateEnergy;
            barrierState = candidate;

            for (std::size_t refinement = 0; refinement < 32; ++refinement)
            {
                const double barrierDistance =
                    std::abs(barrierEnergy - previousEnergy);
                if (barrierDistance <= minimumStep) break;
                candidateEnergy = 0.5*(previousEnergy + barrierEnergy);
                candidate = rhoEFlash_.flashMassFractions(
                    targetDensity, candidateEnergy, massFractions,
                    {previousState.pressure, previousState.temperature,
                     true, requiredRegion});
                candidateValid =
                    candidate.success() && std::isfinite(candidate.entropy);
                if (candidateValid) break;
                ++result.invalidSamples;
                if (candidate.status
                    == RhoEFlashStatus::REQUIRED_REGION_MISMATCH)
                {
                    result.state = candidate;
                    result.entropyResidual = previousResidual;
                    result.iterations = probe + 1;
                    result.status =
                        PhaseSoundSpeedStatus::PHASE_TOPOLOGY_CHANGE;
                    return result;
                }
                result.state = candidate;
                barrierEnergy = candidateEnergy;
                barrierState = candidate;
            }
            if (!candidateValid)
            {
                result.state = candidate;
                result.entropyResidual = previousResidual;
                result.iterations = probe + 1;
                result.status = PhaseSoundSpeedStatus::FLASH_FAILURE;
                return result;
            }
            trialStep = std::abs(candidateEnergy - previousEnergy);
        }

        ++result.validSamples;
        const double residual = candidate.entropy - targetEntropy;
        if (entropyConverged(candidate.entropy, targetEntropy))
        {
            result.status = PhaseSoundSpeedStatus::SUCCESS;
            result.state = candidate;
            result.entropyResidual = residual;
            result.iterations = probe + 1;
            return result;
        }
        if (oppositeSigns(previousResidual, residual))
            return solveBracket(
                previousEnergy, previousState, previousResidual,
                candidateEnergy, candidate, residual);

        previousEnergy = candidateEnergy;
        previousResidual = residual;
        previousState = candidate;
        if (!haveBarrier)
            step = std::max(2.0*trialStep, minimumStep);
    }

    result.state = haveBarrier ? barrierState : previousState;
    result.entropyResidual = previousResidual;
    result.iterations = probesUsed;
    result.status = haveBarrier
        ? PhaseSoundSpeedStatus::FLASH_FAILURE
        : PhaseSoundSpeedStatus::NO_ENTROPY_BRACKET;
    return result;
}

PhaseSoundSpeedResult PRPhaseSoundSpeed::evaluateVerified(
    const RhoEFlashResult& base,
    const std::vector<double>& overallMassFractions) const
{
    PhaseSoundSpeedResult result;
    result.baseFlashStatus = base.status;
    result.baseRegion = base.tpFlash.region;
    result.targetDensity = base.targetDensity;
    result.targetEnergy = base.targetEnergy;
    result.requestedRelativeDensityPerturbation =
        options_.relativeDensityPerturbation;
    if (configurationStatus_ != PhaseSoundSpeedStatus::SUCCESS)
    {
        result.status = configurationStatus_;
        return result;
    }
    if (!base.success() || !finitePositive(base.pressure)
        || !finitePositive(base.temperature) || !finitePositive(base.density)
        || !std::isfinite(base.energy))
    {
        result.status = PhaseSoundSpeedStatus::BASE_STATE_INVALID;
        return result;
    }
    if (!std::isfinite(base.entropy))
    {
        result.status = PhaseSoundSpeedStatus::NONFINITE_ENTROPY;
        return result;
    }
    if (overallMassFractions.empty())
    {
        result.status = PhaseSoundSpeedStatus::INVALID_INPUT;
        return result;
    }
    result.baseState = base;
    result.targetEntropy = base.entropy;
    double pureEndpointTolerance = 1.0e-12;
    if (base.liquidPhase && base.vaporPhase)
    {
        const double phaseVolumeSpan =
            1.0/base.vaporPhase->rho - 1.0/base.liquidPhase->rho;
        const double acceptedDensityError =
            rhoEFlash_.options().densityAbsoluteTolerance
            + rhoEFlash_.options().densityRelativeTolerance
                *std::max(base.density, 1.0);
        if (finitePositive(phaseVolumeSpan)
            && finitePositive(acceptedDensityError))
        {
            pureEndpointTolerance = std::max(
                pureEndpointTolerance,
                4.0*acceptedDensityError
                    /(base.density*base.density*phaseVolumeSpan));
        }
    }
    if (base.tpFlash.region
            == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE
        && std::isfinite(base.vaporMoleFraction)
        && (base.vaporMoleFraction <= pureEndpointTolerance
            || base.vaporMoleFraction >= 1.0 - pureEndpointTolerance))
    {
        result.usedRelativeDensityPerturbation =
            options_.relativeDensityPerturbation;
        result.status = PhaseSoundSpeedStatus::PHASE_TOPOLOGY_CHANGE;
        return result;
    }

    double delta = options_.relativeDensityPerturbation;
    std::size_t topologyRefinements = 0;
    std::size_t derivativeComparisons = 0;
    std::size_t consecutiveConvergedComparisons = 0;
    std::optional<PhaseSoundSpeedResult> previousDerivative;
    auto stencilCanRefine = [](PhaseSoundSpeedStatus status)
    {
        return status == PhaseSoundSpeedStatus::FLASH_FAILURE
            || status == PhaseSoundSpeedStatus::NO_ENTROPY_BRACKET
            || status == PhaseSoundSpeedStatus::NOT_CONVERGED
            || status == PhaseSoundSpeedStatus::PHASE_TOPOLOGY_CHANGE;
    };
    for (std::size_t topologyRefinement = 0;
         topologyRefinement <= options_.maximumTopologyRefinements;
         ++topologyRefinement)
    {
        PhaseSoundSpeedResult attempt = result;
        attempt.usedRelativeDensityPerturbation = delta;
        attempt.stencilRefinements = topologyRefinement;
        attempt.topologyRefinements = topologyRefinements;
        const double lowerDensity = base.density*(1.0 - delta);
        const double upperDensity = base.density*(1.0 + delta);
        if (!finitePositive(lowerDensity) || !finitePositive(upperDensity))
        {
            attempt.status = PhaseSoundSpeedStatus::INVALID_INPUT;
            return attempt;
        }

        const double baseSpecificVolume = 1.0/base.density;
        const double lowerEnergySeed = base.energy + base.pressure
            *(baseSpecificVolume - 1.0/lowerDensity);
        const double upperEnergySeed = base.energy + base.pressure
            *(baseSpecificVolume - 1.0/upperDensity);
        const DensitySolve lower = solveAtDensity(
            lowerDensity, base.entropy, lowerEnergySeed,
            base.pressure, base.temperature, base.tpFlash.region,
            overallMassFractions);
        attempt.lowerFlashStatus = lower.state.status;
        attempt.lowerIterations = lower.iterations;
        attempt.lowerValidSamples = lower.validSamples;
        attempt.lowerInvalidSamples = lower.invalidSamples;
        attempt.lowerEntropyResidual = lower.entropyResidual;
        if (lower.state.tpFlash.valid())
            attempt.lowerRegion = lower.state.tpFlash.region;
        if (lower.status != PhaseSoundSpeedStatus::SUCCESS)
        {
            if (lower.status == PhaseSoundSpeedStatus::PHASE_TOPOLOGY_CHANGE)
            {
                ++topologyRefinements;
                attempt.topologyRefinements = topologyRefinements;
            }
            const double nextDelta = 0.5*delta;
            if (stencilCanRefine(lower.status)
                && topologyRefinement
                    < options_.maximumTopologyRefinements
                && nextDelta
                    >= options_.minimumRelativeDensityPerturbation)
            {
                delta = nextDelta;
                continue;
            }
            attempt.status = lower.status;
            return attempt;
        }
        const DensitySolve upper = solveAtDensity(
            upperDensity, base.entropy, upperEnergySeed,
            base.pressure, base.temperature, base.tpFlash.region,
            overallMassFractions);
        attempt.upperFlashStatus = upper.state.status;
        attempt.upperIterations = upper.iterations;
        attempt.upperValidSamples = upper.validSamples;
        attempt.upperInvalidSamples = upper.invalidSamples;
        attempt.upperEntropyResidual = upper.entropyResidual;
        if (upper.state.tpFlash.valid())
            attempt.upperRegion = upper.state.tpFlash.region;
        if (upper.status != PhaseSoundSpeedStatus::SUCCESS)
        {
            if (upper.status == PhaseSoundSpeedStatus::PHASE_TOPOLOGY_CHANGE)
            {
                ++topologyRefinements;
                attempt.topologyRefinements = topologyRefinements;
            }
            const double nextDelta = 0.5*delta;
            if (stencilCanRefine(upper.status)
                && topologyRefinement
                    < options_.maximumTopologyRefinements
                && nextDelta
                    >= options_.minimumRelativeDensityPerturbation)
            {
                delta = nextDelta;
                continue;
            }
            attempt.status = upper.status;
            return attempt;
        }

        attempt.lowerState = lower.state;
        attempt.upperState = upper.state;
        attempt.lowerDensity = lower.state.density;
        attempt.upperDensity = upper.state.density;
        attempt.lowerPressure = lower.state.pressure;
        attempt.upperPressure = upper.state.pressure;
        attempt.lowerEntropy = lower.state.entropy;
        attempt.upperEntropy = upper.state.entropy;

        if (attempt.lowerRegion != attempt.baseRegion
            || attempt.upperRegion != attempt.baseRegion)
        {
            ++topologyRefinements;
            attempt.topologyRefinements = topologyRefinements;
            const double nextDelta = 0.5*delta;
            if (topologyRefinement == options_.maximumTopologyRefinements
                || nextDelta < options_.minimumRelativeDensityPerturbation)
            {
                attempt.status = PhaseSoundSpeedStatus::PHASE_TOPOLOGY_CHANGE;
                return attempt;
            }
            delta = nextDelta;
            continue;
        }

        const double densityHalfSpacing = delta*base.density;
        const double maximumDensityError = std::max(
            std::abs(lower.state.density - lowerDensity),
            std::abs(upper.state.density - upperDensity));
        if (!std::isfinite(maximumDensityError)
            || maximumDensityError
                > options_.maximumFlashDensityErrorFraction
                    *densityHalfSpacing)
        {
            attempt.status =
                PhaseSoundSpeedStatus::INSUFFICIENT_FLASH_ACCURACY;
            return attempt;
        }

        attempt.densityDerivative =
            (attempt.upperPressure - attempt.lowerPressure)
                /(attempt.upperDensity - attempt.lowerDensity);
        attempt.soundSpeedSquared = attempt.densityDerivative;
        if (!std::isfinite(attempt.densityDerivative)
            || attempt.densityDerivative <= 0.0)
        {
            attempt.status = PhaseSoundSpeedStatus::NONPOSITIVE_DERIVATIVE;
            return attempt;
        }
        attempt.soundSpeed = std::sqrt(attempt.soundSpeedSquared);
        if (previousDerivative.has_value())
        {
            ++derivativeComparisons;
            attempt.derivativeConvergenceComparisons = derivativeComparisons;
            attempt.derivativeRelativeChange = std::abs(
                attempt.densityDerivative
                    - previousDerivative->densityDerivative)
                /std::max({
                    std::abs(attempt.densityDerivative),
                    std::abs(previousDerivative->densityDerivative),
                    1.0});
            if (attempt.derivativeRelativeChange
                <= options_.derivativeRelativeTolerance)
                ++consecutiveConvergedComparisons;
            else
                consecutiveConvergedComparisons = 0;

            if (consecutiveConvergedComparisons
                >= options_.minimumDerivativeConvergenceComparisons)
            {
                attempt.status = PhaseSoundSpeedStatus::SUCCESS;
                return attempt;
            }
        }

        previousDerivative = attempt;
        const double nextDelta = 0.5*delta;
        if (topologyRefinement == options_.maximumTopologyRefinements
            || nextDelta < options_.minimumRelativeDensityPerturbation)
        {
            attempt.status =
                PhaseSoundSpeedStatus::PERTURBATION_NOT_CONVERGED;
            return attempt;
        }
        delta = nextDelta;
    }

    result.status = previousDerivative.has_value()
        ? PhaseSoundSpeedStatus::PERTURBATION_NOT_CONVERGED
        : PhaseSoundSpeedStatus::PHASE_TOPOLOGY_CHANGE;
    return result;
}

PhaseSoundSpeedResult PRPhaseSoundSpeed::evaluate(
    const RhoEFlashResult& base,
    const std::vector<double>& overallMassFractions) const
{
    PhaseSoundSpeedResult result;
    result.baseFlashStatus = base.status;
    result.baseRegion = base.tpFlash.region;
    result.targetDensity = base.targetDensity;
    result.targetEnergy = base.targetEnergy;
    result.requestedRelativeDensityPerturbation =
        options_.relativeDensityPerturbation;
    if (configurationStatus_ != PhaseSoundSpeedStatus::SUCCESS)
    {
        result.status = configurationStatus_;
        return result;
    }
    if (!base.success() || !finitePositive(base.pressure)
        || !finitePositive(base.temperature) || !finitePositive(base.density)
        || !std::isfinite(base.energy))
    {
        result.status = PhaseSoundSpeedStatus::BASE_STATE_INVALID;
        return result;
    }
    if (!std::isfinite(base.entropy))
    {
        result.status = PhaseSoundSpeedStatus::NONFINITE_ENTROPY;
        return result;
    }
    if (overallMassFractions.empty())
    {
        result.status = PhaseSoundSpeedStatus::INVALID_INPUT;
        return result;
    }

    const double compositionTolerance = std::max(
        rhoEFlash_.options().tpOptions.materialTolerance, 1.0e-12);
    if (!sameComposition(
            base.tpFlash.feedMassFractions,
            overallMassFractions,
            compositionTolerance))
    {
        result.status = PhaseSoundSpeedStatus::BASE_STATE_MISMATCH;
        return result;
    }

    const RhoEFlashResult verified = rhoEFlash_.flashMassFractions(
        base.targetDensity, base.targetEnergy, overallMassFractions,
        {base.pressure, base.temperature, true, base.tpFlash.region});
    result.baseFlashStatus = verified.status;
    if (!verified.success()
        || verified.tpFlash.region != base.tpFlash.region
        || !closeScalar(verified.pressure, base.pressure, 1.0e-9)
        || !closeScalar(verified.temperature, base.temperature, 1.0e-9)
        || !closeScalar(verified.density, base.density, 1.0e-9)
        || !closeScalar(verified.energy, base.energy, 1.0e-9)
        || !closeScalar(verified.entropy, base.entropy, 1.0e-9))
    {
        result.status = PhaseSoundSpeedStatus::BASE_STATE_MISMATCH;
        return result;
    }
    return evaluateVerified(verified, overallMassFractions);
}

PhaseSoundSpeedResult PRPhaseSoundSpeed::evaluateRhoE(
    double targetDensity,
    double targetEnergy,
    const std::vector<double>& overallMassFractions) const
{
    PhaseSoundSpeedResult result;
    const RhoEFlashResult base = rhoEFlash_.flashMassFractions(
        targetDensity, targetEnergy, overallMassFractions);
    result = evaluateVerified(base, overallMassFractions);
    result.targetDensity = targetDensity;
    result.targetEnergy = targetEnergy;
    return result;
}

} // End namespace phaseEquilibrium
} // End namespace Foam
