#include "PRTPFlash.H"

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

double maximumAbsoluteMaterialResidual(
    const std::vector<double>& feed,
    const std::vector<double>& liquid,
    const std::vector<double>& vapor,
    double vaporMoleFraction)
{
    double maximumResidual = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i)
    {
        const double reconstructed =
            (1.0 - vaporMoleFraction)*liquid[i]
            + vaporMoleFraction*vapor[i];
        maximumResidual = std::max(
            maximumResidual, std::abs(feed[i] - reconstructed));
    }
    return maximumResidual;
}

bool solveDenseLinearSystem(
    std::vector<std::vector<double>> matrix,
    std::vector<double> rhs,
    std::vector<double>& solution)
{
    const std::size_t size = rhs.size();
    if (size == 0 || matrix.size() != size) return false;

    double matrixScale = 0.0;
    for (const std::vector<double>& row : matrix)
    {
        if (row.size() != size) return false;
        for (double value : row)
        {
            if (!std::isfinite(value)) return false;
            matrixScale = std::max(matrixScale, std::abs(value));
        }
    }
    for (double value : rhs)
    {
        if (!std::isfinite(value)) return false;
    }

    const double pivotTolerance = 64.0
        *std::numeric_limits<double>::epsilon()
        *std::max(1.0, matrixScale)*static_cast<double>(size);
    for (std::size_t column = 0; column < size; ++column)
    {
        std::size_t pivotRow = column;
        for (std::size_t row = column + 1; row < size; ++row)
        {
            if (std::abs(matrix[row][column])
                > std::abs(matrix[pivotRow][column]))
            {
                pivotRow = row;
            }
        }
        if (std::abs(matrix[pivotRow][column]) <= pivotTolerance)
            return false;

        if (pivotRow != column)
        {
            std::swap(matrix[pivotRow], matrix[column]);
            std::swap(rhs[pivotRow], rhs[column]);
        }

        const double pivot = matrix[column][column];
        for (std::size_t entry = column; entry < size; ++entry)
            matrix[column][entry] /= pivot;
        rhs[column] /= pivot;

        for (std::size_t row = 0; row < size; ++row)
        {
            if (row == column) continue;
            const double factor = matrix[row][column];
            for (std::size_t entry = column; entry < size; ++entry)
            {
                matrix[row][entry] -= factor*matrix[column][entry];
            }
            rhs[row] -= factor*rhs[column];
        }
    }

    for (double value : rhs)
    {
        if (!std::isfinite(value)) return false;
    }
    solution = std::move(rhs);
    return true;
}

} // End anonymous namespace

PRTPFlash::PRTPFlash(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    const TPFlashOptions& options)
:
    kernel_(criticalData, kij, speciesNames),
    options_(options),
    optionsStatus_(PhaseEquilibriumStatus::SUCCESS)
{
    if (options_.maxFlashIterations == 0
        || options_.maxStabilityIterations == 0
        || options_.maxRachfordRiceIterations == 0
        || options_.saturationScanIntervals < 2
        || options_.maxSaturationIterations == 0
        || !finitePositive(options_.fugacityTolerance)
        || !finitePositive(options_.materialTolerance)
        || !finitePositive(options_.stabilityTolerance)
        || !finitePositive(options_.iterationTolerance)
        || !finitePositive(options_.saturationTemperatureTolerance)
        || !finitePositive(options_.successiveSubstitutionDamping)
        || options_.successiveSubstitutionDamping > 1.0)
    {
        optionsStatus_ = PhaseEquilibriumStatus::INVALID_CONFIGURATION;
    }
}

PhaseEquilibriumStatus PRTPFlash::configurationStatus() const
{
    if (kernel_.configurationStatus() != PhaseEquilibriumStatus::SUCCESS)
        return kernel_.configurationStatus();
    return optionsStatus_;
}

RachfordRiceResult PRTPFlash::solveRachfordRice(
    const std::vector<double>& feedMoleFractions,
    const std::vector<double>& K) const
{
    RachfordRiceResult result;
    result.status = configurationStatus();
    if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;
    result.status = kernel_.validateMoleFractions(feedMoleFractions);
    if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;
    if (K.size() != feedMoleFractions.size())
    {
        result.status = PhaseEquilibriumStatus::INVALID_COMPOSITION_SIZE;
        return result;
    }
    for (double value : K)
    {
        if (!finitePositive(value))
        {
            result.status = PhaseEquilibriumStatus::WILSON_K_FAILURE;
            return result;
        }
    }

    const auto residual = [&](double vaporFraction, double& value)
    {
        value = 0.0;
        for (std::size_t i = 0; i < K.size(); ++i)
        {
            const double delta = K[i] - 1.0;
            const double denominator = 1.0 + vaporFraction*delta;
            if (!finitePositive(denominator)) return false;
            value += feedMoleFractions[i]*delta/denominator;
        }
        return std::isfinite(value);
    };

    double lowerResidual = 0.0;
    double upperResidual = 0.0;
    if (!residual(0.0, lowerResidual) || !residual(1.0, upperResidual))
    {
        result.status = PhaseEquilibriumStatus::NONFINITE_PARAMETER;
        return result;
    }

    const double tolerance = 0.1*options_.materialTolerance;
    if (std::abs(lowerResidual) <= tolerance)
    {
        result.status = PhaseEquilibriumStatus::SUCCESS;
        result.vaporMoleFraction = 0.0;
        result.residual = lowerResidual;
        return result;
    }
    if (std::abs(upperResidual) <= tolerance)
    {
        result.status = PhaseEquilibriumStatus::SUCCESS;
        result.vaporMoleFraction = 1.0;
        result.residual = upperResidual;
        return result;
    }
    if (!(lowerResidual > 0.0 && upperResidual < 0.0))
    {
        result.status = PhaseEquilibriumStatus::RACHFORD_RICE_NO_BRACKET;
        return result;
    }

    double lower = 0.0;
    double upper = 1.0;
    for (std::size_t iteration = 1;
         iteration <= options_.maxRachfordRiceIterations;
         ++iteration)
    {
        const double middle = 0.5*(lower + upper);
        double middleResidual = 0.0;
        if (!residual(middle, middleResidual))
        {
            result.status = PhaseEquilibriumStatus::NONFINITE_PARAMETER;
            return result;
        }
        result.iterations = iteration;
        result.vaporMoleFraction = middle;
        result.residual = middleResidual;

        if (std::abs(middleResidual) <= tolerance
            && upper - lower <= options_.iterationTolerance)
        {
            result.status = PhaseEquilibriumStatus::SUCCESS;
            return result;
        }
        if (middleResidual > 0.0)
            lower = middle;
        else
            upper = middle;
    }

    result.status = PhaseEquilibriumStatus::RACHFORD_RICE_NOT_CONVERGED;
    return result;
}

PRTPFlash::StabilityTrialResult PRTPFlash::stabilityTrial(
    double pressure,
    double temperature,
    const std::vector<double>& feed,
    const PRPhaseState& feedState,
    const std::vector<double>& initialK,
    bool vaporLike) const
{
    StabilityTrialResult result;
    if (!feedState.valid() || initialK.size() != feed.size())
    {
        result.status = PhaseEquilibriumStatus::INVALID_CONFIGURATION;
        return result;
    }

    std::vector<double> logAmounts(feed.size(), 0.0);
    std::vector<double> feedLogFugacity(feed.size(), 0.0);
    std::vector<std::size_t> activeIndices;
    activeIndices.reserve(feed.size());
    for (std::size_t i = 0; i < feed.size(); ++i)
    {
        if (feed[i] == 0.0) continue;
        const double value = vaporLike
            ? feed[i]*initialK[i]
            : feed[i]/initialK[i];
        if (!finitePositive(value))
        {
            result.status = PhaseEquilibriumStatus::NONFINITE_PARAMETER;
            return result;
        }
        logAmounts[i] = std::log(value);
        feedLogFugacity[i] = std::log(feed[i])
            + feedState.lnFugacityCoefficients[i];
        if (!std::isfinite(logAmounts[i])
            || !std::isfinite(feedLogFugacity[i]))
        {
            result.status = PhaseEquilibriumStatus::NONFINITE_PARAMETER;
            return result;
        }
        activeIndices.push_back(i);
    }
    if (activeIndices.empty())
    {
        result.status = PhaseEquilibriumStatus::INVALID_CONFIGURATION;
        return result;
    }

    const double maximumLog = std::log(std::numeric_limits<double>::max());
    const double minimumLog = std::log(std::numeric_limits<double>::min());
    double bestMaximumLogChange = std::numeric_limits<double>::infinity();
    std::size_t stagnantIterations = 0;
    bool usedExpandedSeedRestart = false;
    const auto evaluateResidual = [&]
    (
        const std::vector<double>& candidateLogs,
        std::vector<double>& composition,
        std::vector<double>& residual,
        double& amountSum
    )
    {
        composition.assign(feed.size(), 0.0);
        residual.assign(feed.size(), 0.0);
        amountSum = 0.0;
        for (std::size_t index : activeIndices)
        {
            const double logAmount = candidateLogs[index];
            if (!std::isfinite(logAmount)
                || logAmount >= maximumLog || logAmount <= minimumLog)
            {
                return PhaseEquilibriumStatus::NONFINITE_PARAMETER;
            }
            composition[index] = std::exp(logAmount);
            if (!finitePositive(composition[index]))
                return PhaseEquilibriumStatus::NONFINITE_PARAMETER;
            amountSum += composition[index];
            if (!finitePositive(amountSum))
                return PhaseEquilibriumStatus::NONFINITE_PARAMETER;
        }
        for (std::size_t index : activeIndices)
            composition[index] /= amountSum;

        const PRPhaseState trialState = kernel_.phaseState(
            pressure,
            temperature,
            composition,
            vaporLike ? PhaseRoot::VAPOR : PhaseRoot::LIQUID);
        if (!trialState.valid()) return trialState.status;

        for (std::size_t index : activeIndices)
        {
            residual[index] = feedLogFugacity[index]
                - trialState.lnFugacityCoefficients[index]
                - candidateLogs[index];
            if (!std::isfinite(residual[index]))
                return PhaseEquilibriumStatus::NONFINITE_PARAMETER;
        }
        return PhaseEquilibriumStatus::SUCCESS;
    };

    for (std::size_t iteration = 1;
         iteration <= options_.maxStabilityIterations;
         ++iteration)
    {
        std::vector<double> trial;
        std::vector<double> residual;
        double amountSum = 0.0;
        const PhaseEquilibriumStatus evaluationStatus = evaluateResidual(
            logAmounts, trial, residual, amountSum);
        if (evaluationStatus != PhaseEquilibriumStatus::SUCCESS)
        {
            result.status = evaluationStatus;
            return result;
        }

        double maximumLogChange = 0.0;
        for (std::size_t index : activeIndices)
        {
            maximumLogChange = std::max(
                maximumLogChange, std::abs(residual[index]));
        }
        result.iterations = iteration;
        if (maximumLogChange <= options_.iterationTolerance)
        {
            result.status = PhaseEquilibriumStatus::SUCCESS;
            result.tpd = -std::log(amountSum);
            result.composition = trial;
            if (!std::isfinite(result.tpd))
                result.status = PhaseEquilibriumStatus::NONFINITE_PARAMETER;
            return result;
        }

        if (maximumLogChange
            < bestMaximumLogChange*(1.0 - 1.0e-3))
        {
            bestMaximumLogChange = maximumLogChange;
            stagnantIterations = 0;
        }
        else
        {
            ++stagnantIterations;
        }
        if (stagnantIterations >= 16 && !usedExpandedSeedRestart)
        {
            // An explicit phase root can terminate before a seed reaches its
            // stationary point. Expand that same seed once to cross the root
            // boundary; a second stall still fails closed.
            std::vector<double> restartedLogs = logAmounts;
            bool restartValid = true;
            for (std::size_t index : activeIndices)
            {
                const double logSeed = std::log(initialK[index]);
                restartedLogs[index] = std::log(feed[index])
                    + (vaporLike ? 2.0 : -2.0)*logSeed;
                restartValid = restartValid
                    && std::isfinite(restartedLogs[index])
                    && restartedLogs[index] < maximumLog
                    && restartedLogs[index] > minimumLog;
            }
            if (restartValid)
            {
                logAmounts.swap(restartedLogs);
                bestMaximumLogChange =
                    std::numeric_limits<double>::infinity();
                stagnantIterations = 0;
                usedExpandedSeedRestart = true;
                continue;
            }
            usedExpandedSeedRestart = true;
        }

        const std::size_t activeCount = activeIndices.size();
        // Newton acceleration acts on the stationarity equations in log
        // amounts. Invalid finite-difference roots or a rejected line search
        // fall back to the configured successive-substitution step.
        std::vector<std::vector<double>> jacobian(
            activeCount, std::vector<double>(activeCount, 0.0));
        bool newtonAvailable = true;
        for (std::size_t column = 0;
             column < activeCount && newtonAvailable;
             ++column)
        {
            const std::size_t species = activeIndices[column];
            double differenceStep = 1.0e-6
                *std::max(1.0, std::abs(logAmounts[species]));
            std::vector<double> perturbedLogs = logAmounts;
            perturbedLogs[species] += differenceStep;
            std::vector<double> perturbedTrial;
            std::vector<double> perturbedResidual;
            double perturbedSum = 0.0;
            PhaseEquilibriumStatus perturbedStatus = evaluateResidual(
                perturbedLogs,
                perturbedTrial,
                perturbedResidual,
                perturbedSum);
            if (perturbedStatus != PhaseEquilibriumStatus::SUCCESS)
            {
                differenceStep = -differenceStep;
                perturbedLogs = logAmounts;
                perturbedLogs[species] += differenceStep;
                perturbedStatus = evaluateResidual(
                    perturbedLogs,
                    perturbedTrial,
                    perturbedResidual,
                    perturbedSum);
            }
            if (perturbedStatus != PhaseEquilibriumStatus::SUCCESS)
            {
                newtonAvailable = false;
                break;
            }

            for (std::size_t row = 0; row < activeCount; ++row)
            {
                const std::size_t rowSpecies = activeIndices[row];
                jacobian[row][column] =
                    (perturbedResidual[rowSpecies] - residual[rowSpecies])
                    /differenceStep;
            }
        }

        std::vector<double> newtonRhs(activeCount, 0.0);
        for (std::size_t row = 0; row < activeCount; ++row)
            newtonRhs[row] = -residual[activeIndices[row]];
        std::vector<double> newtonStep;
        newtonAvailable = newtonAvailable
            && solveDenseLinearSystem(jacobian, newtonRhs, newtonStep);

        bool acceptedNewtonStep = false;
        if (newtonAvailable)
        {
            for (double stepScale = 1.0;
                 stepScale >= 1.0/1024.0;
                 stepScale *= 0.5)
            {
                std::vector<double> candidateLogs = logAmounts;
                for (std::size_t column = 0;
                     column < activeCount;
                     ++column)
                {
                    candidateLogs[activeIndices[column]] +=
                        stepScale*newtonStep[column];
                }

                std::vector<double> candidateTrial;
                std::vector<double> candidateResidual;
                double candidateSum = 0.0;
                if (evaluateResidual(
                        candidateLogs,
                        candidateTrial,
                        candidateResidual,
                        candidateSum)
                    != PhaseEquilibriumStatus::SUCCESS)
                {
                    continue;
                }

                double candidateMaximumLogChange = 0.0;
                for (std::size_t index : activeIndices)
                {
                    candidateMaximumLogChange = std::max(
                        candidateMaximumLogChange,
                        std::abs(candidateResidual[index]));
                }
                if (candidateMaximumLogChange
                    <= maximumLogChange*(1.0 - 1.0e-4*stepScale))
                {
                    logAmounts.swap(candidateLogs);
                    acceptedNewtonStep = true;
                    break;
                }
            }
        }
        if (acceptedNewtonStep) continue;

        const double damping = options_.successiveSubstitutionDamping;
        for (std::size_t index : activeIndices)
        {
            logAmounts[index] += damping*residual[index];
            if (!std::isfinite(logAmounts[index])
                || logAmounts[index] >= maximumLog
                || logAmounts[index] <= minimumLog)
            {
                result.status = PhaseEquilibriumStatus::NONFINITE_PARAMETER;
                return result;
            }
        }
    }

    result.status = PhaseEquilibriumStatus::STABILITY_NOT_CONVERGED;
    return result;
}

StabilityResult PRTPFlash::analyzeStability(
    double pressure,
    double temperature,
    const std::vector<double>& feedMoleFractions) const
{
    StabilityResult result;
    result.status = configurationStatus();
    if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;
    result.status = kernel_.validateMoleFractions(feedMoleFractions);
    if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;

    const PRPhaseState feedState = kernel_.phaseState(
        pressure, temperature, feedMoleFractions, PhaseRoot::STABLE);
    if (!feedState.valid())
    {
        result.status = feedState.status;
        return result;
    }

    std::vector<double> initialK;
    result.status = kernel_.wilsonKValues(pressure, temperature, initialK);
    if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;

    // Deterministic multi-start foundation: Wilson K, sqrt(K), and K^2,
    // each in vapor-like and liquid-like directions. This is auditable and
    // stronger than a single Wilson pair, but it is not a proof that every
    // possible TPD stationary point has been found.
    const double seedExponents[3] = {1.0, 0.5, 2.0};
    const std::vector<double> neutralSeed(initialK.size(), 1.0);
    bool haveVaporLike = false;
    bool haveLiquidLike = false;
    result.vaporLikeTpd = std::numeric_limits<double>::infinity();
    result.liquidLikeTpd = std::numeric_limits<double>::infinity();

    for (double exponent : seedExponents)
    {
        std::vector<double> seed(initialK.size(), 0.0);
        bool seedValid = true;
        for (std::size_t i = 0; i < initialK.size(); ++i)
        {
            seed[i] = std::exp(exponent*std::log(initialK[i]));
            seedValid = seedValid && finitePositive(seed[i]);
        }

        for (int direction = 0; direction < 2; ++direction)
        {
            ++result.attemptedTrials;
            if (!seedValid)
            {
                ++result.failedTrials;
                continue;
            }
            const bool vaporLike = direction == 0;
            StabilityTrialResult trial = stabilityTrial(
                pressure,
                temperature,
                feedMoleFractions,
                feedState,
                seed,
                vaporLike);
            if (trial.status
                == PhaseEquilibriumStatus::STABILITY_NOT_CONVERGED)
            {
                // A strongly separated Wilson seed can approach an explicit
                // root boundary without reaching its stationary point even
                // though the feed stationary point is regular. Retry only
                // exhausted iterations from neutral K=1; other failures stay
                // fail-closed, and successful separated seeds still retain
                // every lower TPD stationary point they discovered.
                trial = stabilityTrial(
                    pressure,
                    temperature,
                    feedMoleFractions,
                    feedState,
                    neutralSeed,
                    vaporLike);
            }
            if (trial.status != PhaseEquilibriumStatus::SUCCESS)
            {
                ++result.failedTrials;
                continue;
            }
            ++result.convergedTrials;

            if (vaporLike && trial.tpd < result.vaporLikeTpd)
            {
                haveVaporLike = true;
                result.vaporLikeTpd = trial.tpd;
                result.vaporLikeComposition = trial.composition;
                result.vaporLikeIterations = trial.iterations;
            }
            if (!vaporLike && trial.tpd < result.liquidLikeTpd)
            {
                haveLiquidLike = true;
                result.liquidLikeTpd = trial.tpd;
                result.liquidLikeComposition = trial.composition;
                result.liquidLikeIterations = trial.iterations;
            }
        }
    }

    if (!haveVaporLike || !haveLiquidLike || result.failedTrials != 0)
    {
        if (haveVaporLike && haveLiquidLike)
            result.minimumTpd = std::min(
                result.vaporLikeTpd, result.liquidLikeTpd);
        result.status = PhaseEquilibriumStatus::STABILITY_NOT_CONVERGED;
        return result;
    }

    result.minimumTpd = std::min(
        result.vaporLikeTpd, result.liquidLikeTpd);
    result.unstable = result.minimumTpd < -options_.stabilityTolerance;
    result.status = PhaseEquilibriumStatus::SUCCESS;
    return result;
}

PhaseEquilibriumStatus PRTPFlash::phaseCompositions(
    const std::vector<double>& feed,
    const std::vector<double>& K,
    double vaporMoleFraction,
    std::vector<double>& liquid,
    std::vector<double>& vapor) const
{
    liquid.assign(feed.size(), 0.0);
    vapor.assign(feed.size(), 0.0);
    if (!std::isfinite(vaporMoleFraction)
        || vaporMoleFraction < 0.0 || vaporMoleFraction > 1.0
        || K.size() != feed.size())
        return PhaseEquilibriumStatus::NONFINITE_PARAMETER;

    double liquidSum = 0.0;
    double vaporSum = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i)
    {
        if (!finitePositive(K[i]))
            return PhaseEquilibriumStatus::NONFINITE_PARAMETER;
        const double denominator =
            1.0 + vaporMoleFraction*(K[i] - 1.0);
        if (!finitePositive(denominator))
            return PhaseEquilibriumStatus::NONFINITE_PARAMETER;
        liquid[i] = feed[i]/denominator;
        vapor[i] = K[i]*liquid[i];
        if (feed[i] > 0.0
            && (!finitePositive(liquid[i]) || !finitePositive(vapor[i])))
            return PhaseEquilibriumStatus::NONFINITE_PARAMETER;
        liquidSum += liquid[i];
        vaporSum += vapor[i];
    }
    if (!finitePositive(liquidSum) || !finitePositive(vaporSum))
        return PhaseEquilibriumStatus::NONFINITE_PARAMETER;

    const double permittedRoundoff = 10.0*options_.materialTolerance;
    if (std::abs(liquidSum - 1.0) > permittedRoundoff
        || std::abs(vaporSum - 1.0) > permittedRoundoff)
        return PhaseEquilibriumStatus::RACHFORD_RICE_NOT_CONVERGED;

    // Rachford-Rice guarantees these sums analytically. This explicit
    // roundoff normalization is accepted only after the residual gate above.
    for (std::size_t i = 0; i < feed.size(); ++i)
    {
        liquid[i] /= liquidSum;
        vapor[i] /= vaporSum;
    }
    return PhaseEquilibriumStatus::SUCCESS;
}

PureSaturationResult PRTPFlash::solvePureSaturationTemperature(
    double pressure,
    double minimumTemperature,
    double maximumTemperature) const
{
    PureSaturationResult result;
    result.pressure = pressure;
    result.status = configurationStatus();
    if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;
    if (kernel_.speciesNames().size() != 1)
    {
        result.status = PhaseEquilibriumStatus::INVALID_CONFIGURATION;
        return result;
    }
    if (!finitePositive(pressure))
    {
        result.status = PhaseEquilibriumStatus::INVALID_PRESSURE;
        return result;
    }
    if (!finitePositive(minimumTemperature)
        || !finitePositive(maximumTemperature)
        || !(minimumTemperature < maximumTemperature))
    {
        result.status = PhaseEquilibriumStatus::INVALID_TEMPERATURE;
        return result;
    }
    if (pressure >= kernel_.properties()[0].Pc)
    {
        result.status = PhaseEquilibriumStatus::PURE_SATURATION_NO_BRACKET;
        return result;
    }

    const std::vector<double> pure = {1.0};
    const auto fugacityDifference = [&](
        double temperature,
        double& difference,
        PRPhaseState& liquid,
        PRPhaseState& vapor)
    {
        const PRRootSet rootSet = kernel_.roots(
            pressure, temperature, pure);
        if (!rootSet.valid() || rootSet.physicalRoots.size() < 2) return false;
        liquid = kernel_.phaseState(
            pressure, temperature, pure, PhaseRoot::LIQUID);
        vapor = kernel_.phaseState(
            pressure, temperature, pure, PhaseRoot::VAPOR);
        if (!liquid.valid() || !vapor.valid()
            || liquid.compressibility == vapor.compressibility)
            return false;
        difference = liquid.lnFugacityCoefficients[0]
            - vapor.lnFugacityCoefficients[0];
        return std::isfinite(difference);
    };

    bool havePrevious = false;
    double previousTemperature = 0.0;
    double previousDifference = 0.0;
    double lowerTemperature = 0.0;
    double upperTemperature = 0.0;
    double lowerDifference = 0.0;
    double upperDifference = 0.0;
    bool bracketed = false;

    for (std::size_t interval = 0;
         interval <= options_.saturationScanIntervals;
         ++interval)
    {
        const double fraction = static_cast<double>(interval)
            /static_cast<double>(options_.saturationScanIntervals);
        const double temperature = minimumTemperature
            + fraction*(maximumTemperature - minimumTemperature);
        double difference = 0.0;
        PRPhaseState liquid;
        PRPhaseState vapor;
        if (!fugacityDifference(temperature, difference, liquid, vapor))
        {
            havePrevious = false;
            continue;
        }
        if (std::abs(difference) <= options_.fugacityTolerance)
        {
            result.status = PhaseEquilibriumStatus::SUCCESS;
            result.temperature = temperature;
            result.liquidCompressibility = liquid.compressibility;
            result.vaporCompressibility = vapor.compressibility;
            result.logFugacityResidual = difference;
            return result;
        }
        if (havePrevious && previousDifference*difference < 0.0)
        {
            lowerTemperature = previousTemperature;
            lowerDifference = previousDifference;
            upperTemperature = temperature;
            upperDifference = difference;
            bracketed = true;
            break;
        }
        havePrevious = true;
        previousTemperature = temperature;
        previousDifference = difference;
    }

    if (!bracketed)
    {
        result.status = PhaseEquilibriumStatus::PURE_SATURATION_NO_BRACKET;
        return result;
    }

    for (std::size_t iteration = 1;
         iteration <= options_.maxSaturationIterations;
         ++iteration)
    {
        const double middleTemperature =
            0.5*(lowerTemperature + upperTemperature);
        double middleDifference = 0.0;
        PRPhaseState liquid;
        PRPhaseState vapor;
        if (!fugacityDifference(
                middleTemperature, middleDifference, liquid, vapor))
        {
            result.status =
                PhaseEquilibriumStatus::PURE_SATURATION_NOT_CONVERGED;
            return result;
        }

        result.iterations = iteration;
        result.temperature = middleTemperature;
        result.liquidCompressibility = liquid.compressibility;
        result.vaporCompressibility = vapor.compressibility;
        result.logFugacityResidual = middleDifference;

        const double temperatureTolerance =
            options_.saturationTemperatureTolerance
            *std::max(1.0, std::abs(middleTemperature));
        if (std::abs(middleDifference) <= options_.fugacityTolerance
            && upperTemperature - lowerTemperature <= temperatureTolerance)
        {
            result.status = PhaseEquilibriumStatus::SUCCESS;
            return result;
        }
        if (lowerDifference*middleDifference <= 0.0)
        {
            upperTemperature = middleTemperature;
            upperDifference = middleDifference;
        }
        else
        {
            lowerTemperature = middleTemperature;
            lowerDifference = middleDifference;
        }
        (void)upperDifference;
    }

    result.status = PhaseEquilibriumStatus::PURE_SATURATION_NOT_CONVERGED;
    return result;
}

TPFlashResult PRTPFlash::flashPure(
    double pressure,
    double temperature,
    const std::vector<double>& feed) const
{
    TPFlashResult result;
    result.pressure = pressure;
    result.temperature = temperature;
    result.feedMoleFractions = feed;
    result.status = kernel_.moleToMassFractions(feed, result.feedMassFractions);
    if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;

    const CriticalProperties& properties = kernel_.properties()[0];
    const PRRootSet rootSet = kernel_.roots(pressure, temperature, feed);
    if (!rootSet.valid())
    {
        result.status = rootSet.status;
        return result;
    }

    const PRPhaseState liquid = kernel_.phaseState(
        pressure, temperature, feed, PhaseRoot::LIQUID);
    const PRPhaseState vapor = kernel_.phaseState(
        pressure, temperature, feed, PhaseRoot::VAPOR);
    if (!liquid.valid() || !vapor.valid())
    {
        result.status = !liquid.valid() ? liquid.status : vapor.status;
        return result;
    }

    if (rootSet.physicalRoots.size() >= 2)
    {
        const double logResidual = liquid.lnFugacityCoefficients[0]
            - vapor.lnFugacityCoefficients[0];
        if (std::abs(logResidual) <= options_.fugacityTolerance)
        {
            result.status = PhaseEquilibriumStatus::SUCCESS;
            result.region = EquilibriumRegion::PURE_COMPONENT_COEXISTENCE;
            result.phaseFractionDetermined = false;
            result.vaporMoleFraction =
                std::numeric_limits<double>::quiet_NaN();
            result.liquidMoleFractions = feed;
            result.vaporMoleFractions = feed;
            result.liquidMassFractions = result.feedMassFractions;
            result.vaporMassFractions = result.feedMassFractions;
            result.liquidState = liquid;
            result.vaporState = vapor;
            result.maximumLogFugacityResidual = std::abs(logResidual);
            result.maximumMoleMaterialResidual = 0.0;
            return result;
        }
    }

    if (temperature >= properties.Tc && pressure >= properties.Pc)
    {
        const PRPhaseState stable = kernel_.phaseState(
            pressure, temperature, feed, PhaseRoot::STABLE);
        if (!stable.valid())
        {
            result.status = stable.status;
            return result;
        }
        result.status = PhaseEquilibriumStatus::SUCCESS;
        result.region = EquilibriumRegion::SINGLE_SUPERCRITICAL;
        result.phaseFractionDetermined = false;
        result.vaporMoleFraction =
            std::numeric_limits<double>::quiet_NaN();
        result.singlePhaseMoleFractions = feed;
        result.singlePhaseMassFractions = result.feedMassFractions;
        result.singlePhaseState = stable;
        result.maximumMoleMaterialResidual = 0.0;
        return result;
    }

    bool liquidRegion = false;
    if (pressure >= properties.Pc)
    {
        liquidRegion = temperature < properties.Tc;
    }
    else if (temperature >= properties.Tc)
    {
        liquidRegion = false;
    }
    else
    {
        const double minimumSaturationTemperature =
            std::max(1.0, 0.2*properties.Tc);
        const double maximumSaturationTemperature =
            properties.Tc*(1.0 - 1.0e-10);
        const PureSaturationResult saturation =
            solvePureSaturationTemperature(
                pressure,
                minimumSaturationTemperature,
                maximumSaturationTemperature);
        if (!saturation.valid())
        {
            result.status = saturation.status;
            return result;
        }
        const double temperatureTolerance =
            options_.saturationTemperatureTolerance
            *std::max(1.0, saturation.temperature);
        if (std::abs(temperature - saturation.temperature)
            <= temperatureTolerance)
        {
            const PRPhaseState saturationLiquid = kernel_.phaseState(
                pressure, saturation.temperature, feed, PhaseRoot::LIQUID);
            const PRPhaseState saturationVapor = kernel_.phaseState(
                pressure, saturation.temperature, feed, PhaseRoot::VAPOR);
            if (!saturationLiquid.valid() || !saturationVapor.valid())
            {
                result.status = PhaseEquilibriumStatus::NO_PHYSICAL_ROOT;
                return result;
            }
            result.status = PhaseEquilibriumStatus::SUCCESS;
            result.region = EquilibriumRegion::PURE_COMPONENT_COEXISTENCE;
            result.phaseFractionDetermined = false;
            result.vaporMoleFraction =
                std::numeric_limits<double>::quiet_NaN();
            result.liquidMoleFractions = feed;
            result.vaporMoleFractions = feed;
            result.liquidMassFractions = result.feedMassFractions;
            result.vaporMassFractions = result.feedMassFractions;
            result.liquidState = saturationLiquid;
            result.vaporState = saturationVapor;
            result.maximumLogFugacityResidual = std::abs(
                saturationLiquid.lnFugacityCoefficients[0]
                - saturationVapor.lnFugacityCoefficients[0]);
            result.maximumMoleMaterialResidual = 0.0;
            return result;
        }
        liquidRegion = temperature < saturation.temperature;
    }

    result.status = PhaseEquilibriumStatus::SUCCESS;
    result.phaseFractionDetermined = true;
    result.maximumMoleMaterialResidual = 0.0;
    if (liquidRegion)
    {
        result.region = EquilibriumRegion::SINGLE_LIQUID;
        result.vaporMoleFraction = 0.0;
        result.liquidMoleFractions = feed;
        result.liquidMassFractions = result.feedMassFractions;
        result.liquidState = liquid;
        result.singlePhaseMoleFractions = feed;
        result.singlePhaseMassFractions = result.feedMassFractions;
        result.singlePhaseState = liquid;
    }
    else
    {
        result.region = EquilibriumRegion::SINGLE_VAPOR;
        result.vaporMoleFraction = 1.0;
        result.vaporMoleFractions = feed;
        result.vaporMassFractions = result.feedMassFractions;
        result.vaporState = vapor;
        result.singlePhaseMoleFractions = feed;
        result.singlePhaseMassFractions = result.feedMassFractions;
        result.singlePhaseState = vapor;
    }
    return result;
}

TPFlashResult PRTPFlash::flashMoleFractions(
    double pressure,
    double temperature,
    const std::vector<double>& feedMoleFractions) const
{
    TPFlashResult result;
    result.pressure = pressure;
    result.temperature = temperature;
    result.feedMoleFractions = feedMoleFractions;
    result.status = configurationStatus();
    if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;
    if (!finitePositive(pressure))
    {
        result.status = PhaseEquilibriumStatus::INVALID_PRESSURE;
        return result;
    }
    if (!finitePositive(temperature))
    {
        result.status = PhaseEquilibriumStatus::INVALID_TEMPERATURE;
        return result;
    }
    result.status = kernel_.validateMoleFractions(feedMoleFractions);
    if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;
    if (feedMoleFractions.size() == 1)
        return flashPure(pressure, temperature, feedMoleFractions);

    result.status = kernel_.moleToMassFractions(
        feedMoleFractions, result.feedMassFractions);
    if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;

    const StabilityResult feedStability = analyzeStability(
        pressure, temperature, feedMoleFractions);
    if (!feedStability.valid())
    {
        result.status = feedStability.status;
        return result;
    }
    result.feedMinimumTpd = feedStability.minimumTpd;

    std::vector<double> K;
    result.status = kernel_.wilsonKValues(pressure, temperature, K);
    if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;
    std::vector<std::size_t> flashActiveIndices;
    flashActiveIndices.reserve(feedMoleFractions.size());
    for (std::size_t i = 0; i < feedMoleFractions.size(); ++i)
    {
        if (feedMoleFractions[i] > 0.0) flashActiveIndices.push_back(i);
    }

    if (!feedStability.unstable)
    {
        const PRPhaseState stableState = kernel_.phaseState(
            pressure, temperature, feedMoleFractions, PhaseRoot::STABLE);
        const PRRootSet rootSet = kernel_.roots(
            pressure, temperature, feedMoleFractions);
        if (!stableState.valid() || !rootSet.valid())
        {
            result.status = !stableState.valid()
                ? stableState.status : rootSet.status;
            return result;
        }

        double rrAtZero = 0.0;
        double rrAtOne = 0.0;
        for (std::size_t i = 0; i < K.size(); ++i)
        {
            rrAtZero += feedMoleFractions[i]*(K[i] - 1.0);
            rrAtOne += feedMoleFractions[i]*(1.0 - 1.0/K[i]);
        }

        result.status = PhaseEquilibriumStatus::SUCCESS;
        result.maximumMoleMaterialResidual = 0.0;
        if (rootSet.physicalRoots.size() > 1
            && stableState.compressibility
                == rootSet.physicalRoots.front())
        {
            result.region = EquilibriumRegion::SINGLE_LIQUID;
        }
        else if (rootSet.physicalRoots.size() > 1
            && stableState.compressibility
                == rootSet.physicalRoots.back())
        {
            result.region = EquilibriumRegion::SINGLE_VAPOR;
        }
        else if (rrAtZero < 0.0 && rrAtOne < 0.0)
        {
            result.region = EquilibriumRegion::SINGLE_LIQUID;
        }
        else if (rrAtZero > 0.0 && rrAtOne > 0.0)
        {
            result.region = EquilibriumRegion::SINGLE_VAPOR;
        }
        else
        {
            result.region = EquilibriumRegion::SINGLE_SUPERCRITICAL;
        }

        if (result.region == EquilibriumRegion::SINGLE_LIQUID)
        {
            result.phaseFractionDetermined = true;
            result.vaporMoleFraction = 0.0;
            result.liquidMoleFractions = feedMoleFractions;
            result.liquidMassFractions = result.feedMassFractions;
            result.liquidState = kernel_.phaseState(
                pressure, temperature, feedMoleFractions, PhaseRoot::LIQUID);
            result.singlePhaseMoleFractions = feedMoleFractions;
            result.singlePhaseMassFractions = result.feedMassFractions;
            result.singlePhaseState = result.liquidState;
        }
        else if (result.region == EquilibriumRegion::SINGLE_VAPOR)
        {
            result.phaseFractionDetermined = true;
            result.vaporMoleFraction = 1.0;
            result.vaporMoleFractions = feedMoleFractions;
            result.vaporMassFractions = result.feedMassFractions;
            result.vaporState = kernel_.phaseState(
                pressure, temperature, feedMoleFractions, PhaseRoot::VAPOR);
            result.singlePhaseMoleFractions = feedMoleFractions;
            result.singlePhaseMassFractions = result.feedMassFractions;
            result.singlePhaseState = result.vaporState;
        }
        else
        {
            result.singlePhaseMoleFractions = feedMoleFractions;
            result.singlePhaseMassFractions = result.feedMassFractions;
            result.singlePhaseState = stableState;
        }
        return result;
    }

    std::vector<std::vector<double>> candidates;
    const auto appendCandidateFamily = [&]
    (
        const std::vector<double>& seed
    )
    {
        candidates.push_back(seed);
        if (flashActiveIndices.size() != 2) return;

        const double seedExponents[2] = {0.5, 2.0};
        for (double exponent : seedExponents)
        {
            std::vector<double> candidate(seed.size(), 1.0);
            bool valid = true;
            for (std::size_t i = 0; i < seed.size(); ++i)
            {
                candidate[i] = std::exp(exponent*std::log(seed[i]));
                valid = valid && finitePositive(candidate[i]);
            }
            if (valid) candidates.push_back(candidate);
        }
    };
    appendCandidateFamily(K);
    if (feedStability.vaporLikeTpd < -options_.stabilityTolerance)
    {
        std::vector<double> candidate(K.size(), 1.0);
        bool valid = true;
        for (std::size_t i = 0; i < K.size(); ++i)
        {
            if (feedMoleFractions[i] == 0.0) continue;
            candidate[i] = feedStability.vaporLikeComposition[i]
                /feedMoleFractions[i];
            valid = valid && finitePositive(candidate[i]);
        }
        if (valid) appendCandidateFamily(candidate);
    }
    if (feedStability.liquidLikeTpd < -options_.stabilityTolerance)
    {
        std::vector<double> candidate(K.size(), 1.0);
        bool valid = true;
        for (std::size_t i = 0; i < K.size(); ++i)
        {
            if (feedMoleFractions[i] == 0.0) continue;
            candidate[i] = feedMoleFractions[i]
                /feedStability.liquidLikeComposition[i];
            valid = valid && finitePositive(candidate[i]);
        }
        if (valid) appendCandidateFamily(candidate);
    }
    if (feedStability.vaporLikeTpd < -options_.stabilityTolerance
        && feedStability.liquidLikeTpd < -options_.stabilityTolerance)
    {
        std::vector<double> candidate(K.size(), 1.0);
        bool valid = true;
        for (std::size_t i = 0; i < K.size(); ++i)
        {
            if (feedMoleFractions[i] == 0.0) continue;
            candidate[i] = feedStability.vaporLikeComposition[i]
                /feedStability.liquidLikeComposition[i];
            valid = valid && finitePositive(candidate[i]);
        }
        if (valid) appendCandidateFamily(candidate);
    }

    const TPFlashResult baseResult = result;
    const auto solveCandidate = [&]
    (
        const std::vector<double>& initialK
    )
    {
        TPFlashResult result = baseResult;
        std::vector<double> K = initialK;
        RachfordRiceResult rachfordRice =
            solveRachfordRice(feedMoleFractions, K);
        if (!rachfordRice.valid())
        {
            result.status = rachfordRice.status;
            return result;
        }

        std::vector<double> logK(K.size(), 0.0);
        for (std::size_t i = 0; i < K.size(); ++i)
            logK[i] = std::log(K[i]);

    const auto evaluateFlashResidual = [&]
    (
        const std::vector<double>& candidateLogK,
        std::vector<double>& correction
    )
    {
        if (candidateLogK.size() != K.size()) return false;
        std::vector<double> candidateK(candidateLogK.size(), 0.0);
        for (std::size_t i = 0; i < candidateLogK.size(); ++i)
        {
            candidateK[i] = std::exp(candidateLogK[i]);
            if (!finitePositive(candidateK[i])) return false;
        }

        const RachfordRiceResult candidateRachfordRice =
            solveRachfordRice(feedMoleFractions, candidateK);
        if (!candidateRachfordRice.valid()) return false;
        std::vector<double> candidateLiquid;
        std::vector<double> candidateVapor;
        if (phaseCompositions(
                feedMoleFractions,
                candidateK,
                candidateRachfordRice.vaporMoleFraction,
                candidateLiquid,
                candidateVapor)
            != PhaseEquilibriumStatus::SUCCESS)
        {
            return false;
        }
        const PRPhaseState candidateLiquidState = kernel_.phaseState(
            pressure, temperature, candidateLiquid, PhaseRoot::LIQUID);
        const PRPhaseState candidateVaporState = kernel_.phaseState(
            pressure, temperature, candidateVapor, PhaseRoot::VAPOR);
        if (!candidateLiquidState.valid() || !candidateVaporState.valid())
            return false;

        correction.assign(candidateLogK.size(), 0.0);
        for (std::size_t i = 0; i < candidateLogK.size(); ++i)
        {
            correction[i] = candidateLiquidState.lnFugacityCoefficients[i]
                - candidateVaporState.lnFugacityCoefficients[i]
                - candidateLogK[i];
            if (!std::isfinite(correction[i])) return false;
        }
        return true;
    };

    for (std::size_t iteration = 1;
         iteration <= options_.maxFlashIterations;
         ++iteration)
    {
        rachfordRice = solveRachfordRice(feedMoleFractions, K);
        if (!rachfordRice.valid())
        {
            result.status = rachfordRice.status;
            return result;
        }

        std::vector<double> liquid;
        std::vector<double> vapor;
        result.status = phaseCompositions(
            feedMoleFractions,
            K,
            rachfordRice.vaporMoleFraction,
            liquid,
            vapor);
        if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;

        const PRPhaseState liquidState = kernel_.phaseState(
            pressure, temperature, liquid, PhaseRoot::LIQUID);
        const PRPhaseState vaporState = kernel_.phaseState(
            pressure, temperature, vapor, PhaseRoot::VAPOR);
        if (!liquidState.valid() || !vaporState.valid())
        {
            result.status = !liquidState.valid()
                ? liquidState.status : vaporState.status;
            return result;
        }

        double maximumLogFugacityResidual = 0.0;
        double maximumLogKCorrection = 0.0;
        std::vector<double> targetLogK(logK.size(), 0.0);
        for (std::size_t i = 0; i < K.size(); ++i)
        {
            targetLogK[i] = liquidState.lnFugacityCoefficients[i]
                - vaporState.lnFugacityCoefficients[i];
            const double correction = targetLogK[i] - logK[i];
            if (feedMoleFractions[i] > 0.0
                || flashActiveIndices.size() != 2)
            {
                maximumLogKCorrection = std::max(
                    maximumLogKCorrection, std::abs(correction));
            }
            if (feedMoleFractions[i] > 0.0)
            {
                const double fugacityResidual = std::log(liquid[i])
                    + liquidState.lnFugacityCoefficients[i]
                    - std::log(vapor[i])
                    - vaporState.lnFugacityCoefficients[i];
                maximumLogFugacityResidual = std::max(
                    maximumLogFugacityResidual,
                    std::abs(fugacityResidual));
            }
        }
        const double materialResidual = maximumAbsoluteMaterialResidual(
            feedMoleFractions,
            liquid,
            vapor,
            rachfordRice.vaporMoleFraction);

        result.flashIterations = iteration;
        result.rachfordRiceResidual = std::abs(rachfordRice.residual);
        result.maximumLogFugacityResidual = maximumLogFugacityResidual;
        result.maximumMoleMaterialResidual = materialResidual;

        if (maximumLogFugacityResidual <= options_.fugacityTolerance
            && maximumLogKCorrection <= options_.iterationTolerance
            && materialResidual <= options_.materialTolerance
            && std::abs(rachfordRice.residual) <= options_.materialTolerance)
        {
            const StabilityResult liquidStability = analyzeStability(
                pressure, temperature, liquid);
            const StabilityResult vaporStability = analyzeStability(
                pressure, temperature, vapor);
            if (!liquidStability.valid() || !vaporStability.valid())
            {
                result.status = PhaseEquilibriumStatus::STABILITY_NOT_CONVERGED;
                return result;
            }
            result.liquidMinimumTpd = liquidStability.minimumTpd;
            result.vaporMinimumTpd = vaporStability.minimumTpd;
            if (liquidStability.minimumTpd
                    < -5.0*options_.stabilityTolerance
                || vaporStability.minimumTpd
                    < -5.0*options_.stabilityTolerance)
            {
                result.status = PhaseEquilibriumStatus::FLASH_NOT_CONVERGED;
                return result;
            }

            result.status = kernel_.moleToMassFractions(
                liquid, result.liquidMassFractions);
            if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;
            result.status = kernel_.moleToMassFractions(
                vapor, result.vaporMassFractions);
            if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;
            result.status = PhaseEquilibriumStatus::SUCCESS;
            result.region = EquilibriumRegion::TWO_PHASE;
            result.phaseFractionDetermined = true;
            result.vaporMoleFraction = rachfordRice.vaporMoleFraction;
            result.liquidMoleFractions = liquid;
            result.vaporMoleFractions = vapor;
            result.liquidState = liquidState;
            result.vaporState = vaporState;
            return result;
        }

        if (flashActiveIndices.size() != 2)
        {
            // The safeguarded log-K Newton path below is currently closed by
            // the binary topology oracle. Preserve the established generic
            // multicomponent substitution path until it has its own oracle.
            const double damping = options_.successiveSubstitutionDamping;
            for (std::size_t i = 0; i < K.size(); ++i)
            {
                logK[i] += damping*(targetLogK[i] - logK[i]);
                K[i] = std::exp(logK[i]);
                if (!finitePositive(K[i]))
                {
                    result.status =
                        PhaseEquilibriumStatus::NONFINITE_PARAMETER;
                    return result;
                }
            }
            continue;
        }

        const std::size_t componentCount = flashActiveIndices.size();
        std::vector<std::vector<double>> jacobian(
            componentCount, std::vector<double>(componentCount, 0.0));
        bool newtonAvailable = true;
        for (std::size_t column = 0;
             column < componentCount && newtonAvailable;
             ++column)
        {
            const std::size_t species = flashActiveIndices[column];
            double differenceStep = 1.0e-6
                *std::max(1.0, std::abs(logK[species]));
            std::vector<double> perturbedLogK = logK;
            perturbedLogK[species] += differenceStep;
            std::vector<double> perturbedCorrection;
            if (!evaluateFlashResidual(
                    perturbedLogK, perturbedCorrection))
            {
                differenceStep = -differenceStep;
                perturbedLogK = logK;
                perturbedLogK[species] += differenceStep;
                if (!evaluateFlashResidual(
                        perturbedLogK, perturbedCorrection))
                {
                    newtonAvailable = false;
                    break;
                }
            }
            for (std::size_t row = 0; row < componentCount; ++row)
            {
                const std::size_t rowSpecies = flashActiveIndices[row];
                const double baseCorrection =
                    targetLogK[rowSpecies] - logK[rowSpecies];
                jacobian[row][column] =
                    (perturbedCorrection[rowSpecies] - baseCorrection)
                    /differenceStep;
            }
        }

        std::vector<double> newtonRhs(componentCount, 0.0);
        for (std::size_t row = 0; row < componentCount; ++row)
        {
            const std::size_t species = flashActiveIndices[row];
            newtonRhs[row] = logK[species] - targetLogK[species];
        }
        std::vector<double> newtonStep;
        newtonAvailable = newtonAvailable
            && solveDenseLinearSystem(jacobian, newtonRhs, newtonStep);

        bool acceptedNewtonStep = false;
        if (newtonAvailable)
        {
            for (double stepScale = 1.0;
                 stepScale >= 1.0/1024.0;
                 stepScale *= 0.5)
            {
                std::vector<double> candidateLogK = logK;
                for (std::size_t component = 0;
                     component < componentCount;
                     ++component)
                {
                    candidateLogK[flashActiveIndices[component]] +=
                        stepScale*newtonStep[component];
                }
                std::vector<double> candidateCorrection;
                if (!evaluateFlashResidual(
                        candidateLogK, candidateCorrection))
                {
                    continue;
                }
                double candidateMaximumCorrection = 0.0;
                for (std::size_t species : flashActiveIndices)
                {
                    candidateMaximumCorrection = std::max(
                        candidateMaximumCorrection,
                        std::abs(candidateCorrection[species]));
                }
                if (candidateMaximumCorrection
                    <= maximumLogKCorrection*(1.0 - 1.0e-4*stepScale))
                {
                    logK.swap(candidateLogK);
                    for (std::size_t species : flashActiveIndices)
                    {
                        K[species] = std::exp(logK[species]);
                    }
                    acceptedNewtonStep = true;
                    break;
                }
            }
        }
        if (acceptedNewtonStep) continue;

        bool acceptedSubstitutionStep = false;
        for (double stepScale = options_.successiveSubstitutionDamping;
             stepScale >= options_.successiveSubstitutionDamping/1024.0;
             stepScale *= 0.5)
        {
            std::vector<double> candidateLogK = logK;
            for (std::size_t species : flashActiveIndices)
            {
                candidateLogK[species] += stepScale
                    *(targetLogK[species] - logK[species]);
            }
            std::vector<double> candidateCorrection;
            if (!evaluateFlashResidual(
                    candidateLogK, candidateCorrection))
            {
                continue;
            }
            double candidateMaximumCorrection = 0.0;
            for (std::size_t species : flashActiveIndices)
            {
                candidateMaximumCorrection = std::max(
                    candidateMaximumCorrection,
                    std::abs(candidateCorrection[species]));
            }
            if (candidateMaximumCorrection
                <= maximumLogKCorrection*(1.0 - 1.0e-4*stepScale))
            {
                logK.swap(candidateLogK);
                for (std::size_t species : flashActiveIndices)
                    K[species] = std::exp(logK[species]);
                acceptedSubstitutionStep = true;
                break;
            }
        }
        if (!acceptedSubstitutionStep)
        {
            result.status = PhaseEquilibriumStatus::FLASH_NOT_CONVERGED;
            return result;
        }
    }

        result.status = PhaseEquilibriumStatus::FLASH_NOT_CONVERGED;
        return result;
    };

    if (flashActiveIndices.size() != 2)
    {
        // Preserve the established generic path: select the first seed with
        // an interior RR bracket, then run one complete substitution solve.
        TPFlashResult noBracket = baseResult;
        for (const std::vector<double>& candidate : candidates)
        {
            const RachfordRiceResult candidateRachfordRice =
                solveRachfordRice(feedMoleFractions, candidate);
            if (candidateRachfordRice.valid())
                return solveCandidate(candidate);
            noBracket.status = candidateRachfordRice.status;
        }
        return noBracket;
    }

    // Binary topology is closed by the dedicated global-tie oracle. Run the
    // complete safeguarded solve from Wilson first, then each independent
    // stability-derived seed. A seed is accepted only after fugacity,
    // material, RR, and post-split stability checks all pass.
    TPFlashResult bestFailure = baseResult;
    bool haveFailure = false;
    for (const std::vector<double>& candidate : candidates)
    {
        TPFlashResult attempt = solveCandidate(candidate);
        if (attempt.valid()) return attempt;

        const bool attemptFinite =
            std::isfinite(attempt.maximumLogFugacityResidual);
        const bool bestFinite =
            std::isfinite(bestFailure.maximumLogFugacityResidual);
        const bool betterFailure = !haveFailure
            || (attemptFinite && !bestFinite)
            || (attemptFinite && bestFinite
                && attempt.maximumLogFugacityResidual
                    < bestFailure.maximumLogFugacityResidual)
            || (attemptFinite == bestFinite
                && attempt.flashIterations > bestFailure.flashIterations);
        if (betterFailure)
        {
            bestFailure = std::move(attempt);
            haveFailure = true;
        }
    }
    if (!haveFailure)
        bestFailure.status = PhaseEquilibriumStatus::FLASH_NOT_CONVERGED;
    return bestFailure;
}

TPFlashResult PRTPFlash::flashMassFractions(
    double pressure,
    double temperature,
    const std::vector<double>& feedMassFractions) const
{
    std::vector<double> feedMoleFractions;
    const PhaseEquilibriumStatus conversionStatus =
        kernel_.massToMoleFractions(feedMassFractions, feedMoleFractions);
    if (conversionStatus != PhaseEquilibriumStatus::SUCCESS)
    {
        TPFlashResult result;
        result.status = conversionStatus;
        result.pressure = pressure;
        result.temperature = temperature;
        result.feedMassFractions = feedMassFractions;
        return result;
    }

    TPFlashResult result = flashMoleFractions(
        pressure, temperature, feedMoleFractions);
    if (result.valid()) result.feedMassFractions = feedMassFractions;
    return result;
}

} // End namespace phaseEquilibrium
} // End namespace Foam
