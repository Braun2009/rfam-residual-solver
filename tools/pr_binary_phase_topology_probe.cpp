#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "PhaseEquilibrium/PRMixtureKernel.H"
#include "PhaseEquilibrium/PRTPFlash.H"
#include "MixingRules/BinaryInteractionDB.H"
#include "MixingRules/SpeciesCriticalDB.H"

using namespace Foam;
using namespace Foam::phaseEquilibrium;

namespace
{

constexpr const char* schema = "rfam.pr-binary-phase-topology-screen/v1";
constexpr double slopeTolerance = 1.0e-11;
constexpr double contactTolerance = 1.0e-8;
constexpr double tpEndpointTolerance = 1.0e-7;
constexpr double sharedPhaseTolerance = 1.0e-7;
constexpr double commonTangentTolerance = 1.0e-8;
constexpr std::size_t minimumSkippedSamples = 3;

enum class StableRootKind
{
    SINGLE,
    LIQUID,
    VAPOR,
    INTERMEDIATE,
    INVALID
};

enum class TPTieDisposition
{
    MATCH,
    FLASH_FAILURE,
    NON_TWO_PHASE,
    DIFFERENT_SPLIT
};

struct Sample
{
    double oxygenMoleFraction = std::numeric_limits<double>::quiet_NaN();
    double reducedGibbs = std::numeric_limits<double>::quiet_NaN();
    double compressibility = std::numeric_limits<double>::quiet_NaN();
    std::array<double, 2> reducedChemicalPotential = {
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};
    StableRootKind rootKind = StableRootKind::INVALID;
    bool valid = false;
};

struct TopologyRow
{
    std::size_t validSamples = 0;
    std::size_t invalidSamples = 0;
    std::size_t hullVertices = 0;
    std::size_t tieSegments = 0;
    std::size_t lleSegments = 0;
    std::size_t vleSegments = 0;
    std::size_t otherSegments = 0;
    std::size_t maximumContactClusters = 0;
    std::size_t threePhaseTangentPairs = 0;
    std::size_t threePhaseTieSegments = 0;
    double maximumThreePhaseSharedPhaseError =
        std::numeric_limits<double>::quiet_NaN();
    double maximumThreePhaseSlopeError =
        std::numeric_limits<double>::quiet_NaN();
    double maximumThreePhaseInterceptError =
        std::numeric_limits<double>::quiet_NaN();
    std::size_t maximumSkippedSamples = 0;
    std::size_t refinedTieSegments = 0;
    std::size_t unrefinedTieSegments = 0;
    std::size_t tpSupportedTieSegments = 0;
    std::size_t tpDisagreements = 0;
    std::size_t tpFlashFailures = 0;
    std::size_t tpNonTwoPhase = 0;
    std::size_t tpDifferentSplits = 0;
    std::size_t tpStabilityAttempted = 0;
    std::size_t tpStabilityConverged = 0;
    std::size_t tpStabilityFailed = 0;
    std::size_t tpFlashIterations = 0;
    double tpMaximumFugacityResidual =
        std::numeric_limits<double>::quiet_NaN();
    double tpMaximumMaterialResidual =
        std::numeric_limits<double>::quiet_NaN();
    double tpLiquidMinimumTpd = std::numeric_limits<double>::quiet_NaN();
    double tpVaporMinimumTpd = std::numeric_limits<double>::quiet_NaN();
    double widestTieSpan = 0.0;
    double widestTieLeftOxygen = std::numeric_limits<double>::quiet_NaN();
    double widestTieRightOxygen = std::numeric_limits<double>::quiet_NaN();
    double widestTieLeftZ = std::numeric_limits<double>::quiet_NaN();
    double widestTieRightZ = std::numeric_limits<double>::quiet_NaN();
    double widestTieFugacityResidual =
        std::numeric_limits<double>::infinity();
    std::string tieLeftStabilityStatus = "NONE";
    std::size_t tieLeftStabilityAttempted = 0;
    std::size_t tieLeftStabilityConverged = 0;
    std::size_t tieLeftStabilityFailed = 0;
    double tieLeftMinimumTpd = std::numeric_limits<double>::quiet_NaN();
    std::string tieRightStabilityStatus = "NONE";
    std::size_t tieRightStabilityAttempted = 0;
    std::size_t tieRightStabilityConverged = 0;
    std::size_t tieRightStabilityFailed = 0;
    double tieRightMinimumTpd = std::numeric_limits<double>::quiet_NaN();
    double maximumTpEndpointError = 0.0;
    std::string tpFailureStatus = "NONE";
    bool possibleThreePhase = false;
};

struct RefinedTie
{
    bool success = false;
    Sample left;
    Sample right;
    double fugacityResidual = std::numeric_limits<double>::infinity();
    std::size_t iterations = 0;
};

struct TieEvidence
{
    Sample left;
    Sample right;
    std::size_t contactClusters = 0;
    bool refined = false;
};

struct SingleState
{
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double kij = std::numeric_limits<double>::quiet_NaN();
};

bool close(double lhs, double rhs)
{
    return std::abs(lhs - rhs)
        <= 1.0e-10*std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

double idealMixingGibbs(const std::vector<double>& moleFractions)
{
    double result = 0.0;
    for (double fraction : moleFractions)
    {
        if (fraction > 0.0) result += fraction*std::log(fraction);
    }
    return result;
}

StableRootKind classifyRoot(
    const PRRootSet& roots,
    const PRPhaseState& stable)
{
    if (!roots.valid() || !stable.valid() || roots.physicalRoots.empty())
        return StableRootKind::INVALID;
    if (roots.physicalRoots.size() == 1) return StableRootKind::SINGLE;
    if (close(stable.compressibility, roots.physicalRoots.front()))
        return StableRootKind::LIQUID;
    if (close(stable.compressibility, roots.physicalRoots.back()))
        return StableRootKind::VAPOR;
    return StableRootKind::INTERMEDIATE;
}

Sample evaluateSample(
    const PRMixtureKernel& kernel,
    double pressure,
    double temperature,
    double oxygenMoleFraction)
{
    Sample result;
    result.oxygenMoleFraction = oxygenMoleFraction;
    const std::vector<double> composition = {
        oxygenMoleFraction, 1.0 - oxygenMoleFraction};
    const PRRootSet roots = kernel.roots(
        pressure, temperature, composition);
    const PRPhaseState stable = kernel.phaseState(
        pressure, temperature, composition, PhaseRoot::STABLE);
    if (!roots.valid() || !stable.valid()) return result;

    result.reducedGibbs = idealMixingGibbs(composition)
        + stable.residualGibbsOverRT;
    result.compressibility = stable.compressibility;
    for (std::size_t species = 0; species < composition.size(); ++species)
    {
        if (composition[species] > 0.0)
        {
            result.reducedChemicalPotential[species] =
                std::log(composition[species])
                + stable.lnFugacityCoefficients[species];
        }
    }
    result.rootKind = classifyRoot(roots, stable);
    result.valid = std::isfinite(result.reducedGibbs)
        && std::isfinite(result.compressibility)
        && result.rootKind != StableRootKind::INVALID;
    return result;
}

double slope(const Sample& lhs, const Sample& rhs)
{
    return (rhs.reducedGibbs - lhs.reducedGibbs)
        /(rhs.oxygenMoleFraction - lhs.oxygenMoleFraction);
}

double intercept(const Sample& lhs, const Sample& rhs)
{
    const double tieSlope = slope(lhs, rhs);
    return lhs.reducedGibbs - tieSlope*lhs.oxygenMoleFraction;
}

void updateMaximum(double& current, double value)
{
    if (!std::isfinite(current)) current = value;
    else current = std::max(current, value);
}

double fugacityResidual(const Sample& lhs, const Sample& rhs)
{
    double residual = 0.0;
    for (std::size_t species = 0; species < 2; ++species)
    {
        if (!std::isfinite(lhs.reducedChemicalPotential[species])
            || !std::isfinite(rhs.reducedChemicalPotential[species]))
            return std::numeric_limits<double>::infinity();
        residual = std::max(
            residual,
            std::abs(lhs.reducedChemicalPotential[species]
                - rhs.reducedChemicalPotential[species]));
    }
    return residual;
}

bool chemicalPotentialDerivative(
    const PRMixtureKernel& kernel,
    double pressure,
    double temperature,
    double oxygenMoleFraction,
    std::array<double, 2>& derivative)
{
    const double distance = std::min(
        oxygenMoleFraction, 1.0 - oxygenMoleFraction);
    const double step = std::min(1.0e-5, 0.25*distance);
    if (!std::isfinite(step) || step <= 1.0e-10) return false;
    const Sample lower = evaluateSample(
        kernel, pressure, temperature, oxygenMoleFraction - step);
    const Sample upper = evaluateSample(
        kernel, pressure, temperature, oxygenMoleFraction + step);
    if (!lower.valid || !upper.valid) return false;
    for (std::size_t species = 0; species < 2; ++species)
    {
        derivative[species] =
            (upper.reducedChemicalPotential[species]
                - lower.reducedChemicalPotential[species])/(2.0*step);
        if (!std::isfinite(derivative[species])) return false;
    }
    return true;
}

RefinedTie refineTie(
    const PRMixtureKernel& kernel,
    double pressure,
    double temperature,
    const Sample& sampledLeft,
    const Sample& sampledRight)
{
    RefinedTie result;
    const double initialSpan = sampledRight.oxygenMoleFraction
        - sampledLeft.oxygenMoleFraction;
    if (!sampledLeft.valid || !sampledRight.valid
        || !std::isfinite(initialSpan) || initialSpan <= 0.0)
        return result;

    constexpr double endpointMargin = 1.0e-10;
    const double minimumSpan = 0.5*initialSpan;
    double leftX = std::max(
        endpointMargin, sampledLeft.oxygenMoleFraction);
    double rightX = std::min(
        1.0 - endpointMargin, sampledRight.oxygenMoleFraction);

    for (std::size_t iteration = 0; iteration <= 100; ++iteration)
    {
        const Sample left = evaluateSample(
            kernel, pressure, temperature, leftX);
        const Sample right = evaluateSample(
            kernel, pressure, temperature, rightX);
        if (!left.valid || !right.valid) return result;
        const double residual = fugacityResidual(left, right);
        if (residual <= 1.0e-10)
        {
            result.success = true;
            result.left = left;
            result.right = right;
            result.fugacityResidual = residual;
            result.iterations = iteration;
            return result;
        }

        std::array<double, 2> leftDerivative;
        std::array<double, 2> rightDerivative;
        if (!chemicalPotentialDerivative(
                kernel, pressure, temperature, leftX, leftDerivative)
            || !chemicalPotentialDerivative(
                kernel, pressure, temperature, rightX, rightDerivative))
            return result;

        const double f0 = left.reducedChemicalPotential[0]
            - right.reducedChemicalPotential[0];
        const double f1 = left.reducedChemicalPotential[1]
            - right.reducedChemicalPotential[1];
        const double a = leftDerivative[0];
        const double b = -rightDerivative[0];
        const double c = leftDerivative[1];
        const double d = -rightDerivative[1];
        const double determinant = a*d - b*c;
        if (!std::isfinite(determinant)
            || std::abs(determinant) <= 1.0e-14)
            return result;
        const double deltaLeft = (-f0*d + b*f1)/determinant;
        const double deltaRight = (c*f0 - a*f1)/determinant;
        if (!std::isfinite(deltaLeft) || !std::isfinite(deltaRight))
            return result;

        bool accepted = false;
        double damping = 1.0;
        for (std::size_t lineSearch = 0; lineSearch < 30; ++lineSearch)
        {
            const double candidateLeft = leftX + damping*deltaLeft;
            const double candidateRight = rightX + damping*deltaRight;
            if (candidateLeft > endpointMargin
                && candidateRight < 1.0 - endpointMargin
                && candidateRight - candidateLeft >= minimumSpan)
            {
                const Sample candidateLeftState = evaluateSample(
                    kernel, pressure, temperature, candidateLeft);
                const Sample candidateRightState = evaluateSample(
                    kernel, pressure, temperature, candidateRight);
                const double candidateResidual = fugacityResidual(
                    candidateLeftState, candidateRightState);
                if (candidateLeftState.valid && candidateRightState.valid
                    && candidateResidual < residual)
                {
                    leftX = candidateLeft;
                    rightX = candidateRight;
                    accepted = true;
                    break;
                }
            }
            damping *= 0.5;
        }
        if (!accepted) return result;
    }
    return result;
}

std::vector<std::size_t> lowerConvexHull(
    const std::vector<Sample>& samples)
{
    std::vector<std::size_t> hull;
    hull.reserve(samples.size());
    for (std::size_t index = 0; index < samples.size(); ++index)
    {
        if (!samples[index].valid) continue;
        while (hull.size() >= 2)
        {
            const Sample& first = samples[hull[hull.size() - 2]];
            const Sample& second = samples[hull[hull.size() - 1]];
            const Sample& third = samples[index];
            if (slope(second, third) > slope(first, second) + slopeTolerance)
                break;
            hull.pop_back();
        }
        hull.push_back(index);
    }
    return hull;
}

std::size_t contactClusters(
    const std::vector<Sample>& samples,
    std::size_t left,
    std::size_t right)
{
    const double tieSlope = slope(samples[left], samples[right]);
    const double intercept = samples[left].reducedGibbs
        - tieSlope*samples[left].oxygenMoleFraction;
    std::size_t clusters = 0;
    bool inCluster = false;
    for (std::size_t index = left; index <= right; ++index)
    {
        if (!samples[index].valid)
        {
            inCluster = false;
            continue;
        }
        const double line = tieSlope*samples[index].oxygenMoleFraction
            + intercept;
        const double gap = samples[index].reducedGibbs - line;
        const bool contact = gap <= contactTolerance;
        if (contact && !inCluster) ++clusters;
        inCluster = contact;
    }
    return clusters;
}

void classifyTieFromRoots(
    StableRootKind left,
    StableRootKind right,
    TopologyRow& row)
{
    if (left == StableRootKind::LIQUID
        && right == StableRootKind::LIQUID)
    {
        ++row.lleSegments;
        return;
    }
    if ((left == StableRootKind::LIQUID
            && (right == StableRootKind::VAPOR
                || right == StableRootKind::SINGLE))
        || (right == StableRootKind::LIQUID
            && (left == StableRootKind::VAPOR
                || left == StableRootKind::SINGLE)))
    {
        ++row.vleSegments;
        return;
    }
    ++row.otherSegments;
}

TPTieDisposition tpTieDisposition(
    const PRTPFlash& flash,
    double pressure,
    double temperature,
    double leftOxygen,
    double rightOxygen,
    double& auditedFeedOxygen,
    double& endpointError,
    TPFlashResult& flashResult)
{
    constexpr double feedFractions[] = {0.10, 0.25, 0.50, 0.75, 0.90};
    endpointError = 0.0;
    for (double feedFraction : feedFractions)
    {
        const double feedOxygen = leftOxygen
            + feedFraction*(rightOxygen - leftOxygen);
        auditedFeedOxygen = feedOxygen;
        const TPFlashResult split = flash.flashMoleFractions(
            pressure, temperature, {feedOxygen, 1.0 - feedOxygen});
        flashResult = split;
        if (!split.valid()) return TPTieDisposition::FLASH_FAILURE;
        if (split.region != EquilibriumRegion::TWO_PHASE)
            return TPTieDisposition::NON_TWO_PHASE;
        if (split.liquidMoleFractions.size() != 2
            || split.vaporMoleFractions.size() != 2)
            return TPTieDisposition::FLASH_FAILURE;

        const double splitLeft = std::min(
            split.liquidMoleFractions[0], split.vaporMoleFractions[0]);
        const double splitRight = std::max(
            split.liquidMoleFractions[0], split.vaporMoleFractions[0]);
        const double currentEndpointError = std::max(
            std::abs(splitLeft - leftOxygen),
            std::abs(splitRight - rightOxygen));
        endpointError = std::max(endpointError, currentEndpointError);
        if (currentEndpointError > tpEndpointTolerance)
            return TPTieDisposition::DIFFERENT_SPLIT;
    }
    return TPTieDisposition::MATCH;
}

TopologyRow analyze(
    double pressure,
    double temperature,
    double kijValue,
    std::size_t compositionIntervals)
{
    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    BinaryInteractionDB kij;
    kij.add("O2", "C3H8", kijValue);
    const PRMixtureKernel kernel(
        critical, kij, {"O2", "C3H8"});
    const PRTPFlash flash(
        critical, kij, {"O2", "C3H8"});
    if (!kernel.configured())
        throw std::runtime_error("kernel-configuration");
    if (flash.configurationStatus() != PhaseEquilibriumStatus::SUCCESS)
        throw std::runtime_error("flash-configuration");

    std::vector<Sample> samples;
    samples.reserve(compositionIntervals + 1);
    TopologyRow row;
    for (std::size_t index = 0; index <= compositionIntervals; ++index)
    {
        const double oxygenMoleFraction = static_cast<double>(index)
            / static_cast<double>(compositionIntervals);
        Sample sample = evaluateSample(
            kernel, pressure, temperature, oxygenMoleFraction);
        if (sample.valid) ++row.validSamples;
        else ++row.invalidSamples;
        samples.push_back(std::move(sample));
    }

    const std::vector<std::size_t> hull = lowerConvexHull(samples);
    row.hullVertices = hull.size();
    std::vector<TieEvidence> ties;
    ties.reserve(hull.size());
    for (std::size_t segment = 1; segment < hull.size(); ++segment)
    {
        const std::size_t left = hull[segment - 1];
        const std::size_t right = hull[segment];
        const std::size_t skipped = right - left - 1;
        if (skipped < minimumSkippedSamples) continue;

        ++row.tieSegments;
        row.maximumSkippedSamples = std::max(
            row.maximumSkippedSamples, skipped);
        const RefinedTie refined = refineTie(
            kernel, pressure, temperature, samples[left], samples[right]);
        const Sample& tieLeft = refined.success
            ? refined.left : samples[left];
        const Sample& tieRight = refined.success
            ? refined.right : samples[right];
        if (refined.success) ++row.refinedTieSegments;
        else ++row.unrefinedTieSegments;

        TieEvidence evidence;
        evidence.left = tieLeft;
        evidence.right = tieRight;
        evidence.refined = refined.success;
        evidence.contactClusters = contactClusters(samples, left, right);
        row.maximumContactClusters = std::max(
            row.maximumContactClusters, evidence.contactClusters);
        ties.push_back(evidence);

        const double tieSpan = tieRight.oxygenMoleFraction
            - tieLeft.oxygenMoleFraction;
        if (tieSpan > row.widestTieSpan)
        {
            row.widestTieSpan = tieSpan;
            row.widestTieLeftOxygen = tieLeft.oxygenMoleFraction;
            row.widestTieRightOxygen = tieRight.oxygenMoleFraction;
            row.widestTieLeftZ = tieLeft.compressibility;
            row.widestTieRightZ = tieRight.compressibility;
            row.widestTieFugacityResidual = fugacityResidual(
                tieLeft, tieRight);
            const StabilityResult leftStability = flash.analyzeStability(
                pressure, temperature,
                {tieLeft.oxygenMoleFraction,
                    1.0 - tieLeft.oxygenMoleFraction});
            row.tieLeftStabilityStatus = statusName(leftStability.status);
            row.tieLeftStabilityAttempted = leftStability.attemptedTrials;
            row.tieLeftStabilityConverged = leftStability.convergedTrials;
            row.tieLeftStabilityFailed = leftStability.failedTrials;
            row.tieLeftMinimumTpd = leftStability.minimumTpd;
            const StabilityResult rightStability = flash.analyzeStability(
                pressure, temperature,
                {tieRight.oxygenMoleFraction,
                    1.0 - tieRight.oxygenMoleFraction});
            row.tieRightStabilityStatus = statusName(rightStability.status);
            row.tieRightStabilityAttempted = rightStability.attemptedTrials;
            row.tieRightStabilityConverged = rightStability.convergedTrials;
            row.tieRightStabilityFailed = rightStability.failedTrials;
            row.tieRightMinimumTpd = rightStability.minimumTpd;
        }
    }

    std::vector<bool> threePhaseTie(ties.size(), false);
    for (std::size_t tie = 0; tie < ties.size(); ++tie)
    {
        if (ties[tie].contactClusters >= 3)
        {
            threePhaseTie[tie] = true;
            row.possibleThreePhase = true;
        }
    }
    for (std::size_t tie = 1; tie < ties.size(); ++tie)
    {
        const TieEvidence& previous = ties[tie - 1];
        const TieEvidence& current = ties[tie];
        if (!previous.refined || !current.refined) continue;

        const double sharedPhaseError = std::abs(
            previous.right.oxygenMoleFraction
            - current.left.oxygenMoleFraction);
        const double slopeError = std::abs(
            slope(previous.left, previous.right)
            - slope(current.left, current.right));
        const double interceptError = std::abs(
            intercept(previous.left, previous.right)
            - intercept(current.left, current.right));
        if (sharedPhaseError <= sharedPhaseTolerance
            && slopeError <= commonTangentTolerance
            && interceptError <= commonTangentTolerance)
        {
            ++row.threePhaseTangentPairs;
            threePhaseTie[tie - 1] = true;
            threePhaseTie[tie] = true;
            row.possibleThreePhase = true;
            updateMaximum(
                row.maximumThreePhaseSharedPhaseError, sharedPhaseError);
            updateMaximum(row.maximumThreePhaseSlopeError, slopeError);
            updateMaximum(row.maximumThreePhaseInterceptError, interceptError);
        }
    }

    for (std::size_t tie = 0; tie < ties.size(); ++tie)
    {
        const Sample& tieLeft = ties[tie].left;
        const Sample& tieRight = ties[tie].right;
        if (threePhaseTie[tie])
        {
            ++row.threePhaseTieSegments;
            classifyTieFromRoots(
                tieLeft.rootKind, tieRight.rootKind, row);
            continue;
        }

        double tpFeedOxygen = std::numeric_limits<double>::quiet_NaN();
        double tpEndpointError = std::numeric_limits<double>::infinity();
        TPFlashResult tpResult;
        const TPTieDisposition tpDisposition = tpTieDisposition(
                flash, pressure, temperature,
                tieLeft.oxygenMoleFraction,
                tieRight.oxygenMoleFraction,
                tpFeedOxygen,
                tpEndpointError, tpResult);
        row.tpFlashIterations = tpResult.flashIterations;
        row.tpMaximumFugacityResidual =
            tpResult.maximumLogFugacityResidual;
        row.tpMaximumMaterialResidual =
            tpResult.maximumMoleMaterialResidual;
        row.tpLiquidMinimumTpd = tpResult.liquidMinimumTpd;
        row.tpVaporMinimumTpd = tpResult.vaporMinimumTpd;
        if (std::isfinite(tpEndpointError))
        {
            row.maximumTpEndpointError = std::max(
                row.maximumTpEndpointError, tpEndpointError);
        }
        if (tpDisposition == TPTieDisposition::MATCH)
        {
            ++row.tpSupportedTieSegments;
            ++row.vleSegments;
        }
        else
        {
            ++row.tpDisagreements;
            if (tpDisposition == TPTieDisposition::FLASH_FAILURE)
            {
                ++row.tpFlashFailures;
                row.tpFailureStatus = statusName(tpResult.status);
                const StabilityResult stability = flash.analyzeStability(
                    pressure,
                    temperature,
                    {tpFeedOxygen, 1.0 - tpFeedOxygen});
                row.tpStabilityAttempted = stability.attemptedTrials;
                row.tpStabilityConverged = stability.convergedTrials;
                row.tpStabilityFailed = stability.failedTrials;
            }
            else if (tpDisposition == TPTieDisposition::NON_TWO_PHASE)
                ++row.tpNonTwoPhase;
            else
                ++row.tpDifferentSplits;
            classifyTieFromRoots(
                tieLeft.rootKind, tieRight.rootKind, row);
        }
    }
    return row;
}

std::size_t parseIntervals(const char* text)
{
    const std::string input(text);
    std::size_t consumed = 0;
    const unsigned long parsed = std::stoul(input, &consumed);
    if (consumed != input.size())
        throw std::invalid_argument("composition-intervals");
    if (parsed < 100 || parsed > 20000)
        throw std::out_of_range("composition-intervals");
    return static_cast<std::size_t>(parsed);
}

double parseFinite(const char* text, bool positive)
{
    const std::string input(text);
    std::size_t consumed = 0;
    const double parsed = std::stod(input, &consumed);
    if (consumed != input.size())
        throw std::invalid_argument("finite-number");
    if (!std::isfinite(parsed) || (positive && parsed <= 0.0))
        throw std::out_of_range("finite-number");
    return parsed;
}

void writeRow(
    double pressure,
    double temperature,
    double kijValue,
    std::size_t compositionIntervals,
    const TopologyRow& row)
{
    std::cout << pressure << ',' << temperature << ',' << kijValue << ','
              << compositionIntervals << ',' << row.validSamples << ','
              << row.invalidSamples << ',' << row.hullVertices << ','
              << row.tieSegments << ',' << row.lleSegments << ','
              << row.vleSegments << ',' << row.otherSegments << ','
              << row.maximumContactClusters << ','
              << row.threePhaseTangentPairs << ','
              << row.threePhaseTieSegments << ','
              << row.maximumThreePhaseSharedPhaseError << ','
              << row.maximumThreePhaseSlopeError << ','
              << row.maximumThreePhaseInterceptError << ','
              << row.maximumSkippedSamples << ','
              << row.refinedTieSegments << ','
              << row.unrefinedTieSegments << ','
              << row.tpSupportedTieSegments << ',' << row.tpDisagreements
              << ',' << row.tpFlashFailures << ',' << row.tpNonTwoPhase
              << ',' << row.tpDifferentSplits << ',' << row.tpFailureStatus
              << ',' << row.tpStabilityAttempted
              << ',' << row.tpStabilityConverged
              << ',' << row.tpStabilityFailed
              << ',' << row.tpFlashIterations
              << ',' << row.tpMaximumFugacityResidual
              << ',' << row.tpMaximumMaterialResidual
              << ',' << row.tpLiquidMinimumTpd
              << ',' << row.tpVaporMinimumTpd
              << ',' << row.maximumTpEndpointError << ','
              << row.widestTieSpan << ',' << row.widestTieLeftOxygen << ','
              << row.widestTieRightOxygen << ',' << row.widestTieLeftZ << ','
              << row.widestTieRightZ << ','
              << row.widestTieFugacityResidual << ','
              << row.tieLeftStabilityStatus << ','
              << row.tieLeftStabilityAttempted << ','
              << row.tieLeftStabilityConverged << ','
              << row.tieLeftStabilityFailed << ','
              << row.tieLeftMinimumTpd << ','
              << row.tieRightStabilityStatus << ','
              << row.tieRightStabilityAttempted << ','
              << row.tieRightStabilityConverged << ','
              << row.tieRightStabilityFailed << ','
              << row.tieRightMinimumTpd << ','
              << (row.possibleThreePhase ? 1 : 0) << ",0\n";
}

} // namespace

int main(int argc, char** argv)
{
    std::size_t compositionIntervals = 2000;
    std::optional<SingleState> singleState;
    if (argc == 3 && std::string(argv[1]) == "--composition-intervals")
    {
        try
        {
            compositionIntervals = parseIntervals(argv[2]);
        }
        catch (const std::exception&)
        {
            std::cerr << "PR_BINARY_PHASE_TOPOLOGY_SCREEN_FAIL reason="
                      << "invalid-composition-intervals\n";
            return 2;
        }
    }
    else if ((argc == 5 || argc == 7)
        && std::string(argv[1]) == "--state")
    {
        try
        {
            SingleState state;
            state.pressure = parseFinite(argv[2], true);
            state.temperature = parseFinite(argv[3], true);
            state.kij = parseFinite(argv[4], false);
            if (argc == 7)
            {
                if (std::string(argv[5]) != "--composition-intervals")
                    throw std::invalid_argument("state-option");
                compositionIntervals = parseIntervals(argv[6]);
            }
            singleState = state;
        }
        catch (const std::exception&)
        {
            std::cerr << "PR_BINARY_PHASE_TOPOLOGY_SCREEN_FAIL reason="
                      << "invalid-state\n";
            return 2;
        }
    }
    else if (argc != 1)
    {
        std::cerr << "PR_BINARY_PHASE_TOPOLOGY_SCREEN_FAIL reason="
                  << "invalid-arguments\n";
        return 2;
    }

    const std::array<double, 3> pressures = {4.0e6, 5.0e6, 6.0e6};
    const std::array<double, 5> kijValues = {
        -0.05, 0.0, 0.055, 0.10, 0.15};
    const std::array<double, 18> temperatures = {
        110.0, 120.0, 130.0, 140.0, 150.0, 160.0,
        170.0, 180.0, 190.0, 200.0, 210.0, 220.0,
        230.0, 240.0, 250.0, 260.0, 280.0, 300.0};

    std::cout.imbue(std::locale::classic());
    std::cout << std::setprecision(17);
    std::cout << "# " << schema << '\n';
    std::cout
        << "pressure_pa,temperature_k,kij,composition_intervals,"
        << "valid_samples,invalid_samples,hull_vertices,tie_segments,"
        << "lle_segments,vle_segments,other_segments,"
        << "maximum_contact_clusters,three_phase_tangent_pairs,"
        << "three_phase_tie_segments,maximum_three_phase_shared_phase_error,"
        << "maximum_three_phase_slope_error,"
        << "maximum_three_phase_intercept_error,maximum_skipped_samples,"
        << "refined_tie_segments,unrefined_tie_segments,"
        << "tp_supported_tie_segments,tp_disagreements,"
        << "tp_flash_failures,tp_non_two_phase,tp_different_splits,"
        << "tp_failure_status,"
        << "tp_stability_attempted,tp_stability_converged,"
        << "tp_stability_failed,"
        << "tp_flash_iterations,tp_maximum_fugacity_residual,"
        << "tp_maximum_material_residual,tp_liquid_minimum_tpd,"
        << "tp_vapor_minimum_tpd,"
        << "maximum_tp_endpoint_error,widest_tie_span,"
        << "widest_tie_left_x_o2,widest_tie_right_x_o2,"
        << "widest_tie_left_z,widest_tie_right_z,"
        << "widest_tie_fugacity_residual,"
        << "tie_left_stability_status,tie_left_stability_attempted,"
        << "tie_left_stability_converged,tie_left_stability_failed,"
        << "tie_left_minimum_tpd,tie_right_stability_status,"
        << "tie_right_stability_attempted,tie_right_stability_converged,"
        << "tie_right_stability_failed,tie_right_minimum_tpd,"
        << "possible_three_phase,production_ready\n";

    std::size_t rows = 0;
    std::size_t rowsWithTies = 0;
    std::size_t possibleThreePhaseRows = 0;
    std::size_t invalidRows = 0;
    try
    {
        const auto processState = [&](
            double pressure, double temperature, double kijValue)
        {
            const TopologyRow row = analyze(
                pressure, temperature, kijValue,
                compositionIntervals);
            writeRow(
                pressure, temperature, kijValue,
                compositionIntervals, row);
            ++rows;
            if (row.tieSegments > 0) ++rowsWithTies;
            if (row.possibleThreePhase) ++possibleThreePhaseRows;
            if (row.invalidSamples > 0) ++invalidRows;
        };

        if (singleState.has_value())
        {
            processState(
                singleState->pressure,
                singleState->temperature,
                singleState->kij);
        }
        else
        {
            for (double pressure : pressures)
            {
                for (double temperature : temperatures)
                {
                    for (double kijValue : kijValues)
                        processState(pressure, temperature, kijValue);
                }
            }
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << "PR_BINARY_PHASE_TOPOLOGY_SCREEN_FAIL reason="
                  << error.what() << '\n';
        return 1;
    }

    std::cout << "# PR_BINARY_PHASE_TOPOLOGY_SCREEN_BLOCKED rows="
              << rows << " rows_with_ties=" << rowsWithTies
              << " possible_three_phase_rows=" << possibleThreePhaseRows
              << " rows_with_invalid_samples=" << invalidRows
              << " composition_intervals=" << compositionIntervals
              << " phase_topology_validated=0 production_ready=0\n";
    return 0;
}
