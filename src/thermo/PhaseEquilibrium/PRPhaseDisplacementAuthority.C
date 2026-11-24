#include "PRPhaseDisplacementAuthority.H"

#include "PRPhaseAuthorityIdentity.H"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace Foam
{
namespace phaseEquilibrium
{

namespace
{

constexpr int maximumDyadicDepth = 12;
constexpr double compositionNormalizationTolerance = 1.0e-10;

bool finitePositive(double value)
{
    return std::isfinite(value) && value > 0.0;
}

double relativeResidual(double lhs, double rhs)
{
    return std::abs(lhs - rhs)
        /std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

bool finiteWithin(double value, double lower, double upper)
{
    return std::isfinite(value) && value >= lower && value <= upper;
}

bool validBounds(const PRPhaseDisplacementBounds& bounds)
{
    return finiteWithin(
            bounds.maximumDensityRelativeDisplacement, 0.0, 0.5)
        && bounds.maximumDensityRelativeDisplacement > 0.0
        && finiteWithin(
            bounds.maximumEnergyRelativeDisplacement, 0.0, 0.5)
        && bounds.maximumEnergyRelativeDisplacement > 0.0
        && finiteWithin(
            bounds.maximumCompositionL1Displacement, 0.0, 0.5)
        && bounds.maximumCompositionL1Displacement > 0.0
        && bounds.dyadicDepth >= 1
        && bounds.dyadicDepth <= maximumDyadicDepth
        && finiteWithin(bounds.minimumAdjacentRootGap, 0.0, 1.0)
        && bounds.minimumAdjacentRootGap > 0.0
        && finiteWithin(
            bounds.minimumDimensionlessDpDrhoMargin, 0.0, 1.0e6)
        && bounds.minimumDimensionlessDpDrhoMargin > 0.0
        && finiteWithin(
            bounds.maximumCompressibilityRelativeStep, 0.0, 1.0)
        && bounds.maximumCompressibilityRelativeStep > 0.0
        && finiteWithin(bounds.maximumPressureRelativeStep, 0.0, 1.0)
        && bounds.maximumPressureRelativeStep > 0.0
        && finiteWithin(
            bounds.maximumTemperatureRelativeStep, 0.0, 1.0)
        && bounds.maximumTemperatureRelativeStep > 0.0
        && finiteWithin(
            bounds.densityClosureRelativeTolerance, 0.0, 1.0e-6)
        && bounds.densityClosureRelativeTolerance > 0.0
        && finiteWithin(
            bounds.energyClosureAbsoluteTolerance, 0.0, 1.0)
        && finiteWithin(
            bounds.energyClosureRelativeTolerance, 0.0, 1.0e-6)
        && bounds.energyClosureAbsoluteTolerance
            + bounds.energyClosureRelativeTolerance > 0.0
        && finiteWithin(
            bounds.compositionClosureAbsoluteTolerance, 0.0, 1.0e-8)
        && bounds.compositionClosureAbsoluteTolerance > 0.0;
}

double maximumCompositionResidual(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs)
{
    if (lhs.size() != rhs.size())
        return std::numeric_limits<double>::infinity();
    double maximum = 0.0;
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        if (!std::isfinite(lhs[i]) || !std::isfinite(rhs[i]))
            return std::numeric_limits<double>::infinity();
        maximum = std::max(maximum, std::abs(lhs[i] - rhs[i]));
    }
    return maximum;
}

bool closesRequest(
    const PRPhaseRhoEResult& recovered,
    const PRPhaseRhoERequest& request,
    const PRPhaseDisplacementBounds& bounds,
    double& densityResidual,
    double& energyResidual,
    double& compositionResidual)
{
    if (!recovered.success() || !recovered.state) return false;
    densityResidual = relativeResidual(
        recovered.state->rho, request.density);
    energyResidual = relativeResidual(
        recovered.state->e, request.specificInternalEnergy);
    compositionResidual = maximumCompositionResidual(
        recovered.state->massFractions, request.massFractions);
    const double energyTolerance = bounds.energyClosureAbsoluteTolerance
        + bounds.energyClosureRelativeTolerance
            *std::max(std::abs(request.specificInternalEnergy), 1.0);
    return densityResidual <= bounds.densityClosureRelativeTolerance
        && std::abs(
            recovered.state->e - request.specificInternalEnergy)
            <= energyTolerance
        && compositionResidual
            <= bounds.compositionClosureAbsoluteTolerance;
}

bool admitsSignature(
    const PRPhaseRootSignature& signature,
    const PRPhaseDisplacementBounds& bounds)
{
    return signature.valid()
        && signature.dimensionlessDpDrhoMargin
            >= bounds.minimumDimensionlessDpDrhoMargin
        && (!signature.hasAdjacentRoot
            || signature.adjacentRootGap
                >= bounds.minimumAdjacentRootGap);
}

void observeSignature(
    const PRPhaseRootSignature& signature,
    PRPhaseDisplacementDiagnostics& diagnostics)
{
    diagnostics.minimumObservedDimensionlessDpDrhoMargin = std::min(
        diagnostics.minimumObservedDimensionlessDpDrhoMargin,
        signature.dimensionlessDpDrhoMargin);
    if (signature.hasAdjacentRoot)
    {
        diagnostics.minimumObservedAdjacentRootGap = std::min(
            diagnostics.minimumObservedAdjacentRootGap,
            signature.adjacentRootGap);
    }
}

PRPhaseRhoERequest interpolateRequest(
    const PRPhaseRhoERequest& parent,
    const PRPhaseRhoERequest& requested,
    std::size_t numerator,
    std::size_t denominator,
    double temperatureSeed)
{
    if (numerator == denominator)
    {
        return PRPhaseRhoERequest(
            requested.density,
            requested.specificInternalEnergy,
            requested.massFractions,
            temperatureSeed);
    }
    const double fraction = static_cast<double>(numerator)
        /static_cast<double>(denominator);
    std::vector<double> massFractions(parent.massFractions.size(), 0.0);
    for (std::size_t i = 0; i < massFractions.size(); ++i)
    {
        massFractions[i] = parent.massFractions[i]
            + fraction*(requested.massFractions[i]
                - parent.massFractions[i]);
    }
    return PRPhaseRhoERequest(
        parent.density + fraction*(requested.density - parent.density),
        parent.specificInternalEnergy
            + fraction*(requested.specificInternalEnergy
                - parent.specificInternalEnergy),
        std::move(massFractions),
        temperatureSeed);
}

} // End anonymous namespace

const char* phaseDisplacementStatusName(
    PRPhaseDisplacementStatus status)
{
    switch (status)
    {
        case PRPhaseDisplacementStatus::SUCCESS: return "SUCCESS";
        case PRPhaseDisplacementStatus::INVALID_CONFIGURATION:
            return "INVALID_CONFIGURATION";
        case PRPhaseDisplacementStatus::INVALID_PARENT_REQUEST:
            return "INVALID_PARENT_REQUEST";
        case PRPhaseDisplacementStatus::INVALID_DISPLACED_REQUEST:
            return "INVALID_DISPLACED_REQUEST";
        case PRPhaseDisplacementStatus::PARENT_RECOVERY_FAILED:
            return "PARENT_RECOVERY_FAILED";
        case PRPhaseDisplacementStatus::PARENT_CLOSURE_FAILED:
            return "PARENT_CLOSURE_FAILED";
        case PRPhaseDisplacementStatus::DISPLACEMENT_OUT_OF_BOUNDS:
            return "DISPLACEMENT_OUT_OF_BOUNDS";
        case PRPhaseDisplacementStatus::ROOT_SIGNATURE_FAILED:
            return "ROOT_SIGNATURE_FAILED";
        case PRPhaseDisplacementStatus::ROOT_ADMISSION_FAILED:
            return "ROOT_ADMISSION_FAILED";
        case PRPhaseDisplacementStatus::WAYPOINT_RECOVERY_FAILED:
            return "WAYPOINT_RECOVERY_FAILED";
        case PRPhaseDisplacementStatus::WAYPOINT_CLOSURE_FAILED:
            return "WAYPOINT_CLOSURE_FAILED";
        case PRPhaseDisplacementStatus::ROOT_CONTINUITY_FAILED:
            return "ROOT_CONTINUITY_FAILED";
        case PRPhaseDisplacementStatus::REQUEST_CLOSURE_FAILED:
            return "REQUEST_CLOSURE_FAILED";
    }
    return "UNKNOWN_STATUS";
}

PRPhaseRhoERequest::PRPhaseRhoERequest(
    double requestDensity,
    double requestSpecificInternalEnergy,
    std::vector<double> requestMassFractions,
    double requestTemperatureSeed)
:
    density(requestDensity),
    specificInternalEnergy(requestSpecificInternalEnergy),
    massFractions(std::move(requestMassFractions)),
    temperatureSeed(requestTemperatureSeed)
{}

PRPhaseDisplacementBounds::PRPhaseDisplacementBounds(
    double maximumDensityRelativeDisplacementValue,
    double maximumEnergyRelativeDisplacementValue,
    double maximumCompositionL1DisplacementValue,
    int dyadicDepthValue,
    double minimumAdjacentRootGapValue,
    double minimumDimensionlessDpDrhoMarginValue,
    double maximumCompressibilityRelativeStepValue,
    double maximumPressureRelativeStepValue,
    double maximumTemperatureRelativeStepValue,
    double densityClosureRelativeToleranceValue,
    double energyClosureAbsoluteToleranceValue,
    double energyClosureRelativeToleranceValue,
    double compositionClosureAbsoluteToleranceValue)
:
    maximumDensityRelativeDisplacement(
        maximumDensityRelativeDisplacementValue),
    maximumEnergyRelativeDisplacement(
        maximumEnergyRelativeDisplacementValue),
    maximumCompositionL1Displacement(
        maximumCompositionL1DisplacementValue),
    dyadicDepth(dyadicDepthValue),
    minimumAdjacentRootGap(minimumAdjacentRootGapValue),
    minimumDimensionlessDpDrhoMargin(
        minimumDimensionlessDpDrhoMarginValue),
    maximumCompressibilityRelativeStep(
        maximumCompressibilityRelativeStepValue),
    maximumPressureRelativeStep(maximumPressureRelativeStepValue),
    maximumTemperatureRelativeStep(maximumTemperatureRelativeStepValue),
    densityClosureRelativeTolerance(densityClosureRelativeToleranceValue),
    energyClosureAbsoluteTolerance(energyClosureAbsoluteToleranceValue),
    energyClosureRelativeTolerance(energyClosureRelativeToleranceValue),
    compositionClosureAbsoluteTolerance(
        compositionClosureAbsoluteToleranceValue)
{}

PRPhaseDisplacementAuthority::PRPhaseDisplacementAuthority(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    const PRPhaseRhoEOptions& recoveryOptions,
    const PRPhaseDisplacementBounds& bounds)
:
    criticalData_(criticalData),
    kij_(kij),
    speciesNames_(speciesNames),
    recoveryOptions_(recoveryOptions),
    physicalModelAuthorityFingerprint_(
        makePRPhasePhysicalModelAuthorityFingerprint(
            criticalData, kij, speciesNames)),
    authorityFingerprint_(makePRPhaseRhoEAuthorityFingerprint(
        criticalData, kij, speciesNames, recoveryOptions)),
    displacementAuthorityFingerprint_(
        makePRPhaseDisplacementAuthorityFingerprint(
            authorityFingerprint_, bounds)),
    bounds_(bounds),
    recovery_(criticalData_, kij_, speciesNames_, recoveryOptions_),
    signatureEvaluator_(
        criticalData_, kij_, speciesNames_,
        recoveryOptions_.rootRelativeTolerance),
    configurationStatus_(PRPhaseDisplacementStatus::INVALID_CONFIGURATION)
{
    if (recovery_.configured() && signatureEvaluator_.configured()
        && validBounds(bounds_)
        && recovery_.physicalModelAuthorityFingerprint()
            == physicalModelAuthorityFingerprint_
        && recovery_.authorityFingerprint() == authorityFingerprint_
        && validPRPhaseAuthorityFingerprint(
            physicalModelAuthorityFingerprint_)
        && validPRPhaseAuthorityFingerprint(authorityFingerprint_)
        && validPRPhaseAuthorityFingerprint(
            displacementAuthorityFingerprint_))
    {
        configurationStatus_ = PRPhaseDisplacementStatus::SUCCESS;
    }
}

bool PRPhaseDisplacementAuthority::validRequest(
    const PRPhaseRhoERequest& request) const
{
    if (!finitePositive(request.density)
        || !std::isfinite(request.specificInternalEnergy)
        || !finitePositive(request.temperatureSeed)
        || request.massFractions.size() != speciesNames_.size())
    {
        return false;
    }
    double sum = 0.0;
    for (double fraction : request.massFractions)
    {
        if (!std::isfinite(fraction) || fraction < 0.0) return false;
        sum += fraction;
    }
    return std::isfinite(sum)
        && std::abs(sum - 1.0) <= compositionNormalizationTolerance;
}

PRPhaseDisplacementResult PRPhaseDisplacementAuthority::evaluate(
    const PRPhaseRhoERequest& parent,
    const PRPhaseRhoERequest& requested,
    PRPhaseIdentity phaseIdentity,
    PRPhaseRootIdentity rootIdentity) const
{
    PRPhaseDisplacementResult result;
    result.status = configurationStatus_;
    if (!configured()) return result;
    if (!validRequest(parent))
    {
        result.status = PRPhaseDisplacementStatus::INVALID_PARENT_REQUEST;
        return result;
    }
    if (!validRequest(requested))
    {
        result.status =
            PRPhaseDisplacementStatus::INVALID_DISPLACED_REQUEST;
        return result;
    }

    result.diagnostics.densityRelativeDisplacement = relativeResidual(
        requested.density, parent.density);
    result.diagnostics.energyRelativeDisplacement = relativeResidual(
        requested.specificInternalEnergy,
        parent.specificInternalEnergy);
    result.diagnostics.compositionL1Displacement = 0.0;
    for (std::size_t i = 0; i < parent.massFractions.size(); ++i)
    {
        result.diagnostics.compositionL1Displacement += std::abs(
            requested.massFractions[i] - parent.massFractions[i]);
    }
    if (result.diagnostics.densityRelativeDisplacement
            > bounds_.maximumDensityRelativeDisplacement
        || result.diagnostics.energyRelativeDisplacement
            > bounds_.maximumEnergyRelativeDisplacement
        || result.diagnostics.compositionL1Displacement
            > bounds_.maximumCompositionL1Displacement)
    {
        result.status =
            PRPhaseDisplacementStatus::DISPLACEMENT_OUT_OF_BOUNDS;
        return result;
    }

    ++result.diagnostics.attemptedWaypoints;
    PRPhaseRhoEResult previous = recovery_.recover(
        parent.density,
        parent.specificInternalEnergy,
        parent.massFractions,
        phaseIdentity,
        rootIdentity,
        parent.temperatureSeed);
    result.diagnostics.recoveryStatus = previous.status;
    if (!previous.success())
    {
        result.status =
            PRPhaseDisplacementStatus::PARENT_RECOVERY_FAILED;
        return result;
    }
    double densityResidual = std::numeric_limits<double>::infinity();
    double energyResidual = std::numeric_limits<double>::infinity();
    double compositionResidual = std::numeric_limits<double>::infinity();
    if (!closesRequest(
            previous, parent, bounds_, densityResidual,
            energyResidual, compositionResidual))
    {
        result.status = PRPhaseDisplacementStatus::PARENT_CLOSURE_FAILED;
        return result;
    }
    ++result.diagnostics.acceptedWaypoints;

    PRPhaseRootSignature previousSignature =
        signatureEvaluator_.evaluate(previous);
    result.diagnostics.signatureStatus = previousSignature.status;
    if (!previousSignature.valid())
    {
        result.status =
            PRPhaseDisplacementStatus::ROOT_SIGNATURE_FAILED;
        return result;
    }
    observeSignature(previousSignature, result.diagnostics);
    if (!admitsSignature(previousSignature, bounds_))
    {
        result.status = PRPhaseDisplacementStatus::ROOT_ADMISSION_FAILED;
        return result;
    }

    const std::size_t denominator =
        std::size_t{1} << static_cast<unsigned>(bounds_.dyadicDepth);
    for (std::size_t index = 1; index <= denominator; ++index)
    {
        const double fraction = static_cast<double>(index)
            /static_cast<double>(denominator);
        const PRPhaseRhoERequest waypoint = interpolateRequest(
            parent, requested, index, denominator,
            previous.state->temperature);
        ++result.diagnostics.attemptedWaypoints;
        PRPhaseRhoEResult current = recovery_.recover(
            waypoint.density,
            waypoint.specificInternalEnergy,
            waypoint.massFractions,
            phaseIdentity,
            rootIdentity,
            waypoint.temperatureSeed);
        result.diagnostics.recoveryStatus = current.status;
        if (!current.success())
        {
            result.diagnostics.failedWaypoint = index;
            result.diagnostics.failedFraction = fraction;
            result.status =
                PRPhaseDisplacementStatus::WAYPOINT_RECOVERY_FAILED;
            return result;
        }

        densityResidual = std::numeric_limits<double>::infinity();
        energyResidual = std::numeric_limits<double>::infinity();
        compositionResidual = std::numeric_limits<double>::infinity();
        if (!closesRequest(
                current, waypoint, bounds_, densityResidual,
                energyResidual, compositionResidual))
        {
            result.diagnostics.failedWaypoint = index;
            result.diagnostics.failedFraction = fraction;
            result.status = index == denominator
                ? PRPhaseDisplacementStatus::REQUEST_CLOSURE_FAILED
                : PRPhaseDisplacementStatus::WAYPOINT_CLOSURE_FAILED;
            return result;
        }

        PRPhaseRootSignature signature =
            signatureEvaluator_.evaluate(current);
        result.diagnostics.signatureStatus = signature.status;
        if (!signature.valid())
        {
            result.diagnostics.failedWaypoint = index;
            result.diagnostics.failedFraction = fraction;
            result.status =
                PRPhaseDisplacementStatus::ROOT_SIGNATURE_FAILED;
            return result;
        }
        observeSignature(signature, result.diagnostics);

        const double compressibilityStep = relativeResidual(
            signature.retainedCompressibility,
            previousSignature.retainedCompressibility);
        const double pressureStep = relativeResidual(
            current.state->pressure, previous.state->pressure);
        const double temperatureStep = relativeResidual(
            current.state->temperature, previous.state->temperature);
        result.diagnostics.maximumObservedCompressibilityRelativeStep =
            std::max(
                result.diagnostics
                    .maximumObservedCompressibilityRelativeStep,
                compressibilityStep);
        result.diagnostics.maximumObservedPressureRelativeStep = std::max(
            result.diagnostics.maximumObservedPressureRelativeStep,
            pressureStep);
        result.diagnostics.maximumObservedTemperatureRelativeStep = std::max(
            result.diagnostics.maximumObservedTemperatureRelativeStep,
            temperatureStep);

        if (!admitsSignature(signature, bounds_)
            || signature.physicalRootCount
                != previousSignature.physicalRootCount
            || signature.retainedOrdinal
                != previousSignature.retainedOrdinal
            || signature.hasAdjacentRoot
                != previousSignature.hasAdjacentRoot
            || compressibilityStep
                > bounds_.maximumCompressibilityRelativeStep
            || pressureStep > bounds_.maximumPressureRelativeStep
            || temperatureStep > bounds_.maximumTemperatureRelativeStep)
        {
            result.diagnostics.failedWaypoint = index;
            result.diagnostics.failedFraction = fraction;
            result.status =
                PRPhaseDisplacementStatus::ROOT_CONTINUITY_FAILED;
            return result;
        }

        ++result.diagnostics.acceptedWaypoints;
        previous = std::move(current);
        previousSignature = signature;
    }

    result.diagnostics.requestDensityRelativeResidual = relativeResidual(
        previous.state->rho, requested.density);
    result.diagnostics.requestEnergyRelativeResidual = relativeResidual(
        previous.state->e, requested.specificInternalEnergy);
    result.diagnostics.requestCompositionMaximumAbsoluteResidual =
        maximumCompositionResidual(
            previous.state->massFractions, requested.massFractions);
    if (!closesRequest(
            previous, requested, bounds_, densityResidual,
            energyResidual, compositionResidual))
    {
        result.status = PRPhaseDisplacementStatus::REQUEST_CLOSURE_FAILED;
        return result;
    }

    result.status = PRPhaseDisplacementStatus::SUCCESS;
    result.requested = std::move(previous);
    result.rootSignature = previousSignature;
    return result;
}

PRPhaseDisplacementStatus PRPhaseDisplacementAuthority::tryCommit(
    const PRPhaseRhoERequest& parent,
    const PRPhaseRhoERequest& requested,
    PRPhaseIdentity phaseIdentity,
    PRPhaseRootIdentity rootIdentity,
    PRPhaseRhoEResult& committed,
    PRPhaseDisplacementDiagnostics& diagnostics) const
{
    PRPhaseDisplacementResult candidate = evaluate(
        parent, requested, phaseIdentity, rootIdentity);
    diagnostics = candidate.diagnostics;
    if (!candidate.success()) return candidate.status;
    PRPhaseRhoEResult publication = std::move(*candidate.requested);
    committed = std::move(publication);
    return PRPhaseDisplacementStatus::SUCCESS;
}

} // End namespace phaseEquilibrium
} // End namespace Foam
