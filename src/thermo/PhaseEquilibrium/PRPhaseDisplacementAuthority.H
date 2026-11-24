#ifndef PRPhaseDisplacementAuthority_H
#define PRPhaseDisplacementAuthority_H

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "PRPhaseRootSignature.H"

namespace Foam
{
namespace phaseEquilibrium
{

enum class PRPhaseDisplacementStatus
{
    SUCCESS,
    INVALID_CONFIGURATION,
    INVALID_PARENT_REQUEST,
    INVALID_DISPLACED_REQUEST,
    PARENT_RECOVERY_FAILED,
    PARENT_CLOSURE_FAILED,
    DISPLACEMENT_OUT_OF_BOUNDS,
    ROOT_SIGNATURE_FAILED,
    ROOT_ADMISSION_FAILED,
    WAYPOINT_RECOVERY_FAILED,
    WAYPOINT_CLOSURE_FAILED,
    ROOT_CONTINUITY_FAILED,
    REQUEST_CLOSURE_FAILED
};

const char* phaseDisplacementStatusName(
    PRPhaseDisplacementStatus status);

struct PRPhaseRhoERequest
{
    double density;
    double specificInternalEnergy;
    std::vector<double> massFractions;
    double temperatureSeed;

    PRPhaseRhoERequest(
        double requestDensity,
        double requestSpecificInternalEnergy,
        std::vector<double> requestMassFractions,
        double requestTemperatureSeed);
};

// Every admission and closure bound is explicit. There is intentionally no
// default constructor and no implicit production envelope.
struct PRPhaseDisplacementBounds
{
    const double maximumDensityRelativeDisplacement;
    const double maximumEnergyRelativeDisplacement;
    const double maximumCompositionL1Displacement;
    const int dyadicDepth;
    const double minimumAdjacentRootGap;
    const double minimumDimensionlessDpDrhoMargin;
    const double maximumCompressibilityRelativeStep;
    const double maximumPressureRelativeStep;
    const double maximumTemperatureRelativeStep;
    const double densityClosureRelativeTolerance;
    const double energyClosureAbsoluteTolerance;
    const double energyClosureRelativeTolerance;
    const double compositionClosureAbsoluteTolerance;

    PRPhaseDisplacementBounds() = delete;

    PRPhaseDisplacementBounds(
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
        double compositionClosureAbsoluteToleranceValue);
};

struct PRPhaseDisplacementDiagnostics
{
    std::size_t attemptedWaypoints = 0;
    std::size_t acceptedWaypoints = 0;
    std::size_t failedWaypoint = std::numeric_limits<std::size_t>::max();
    double failedFraction = std::numeric_limits<double>::quiet_NaN();
    PRPhaseRhoEStatus recoveryStatus =
        PRPhaseRhoEStatus::INVALID_CONFIGURATION;
    PRPhaseRootSignatureStatus signatureStatus =
        PRPhaseRootSignatureStatus::INVALID_CONFIGURATION;
    double densityRelativeDisplacement =
        std::numeric_limits<double>::infinity();
    double energyRelativeDisplacement =
        std::numeric_limits<double>::infinity();
    double compositionL1Displacement =
        std::numeric_limits<double>::infinity();
    double maximumObservedCompressibilityRelativeStep = 0.0;
    double maximumObservedPressureRelativeStep = 0.0;
    double maximumObservedTemperatureRelativeStep = 0.0;
    double minimumObservedAdjacentRootGap =
        std::numeric_limits<double>::infinity();
    double minimumObservedDimensionlessDpDrhoMargin =
        std::numeric_limits<double>::infinity();
    double requestDensityRelativeResidual =
        std::numeric_limits<double>::infinity();
    double requestEnergyRelativeResidual =
        std::numeric_limits<double>::infinity();
    double requestCompositionMaximumAbsoluteResidual =
        std::numeric_limits<double>::infinity();
};

struct PRPhaseDisplacementResult
{
    PRPhaseDisplacementStatus status =
        PRPhaseDisplacementStatus::INVALID_CONFIGURATION;
    std::optional<PRPhaseRhoEResult> requested;
    std::optional<PRPhaseRootSignature> rootSignature;
    PRPhaseDisplacementDiagnostics diagnostics;

    bool success() const
    {
        return status == PRPhaseDisplacementStatus::SUCCESS
            && requested.has_value() && requested->success()
            && rootSignature.has_value() && rootSignature->valid();
    }
};

// Offline transactional authority for a bounded displacement from one
// parent-owned phase-local rho-e-Y state. Homotopy requests are always formed
// directly from the immutable parent and final request at k/2^depth; no
// accumulated composition or energy updates are permitted.
class PRPhaseDisplacementAuthority
{
private:
    SpeciesCriticalDB criticalData_;
    BinaryInteractionDB kij_;
    std::vector<std::string> speciesNames_;
    PRPhaseRhoEOptions recoveryOptions_;
    std::string physicalModelAuthorityFingerprint_;
    std::string authorityFingerprint_;
    std::string displacementAuthorityFingerprint_;
    PRPhaseDisplacementBounds bounds_;
    PRPhaseRhoERecovery recovery_;
    PRPhaseRootSignatureEvaluator signatureEvaluator_;
    PRPhaseDisplacementStatus configurationStatus_;

    bool validRequest(const PRPhaseRhoERequest& request) const;

public:
    PRPhaseDisplacementAuthority(
        const SpeciesCriticalDB& criticalData,
        const BinaryInteractionDB& kij,
        const std::vector<std::string>& speciesNames,
        const PRPhaseRhoEOptions& recoveryOptions,
        const PRPhaseDisplacementBounds& bounds);

    bool configured() const
    {
        return configurationStatus_
            == PRPhaseDisplacementStatus::SUCCESS;
    }

    PRPhaseDisplacementStatus configurationStatus() const
    {
        return configurationStatus_;
    }

    const PRPhaseDisplacementBounds& bounds() const
    {
        return bounds_;
    }

    const std::string& physicalModelAuthorityFingerprint() const
    {
        return physicalModelAuthorityFingerprint_;
    }

    const std::string& recoveryAuthorityFingerprint() const
    {
        return authorityFingerprint_;
    }

    const std::string& authorityFingerprint() const
    {
        return displacementAuthorityFingerprint_;
    }

    PRPhaseDisplacementResult evaluate(
        const PRPhaseRhoERequest& parent,
        const PRPhaseRhoERequest& requested,
        PRPhaseIdentity phaseIdentity,
        PRPhaseRootIdentity rootIdentity) const;

    // Atomic publication helper: committed is assigned only after the full
    // homotopy, continuity, and request-relative closure contract succeeds.
    PRPhaseDisplacementStatus tryCommit(
        const PRPhaseRhoERequest& parent,
        const PRPhaseRhoERequest& requested,
        PRPhaseIdentity phaseIdentity,
        PRPhaseRootIdentity rootIdentity,
        PRPhaseRhoEResult& committed,
        PRPhaseDisplacementDiagnostics& diagnostics) const;
};

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
