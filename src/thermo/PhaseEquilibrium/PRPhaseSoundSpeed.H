#ifndef PRPhaseSoundSpeed_H
#define PRPhaseSoundSpeed_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "PRRhoEFlash.H"

namespace Foam
{
namespace phaseEquilibrium
{

enum class PhaseSoundSpeedStatus
{
    SUCCESS,
    INVALID_CONFIGURATION,
    INVALID_INPUT,
    BASE_STATE_INVALID,
    BASE_STATE_MISMATCH,
    NONFINITE_ENTROPY,
    NO_ENTROPY_BRACKET,
    FLASH_FAILURE,
    NOT_CONVERGED,
    PHASE_TOPOLOGY_CHANGE,
    INSUFFICIENT_FLASH_ACCURACY,
    PERTURBATION_NOT_CONVERGED,
    NONPOSITIVE_DERIVATIVE
};

const char* phaseSoundSpeedStatusName(PhaseSoundSpeedStatus status);

struct PhaseSoundSpeedOptions
{
    double relativeDensityPerturbation = 1.0e-4;
    double minimumRelativeDensityPerturbation = 1.0e-7;
    std::size_t entropyScanIntervals = 48;
    std::size_t maximumIterations = 100;
    std::size_t maximumTopologyRefinements = 12;
    std::size_t minimumDerivativeConvergenceComparisons = 2;
    double energyCorrectionSpanFactor = 8.0;
    double entropyAbsoluteTolerance = 1.0e-6;
    double energyIntervalRelativeTolerance = 1.0e-10;
    double maximumFlashDensityErrorFraction = 1.0e-3;
    double derivativeRelativeTolerance = 2.0e-4;
};

struct PhaseSoundSpeedResult
{
    PhaseSoundSpeedStatus status =
        PhaseSoundSpeedStatus::INVALID_CONFIGURATION;
    RhoEFlashStatus baseFlashStatus = RhoEFlashStatus::INVALID_CONFIGURATION;
    RhoEFlashStatus lowerFlashStatus = RhoEFlashStatus::INVALID_CONFIGURATION;
    RhoEFlashStatus upperFlashStatus = RhoEFlashStatus::INVALID_CONFIGURATION;
    EquilibriumRegion baseRegion = EquilibriumRegion::UNKNOWN;
    EquilibriumRegion lowerRegion = EquilibriumRegion::UNKNOWN;
    EquilibriumRegion upperRegion = EquilibriumRegion::UNKNOWN;

    double targetDensity = std::numeric_limits<double>::quiet_NaN();
    double targetEnergy = std::numeric_limits<double>::quiet_NaN();
    double targetEntropy = std::numeric_limits<double>::quiet_NaN();
    double lowerDensity = std::numeric_limits<double>::quiet_NaN();
    double upperDensity = std::numeric_limits<double>::quiet_NaN();
    double lowerPressure = std::numeric_limits<double>::quiet_NaN();
    double upperPressure = std::numeric_limits<double>::quiet_NaN();
    double lowerEntropy = std::numeric_limits<double>::quiet_NaN();
    double upperEntropy = std::numeric_limits<double>::quiet_NaN();
    double densityDerivative = std::numeric_limits<double>::quiet_NaN();
    double soundSpeedSquared = std::numeric_limits<double>::quiet_NaN();
    double soundSpeed = std::numeric_limits<double>::quiet_NaN();
    double requestedRelativeDensityPerturbation =
        std::numeric_limits<double>::quiet_NaN();
    double usedRelativeDensityPerturbation =
        std::numeric_limits<double>::quiet_NaN();
    double lowerEntropyResidual = std::numeric_limits<double>::infinity();
    double upperEntropyResidual = std::numeric_limits<double>::infinity();
    std::size_t lowerIterations = 0;
    std::size_t upperIterations = 0;
    std::size_t lowerValidSamples = 0;
    std::size_t upperValidSamples = 0;
    std::size_t lowerInvalidSamples = 0;
    std::size_t upperInvalidSamples = 0;
    std::size_t stencilRefinements = 0;
    std::size_t topologyRefinements = 0;
    std::size_t derivativeConvergenceComparisons = 0;
    double derivativeRelativeChange = std::numeric_limits<double>::infinity();
    std::optional<RhoEFlashResult> baseState;
    std::optional<RhoEFlashResult> lowerState;
    std::optional<RhoEFlashResult> upperState;

    bool success() const
    {
        return status == PhaseSoundSpeedStatus::SUCCESS
            && baseState.has_value() && lowerState.has_value()
            && upperState.has_value()
            && std::isfinite(soundSpeed) && soundSpeed > 0.0;
    }
};

// Equilibrium HEM acoustic closure.  It computes (dp/drho)_s by solving two
// nearby rho-E flashes at the same overall composition and target entropy.
// This is deliberately independent of frozen-phase/Wood approximations; an
// invalid phase interval or an unresolved entropy bracket is a hard failure.
class PRPhaseSoundSpeed
{
private:
    struct DensitySolve
    {
        PhaseSoundSpeedStatus status =
            PhaseSoundSpeedStatus::NO_ENTROPY_BRACKET;
        RhoEFlashResult state;
        double entropyResidual = std::numeric_limits<double>::infinity();
        std::size_t iterations = 0;
        std::size_t validSamples = 0;
        std::size_t invalidSamples = 0;
    };

    PRRhoEFlash rhoEFlash_;
    PhaseSoundSpeedOptions options_;
    PhaseSoundSpeedStatus configurationStatus_;

    bool entropyConverged(double value, double target) const;
    DensitySolve solveAtDensity(
        double targetDensity,
        double targetEntropy,
        double energySeed,
        double pressureSeed,
        double temperatureSeed,
        EquilibriumRegion requiredRegion,
        const std::vector<double>& massFractions) const;
    PhaseSoundSpeedResult evaluateVerified(
        const RhoEFlashResult& base,
        const std::vector<double>& overallMassFractions) const;

public:
    PRPhaseSoundSpeed(
        const SpeciesCriticalDB& criticalData,
        const BinaryInteractionDB& kij,
        const std::vector<std::string>& speciesNames,
        const RhoEFlashOptions& rhoEOptions,
        const PhaseSoundSpeedOptions& options = PhaseSoundSpeedOptions());

    PhaseSoundSpeedStatus configurationStatus() const
    {
        return configurationStatus_;
    }

    const PhaseSoundSpeedOptions& options() const
    {
        return options_;
    }

    const PRRhoEFlash& rhoEFlash() const
    {
        return rhoEFlash_;
    }

    PhaseSoundSpeedResult evaluate(
        const RhoEFlashResult& base,
        const std::vector<double>& overallMassFractions) const;

    PhaseSoundSpeedResult evaluateRhoE(
        double targetDensity,
        double targetEnergy,
        const std::vector<double>& overallMassFractions) const;
};

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
