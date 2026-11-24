#ifndef PRTPFlash_H
#define PRTPFlash_H

#include <cstddef>
#include <limits>
#include <vector>

#include "PRMixtureKernel.H"

namespace Foam
{
namespace phaseEquilibrium
{

struct TPFlashOptions
{
    std::size_t maxFlashIterations = 250;
    std::size_t maxStabilityIterations = 150;
    std::size_t maxRachfordRiceIterations = 200;
    std::size_t saturationScanIntervals = 2048;
    std::size_t maxSaturationIterations = 200;
    double fugacityTolerance = 1.0e-9;
    double materialTolerance = 1.0e-10;
    double stabilityTolerance = 1.0e-8;
    double iterationTolerance = 1.0e-10;
    double saturationTemperatureTolerance = 1.0e-9;
    double successiveSubstitutionDamping = 0.5;
};

struct RachfordRiceResult
{
    PhaseEquilibriumStatus status =
        PhaseEquilibriumStatus::INVALID_CONFIGURATION;
    double vaporMoleFraction = std::numeric_limits<double>::quiet_NaN();
    double residual = std::numeric_limits<double>::quiet_NaN();
    std::size_t iterations = 0;

    bool valid() const
    {
        return status == PhaseEquilibriumStatus::SUCCESS;
    }

};

struct StabilityResult
{
    PhaseEquilibriumStatus status =
        PhaseEquilibriumStatus::INVALID_CONFIGURATION;
    bool unstable = false;
    double minimumTpd = std::numeric_limits<double>::quiet_NaN();
    double vaporLikeTpd = std::numeric_limits<double>::quiet_NaN();
    double liquidLikeTpd = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> vaporLikeComposition;
    std::vector<double> liquidLikeComposition;
    std::size_t vaporLikeIterations = 0;
    std::size_t liquidLikeIterations = 0;
    std::size_t attemptedTrials = 0;
    std::size_t convergedTrials = 0;
    std::size_t failedTrials = 0;

    bool valid() const
    {
        return status == PhaseEquilibriumStatus::SUCCESS;
    }

    bool allTrialsConverged() const
    {
        return valid() && attemptedTrials > 0
            && convergedTrials == attemptedTrials && failedTrials == 0;
    }
};

struct PureSaturationResult
{
    PhaseEquilibriumStatus status =
        PhaseEquilibriumStatus::INVALID_CONFIGURATION;
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double liquidCompressibility = std::numeric_limits<double>::quiet_NaN();
    double vaporCompressibility = std::numeric_limits<double>::quiet_NaN();
    double logFugacityResidual = std::numeric_limits<double>::quiet_NaN();
    std::size_t iterations = 0;

    bool valid() const
    {
        return status == PhaseEquilibriumStatus::SUCCESS;
    }
};

struct TPFlashResult
{
    PhaseEquilibriumStatus status =
        PhaseEquilibriumStatus::INVALID_CONFIGURATION;
    EquilibriumRegion region = EquilibriumRegion::UNKNOWN;
    bool phaseFractionDetermined = false;
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double vaporMoleFraction = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> feedMoleFractions;
    std::vector<double> feedMassFractions;
    std::vector<double> liquidMoleFractions;
    std::vector<double> vaporMoleFractions;
    std::vector<double> singlePhaseMoleFractions;
    std::vector<double> liquidMassFractions;
    std::vector<double> vaporMassFractions;
    std::vector<double> singlePhaseMassFractions;
    PRPhaseState liquidState;
    PRPhaseState vaporState;
    PRPhaseState singlePhaseState;
    double feedMinimumTpd = std::numeric_limits<double>::quiet_NaN();
    double liquidMinimumTpd = std::numeric_limits<double>::quiet_NaN();
    double vaporMinimumTpd = std::numeric_limits<double>::quiet_NaN();
    double maximumLogFugacityResidual =
        std::numeric_limits<double>::quiet_NaN();
    double maximumMoleMaterialResidual =
        std::numeric_limits<double>::quiet_NaN();
    double rachfordRiceResidual = std::numeric_limits<double>::quiet_NaN();
    std::size_t flashIterations = 0;

    bool valid() const
    {
        return status == PhaseEquilibriumStatus::SUCCESS;
    }

    // This confirms only a complete phase split/root payload. It does not
    // assert that caloric rho/h/e closure has been assembled or validated.
    bool phaseSplitReady() const
    {
        if (!valid()) return false;
        if (region == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE
            || region == EquilibriumRegion::UNKNOWN)
            return false;
        if (region == EquilibriumRegion::TWO_PHASE)
            return phaseFractionDetermined
                && liquidState.valid() && vaporState.valid();
        return singlePhaseState.valid()
            && !singlePhaseMoleFractions.empty()
            && !singlePhaseMassFractions.empty();
    }
};

class PRTPFlash
{
private:
    struct StabilityTrialResult
    {
        PhaseEquilibriumStatus status =
            PhaseEquilibriumStatus::INVALID_CONFIGURATION;
        double tpd = std::numeric_limits<double>::quiet_NaN();
        std::vector<double> composition;
        std::size_t iterations = 0;
    };

    PRMixtureKernel kernel_;
    TPFlashOptions options_;
    PhaseEquilibriumStatus optionsStatus_;

    StabilityTrialResult stabilityTrial(
        double pressure,
        double temperature,
        const std::vector<double>& feed,
        const PRPhaseState& feedState,
        const std::vector<double>& initialK,
        bool vaporLike) const;

    TPFlashResult flashPure(
        double pressure,
        double temperature,
        const std::vector<double>& feed) const;

    PhaseEquilibriumStatus phaseCompositions(
        const std::vector<double>& feed,
        const std::vector<double>& K,
        double vaporMoleFraction,
        std::vector<double>& liquid,
        std::vector<double>& vapor) const;

public:
    PRTPFlash(
        const SpeciesCriticalDB& criticalData,
        const BinaryInteractionDB& kij,
        const std::vector<std::string>& speciesNames,
        const TPFlashOptions& options = TPFlashOptions());

    const PRMixtureKernel& kernel() const
    {
        return kernel_;
    }

    const TPFlashOptions& options() const
    {
        return options_;
    }

    PhaseEquilibriumStatus configurationStatus() const;

    RachfordRiceResult solveRachfordRice(
        const std::vector<double>& feedMoleFractions,
        const std::vector<double>& K) const;

    StabilityResult analyzeStability(
        double pressure,
        double temperature,
        const std::vector<double>& feedMoleFractions) const;

    PureSaturationResult solvePureSaturationTemperature(
        double pressure,
        double minimumTemperature,
        double maximumTemperature) const;

    TPFlashResult flashMoleFractions(
        double pressure,
        double temperature,
        const std::vector<double>& feedMoleFractions) const;

    TPFlashResult flashMassFractions(
        double pressure,
        double temperature,
        const std::vector<double>& feedMassFractions) const;
};

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
