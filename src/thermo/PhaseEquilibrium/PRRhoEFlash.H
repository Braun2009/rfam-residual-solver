#ifndef PRRhoEFlash_H
#define PRRhoEFlash_H

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "ConservedPhaseMixture.H"
#include "PRPhaseCaloric.H"
#include "PRTPFlash.H"

namespace Foam
{
namespace phaseEquilibrium
{

enum class RhoEFlashStatus
{
    SUCCESS,
    INVALID_CONFIGURATION,
    INVALID_OPTIONS,
    INVALID_DENSITY,
    NONFINITE_TARGET_ENERGY,
    INVALID_COMPOSITION,
    OUT_OF_CONFIGURED_RANGE,
    REQUIRED_REGION_MISMATCH,
    TP_FLASH_FAILURE,
    TP_CLOSURE_FAILURE,
    CALORIC_FAILURE,
    UNSUPPORTED_CALORIC_RANGE,
    CONSERVED_MIXTURE_FAILURE,
    PURE_SATURATION_NO_BRACKET,
    NO_VALID_SATURATION_INTERVAL,
    NO_RHOE_BRACKET,
    NOT_CONVERGED,
    RESIDUAL_FAILURE
};

const char* rhoEFlashStatusName(RhoEFlashStatus status);

struct RhoEFlashOptions
{
    // A physical search box is mandatory. NaN defaults intentionally fail
    // closed instead of inventing an EOS domain.
    double minimumPressure = std::numeric_limits<double>::quiet_NaN();
    double maximumPressure = std::numeric_limits<double>::quiet_NaN();
    double minimumTemperature = std::numeric_limits<double>::quiet_NaN();
    double maximumTemperature = std::numeric_limits<double>::quiet_NaN();
    double initialPressure = std::numeric_limits<double>::quiet_NaN();
    double initialTemperature = std::numeric_limits<double>::quiet_NaN();

    std::size_t pressureScanIntervals = 64;
    std::size_t maximumIterations = 100;
    double finiteDifferenceStep = 1.0e-4;
    // The finite-difference perturbation is not a root-bracketing tolerance.
    // Keep a separate pressure interval criterion so a pure saturation rho-E
    // solve cannot report a state with a visible energy residual merely because
    // the Jacobian step is comparatively large.
    double saturationLogPressureTolerance = 1.0e-10;
    double logStepLimit = 1.0;
    double densityAbsoluteTolerance = 1.0e-8;
    double densityRelativeTolerance = 1.0e-8;
    double energyAbsoluteTolerance = 1.0e-5;
    double energyRelativeTolerance = 1.0e-10;
    double tpRootRelativeTolerance = 1.0e-10;
    TPFlashOptions tpOptions;
    ConservedMixtureTolerances conservedTolerances;
};

struct RhoEFlashSeed
{
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double temperature = std::numeric_limits<double>::quiet_NaN();
    // A strict local continuation is opt-in. When false, the canonical cold
    // search is used so a caller-provided seed cannot select a different root.
    bool localOnly = false;
    EquilibriumRegion requiredRegion = EquilibriumRegion::UNKNOWN;
};

struct RhoEFlashResult
{
    RhoEFlashStatus status = RhoEFlashStatus::INVALID_CONFIGURATION;
    PhaseEquilibriumStatus phaseStatus =
        PhaseEquilibriumStatus::INVALID_CONFIGURATION;
    PRPhaseCaloricStatus caloricStatus =
        PRPhaseCaloricStatus::INVALID_CONFIGURATION;
    ConservedMixtureStatus conservedStatus =
        ConservedMixtureStatus::INVALID_INPUT;

    double targetDensity = std::numeric_limits<double>::quiet_NaN();
    double targetEnergy = std::numeric_limits<double>::quiet_NaN();
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double density = std::numeric_limits<double>::quiet_NaN();
    double energy = std::numeric_limits<double>::quiet_NaN();
    double entropy = std::numeric_limits<double>::quiet_NaN();
    double densityResidual = std::numeric_limits<double>::infinity();
    double densityRelativeResidual = std::numeric_limits<double>::infinity();
    double energyResidual = std::numeric_limits<double>::infinity();
    double energyRelativeResidual = std::numeric_limits<double>::infinity();
    double tpFugacityResidual = std::numeric_limits<double>::infinity();
    double tpMaterialResidual = std::numeric_limits<double>::infinity();
    double phaseEnergyIdentityResidual =
        std::numeric_limits<double>::infinity();
    double mixtureMaterialResidual = std::numeric_limits<double>::infinity();
    double mixtureEnergyIdentityResidual =
        std::numeric_limits<double>::infinity();
    double vaporMoleFraction = std::numeric_limits<double>::quiet_NaN();

    bool pureVolumeEnergyLeverRule = false;
    std::size_t iterations = 0;
    std::size_t validSamples = 0;
    std::size_t invalidSamples = 0;

    TPFlashResult tpFlash;
    std::optional<PRPhaseCaloricProperties> liquidPhase;
    std::optional<PRPhaseCaloricProperties> vaporPhase;
    std::optional<PRPhaseCaloricProperties> singlePhase;
    std::optional<ConservedPhaseMixture> conservedMixture;

    bool success() const
    {
        if (status != RhoEFlashStatus::SUCCESS || !tpFlash.valid())
            return false;
        if (singlePhase.has_value())
            return !liquidPhase.has_value() && !vaporPhase.has_value()
                && !conservedMixture.has_value();
        return liquidPhase.has_value() && vaporPhase.has_value()
            && conservedMixture.has_value();
    }
};

// Strict C2P inversion for the phase-aware PR foundation. This class is not a
// production CFD hook yet: callers must still pass the model-capability gate.
class PRRhoEFlash
{
private:
    struct Evaluation
    {
        bool valid = false;
        RhoEFlashStatus status = RhoEFlashStatus::TP_FLASH_FAILURE;
        PhaseEquilibriumStatus phaseStatus =
            PhaseEquilibriumStatus::INVALID_CONFIGURATION;
        PRPhaseCaloricStatus caloricStatus =
            PRPhaseCaloricStatus::INVALID_CONFIGURATION;
        ConservedMixtureStatus conservedStatus =
            ConservedMixtureStatus::INVALID_INPUT;
        double pressure = std::numeric_limits<double>::quiet_NaN();
        double temperature = std::numeric_limits<double>::quiet_NaN();
        double density = std::numeric_limits<double>::quiet_NaN();
        double energy = std::numeric_limits<double>::quiet_NaN();
        double entropy = std::numeric_limits<double>::quiet_NaN();
        double densityResidual = std::numeric_limits<double>::infinity();
        double densityRelativeResidual =
            std::numeric_limits<double>::infinity();
        double energyResidual = std::numeric_limits<double>::infinity();
        double energyRelativeResidual =
            std::numeric_limits<double>::infinity();
        double tpFugacityResidual = std::numeric_limits<double>::infinity();
        double tpMaterialResidual = std::numeric_limits<double>::infinity();
        double phaseEnergyIdentityResidual =
            std::numeric_limits<double>::infinity();
        double mixtureMaterialResidual =
            std::numeric_limits<double>::infinity();
        double mixtureEnergyIdentityResidual =
            std::numeric_limits<double>::infinity();
        double vaporMoleFraction =
            std::numeric_limits<double>::quiet_NaN();
        bool pureVolumeEnergyLeverRule = false;
        TPFlashResult tp;
        std::optional<PRPhaseCaloricProperties> liquid;
        std::optional<PRPhaseCaloricProperties> vapor;
        std::optional<PRPhaseCaloricProperties> single;
        std::optional<ConservedPhaseMixture> mixture;
    };

    std::vector<std::string> speciesNames_;
    PRTPFlash tpFlash_;
    PRPhaseCaloric caloric_;
    RhoEFlashOptions options_;
    RhoEFlashStatus configurationStatus_;

    Evaluation assemblePT(
        double pressure,
        double temperature,
        const std::vector<double>& massFractions) const;

    Evaluation evaluatePT(
        double pressure,
        double temperature,
        double targetDensity,
        double targetEnergy,
        const std::vector<double>& massFractions) const;

    Evaluation evaluatePureSaturationPressure(
        double pressure,
        double targetDensity,
        double targetEnergy,
        const std::vector<double>& massFractions) const;

    RhoEFlashResult commit(
        double targetDensity,
        double targetEnergy,
        const Evaluation& evaluation,
        std::size_t iterations,
        std::size_t validSamples,
        std::size_t invalidSamples) const;
    std::optional<RhoEFlashResult> refinePureSaturation(
        double targetDensity,
        double targetEnergy,
        const std::vector<double>& massFractions,
        Evaluation& best,
        std::size_t& iterationsPerformed,
        std::size_t& validSamples,
        std::size_t& invalidSamples) const;
    std::optional<RhoEFlashResult> solveGeneralSeed(
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
        std::size_t& invalidSamples) const;

public:
    PRRhoEFlash(
        const SpeciesCriticalDB& criticalData,
        const BinaryInteractionDB& kij,
        const std::vector<std::string>& speciesNames,
        const RhoEFlashOptions& options);

    RhoEFlashStatus configurationStatus() const
    {
        return configurationStatus_;
    }

    const RhoEFlashOptions& options() const
    {
        return options_;
    }

    const PRTPFlash& tpFlash() const
    {
        return tpFlash_;
    }

    RhoEFlashResult flashMassFractions(
        double targetDensity,
        double targetEnergy,
        const std::vector<double>& overallMassFractions) const;

    RhoEFlashResult flashMassFractions(
        double targetDensity,
        double targetEnergy,
        const std::vector<double>& overallMassFractions,
        const RhoEFlashSeed& seed) const;

    // Construct a complete equilibrium state at prescribed p,T,Y. A pure
    // component saturation point remains fail-closed because p,T,Y does not
    // determine its phase fraction.
    RhoEFlashResult initializePTYMassFractions(
        double pressure,
        double temperature,
        const std::vector<double>& overallMassFractions) const;
};

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
