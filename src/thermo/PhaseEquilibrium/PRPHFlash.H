#ifndef PRPHFlash_H
#define PRPHFlash_H

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

enum class PHFlashStatus
{
    SUCCESS,
    INVALID_CONFIGURATION,
    INVALID_OPTIONS,
    INVALID_PRESSURE,
    NONFINITE_TARGET_ENTHALPY,
    INVALID_COMPOSITION,
    INVALID_TEMPERATURE_INTERVAL,
    TP_FLASH_FAILURE,
    TP_CLOSURE_FAILURE,
    UNSUPPORTED_CALORIC_RANGE,
    CALORIC_FAILURE,
    NONPOSITIVE_LATENT_ENTHALPY,
    CONSERVED_MIXTURE_FAILURE,
    NO_VALID_CALORIC_INTERVAL,
    NO_ENTHALPY_BRACKET,
    NOT_CONVERGED,
    ENTHALPY_RESIDUAL_FAILURE
};

const char* phFlashStatusName(PHFlashStatus status);

struct PHFlashOptions
{
    // The first PH implementation deliberately requires an explicit search
    // interval. NaN defaults make an omitted physical domain fail closed.
    double minimumTemperature =
        std::numeric_limits<double>::quiet_NaN();
    double maximumTemperature =
        std::numeric_limits<double>::quiet_NaN();
    std::size_t scanIntervals = 64;
    std::size_t maximumIterations = 200;
    double temperatureAbsoluteTolerance = 1.0e-10;
    double temperatureRelativeTolerance = 1.0e-12;
    double enthalpyAbsoluteTolerance = 1.0e-6;
    double enthalpyRelativeTolerance = 1.0e-10;
    double rootRelativeTolerance = 1.0e-10;
    ConservedMixtureTolerances conservedTolerances;
};

struct PHFlashResult
{
    PHFlashStatus status = PHFlashStatus::INVALID_CONFIGURATION;
    PhaseEquilibriumStatus phaseStatus =
        PhaseEquilibriumStatus::INVALID_CONFIGURATION;
    PRPhaseCaloricStatus caloricStatus =
        PRPhaseCaloricStatus::INVALID_CONFIGURATION;
    ConservedMixtureStatus conservedStatus =
        ConservedMixtureStatus::INVALID_INPUT;

    double pressure = std::numeric_limits<double>::quiet_NaN();
    double targetEnthalpy = std::numeric_limits<double>::quiet_NaN();
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double vaporMoleFraction = std::numeric_limits<double>::quiet_NaN();
    double computedEnthalpy = std::numeric_limits<double>::quiet_NaN();
    double enthalpyResidual = std::numeric_limits<double>::infinity();
    double enthalpyRelativeResidual =
        std::numeric_limits<double>::infinity();
    double entropy = std::numeric_limits<double>::quiet_NaN();
    double tpFugacityResidual = std::numeric_limits<double>::infinity();
    double tpMaterialResidual = std::numeric_limits<double>::infinity();
    double phaseEnergyIdentityResidual =
        std::numeric_limits<double>::infinity();
    double mixtureMaterialResidual =
        std::numeric_limits<double>::infinity();
    double mixtureEnergyIdentityResidual =
        std::numeric_limits<double>::infinity();

    bool pureComponentEnthalpyLeverRule = false;
    std::size_t iterations = 0;
    std::size_t validTemperatureSamples = 0;
    std::size_t invalidTemperatureSamples = 0;

    // tpFlash preserves the TP problem's semantics. In particular, a pure
    // coexistence TP payload keeps beta=NaN; vaporMoleFraction above is the
    // additional fraction determined by the PH enthalpy constraint.
    TPFlashResult tpFlash;
    std::optional<PRPhaseCaloricProperties> liquidPhase;
    std::optional<PRPhaseCaloricProperties> vaporPhase;
    std::optional<PRPhaseCaloricProperties> singlePhase;
    std::optional<ConservedPhaseMixture> conservedMixture;

    bool success() const
    {
        if (status != PHFlashStatus::SUCCESS || !tpFlash.valid())
            return false;
        if (singlePhase.has_value())
            return !liquidPhase.has_value() && !vaporPhase.has_value()
                && !conservedMixture.has_value();
        return liquidPhase.has_value() && vaporPhase.has_value()
            && conservedMixture.has_value();
    }
};

// Pressure-enthalpy equilibrium inversion. TP equilibrium remains on a mole
// basis, while all exposed phase and conserved properties are on a mass basis.
// The ascending scan returns the first strictly converged valid bracket; this
// is not a proof that a non-monotone h(T) surface has one globally unique root.
class PRPHFlash
{
private:
    struct TemperatureEvaluation
    {
        bool valid = false;
        PHFlashStatus status = PHFlashStatus::TP_FLASH_FAILURE;
        PhaseEquilibriumStatus phaseStatus =
            PhaseEquilibriumStatus::INVALID_CONFIGURATION;
        PRPhaseCaloricStatus caloricStatus =
            PRPhaseCaloricStatus::INVALID_CONFIGURATION;
        ConservedMixtureStatus conservedStatus =
            ConservedMixtureStatus::INVALID_INPUT;
        double temperature = std::numeric_limits<double>::quiet_NaN();
        double vaporMoleFraction =
            std::numeric_limits<double>::quiet_NaN();
        double enthalpy = std::numeric_limits<double>::quiet_NaN();
        double entropy = std::numeric_limits<double>::quiet_NaN();
        double residual = std::numeric_limits<double>::infinity();
        double tpFugacityResidual =
            std::numeric_limits<double>::infinity();
        double tpMaterialResidual =
            std::numeric_limits<double>::infinity();
        double phaseEnergyIdentityResidual =
            std::numeric_limits<double>::infinity();
        double mixtureMaterialResidual =
            std::numeric_limits<double>::infinity();
        double mixtureEnergyIdentityResidual =
            std::numeric_limits<double>::infinity();
        bool pureLeverRule = false;
        TPFlashResult tp;
        std::optional<PRPhaseCaloricProperties> liquid;
        std::optional<PRPhaseCaloricProperties> vapor;
        std::optional<PRPhaseCaloricProperties> single;
        std::optional<ConservedPhaseMixture> mixture;
    };

    std::vector<std::string> speciesNames_;
    PRTPFlash tpFlash_;
    PRPhaseCaloric caloric_;
    PHFlashOptions options_;
    PHFlashStatus configurationStatus_;

    TemperatureEvaluation evaluateTemperature(
        double pressure,
        double targetEnthalpy,
        const std::vector<double>& massFractions,
        double temperature) const;

    TemperatureEvaluation evaluatePureCoexistence(
        double pressure,
        double targetEnthalpy,
        const std::vector<double>& massFractions,
        const PureSaturationResult& saturation) const;

    PHFlashResult commit(
        double pressure,
        double targetEnthalpy,
        const TemperatureEvaluation& evaluation,
        std::size_t iterations,
        std::size_t validSamples,
        std::size_t invalidSamples) const;

public:
    PRPHFlash(
        const SpeciesCriticalDB& criticalData,
        const BinaryInteractionDB& kij,
        const std::vector<std::string>& speciesNames,
        const PHFlashOptions& options,
        const TPFlashOptions& tpOptions = TPFlashOptions());

    PHFlashStatus configurationStatus() const
    {
        return configurationStatus_;
    }

    const PHFlashOptions& options() const
    {
        return options_;
    }

    const PRTPFlash& tpFlash() const
    {
        return tpFlash_;
    }

    PHFlashResult flashMassFractions(
        double pressure,
        double targetEnthalpy,
        const std::vector<double>& overallMassFractions) const;
};

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
