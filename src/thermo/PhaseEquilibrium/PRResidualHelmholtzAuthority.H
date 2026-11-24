#ifndef PRResidualHelmholtzAuthority_H
#define PRResidualHelmholtzAuthority_H

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "BinaryInteractionDB.H"
#include "SpeciesCriticalDB.H"

namespace Foam
{
namespace phaseEquilibrium
{

enum class PRResidualHelmholtzStatus
{
    SUCCESS,
    CAPABILITY_DISABLED,
    INVALID_CONFIGURATION,
    EMPTY_SPECIES_SET,
    DUPLICATE_SPECIES,
    UNKNOWN_SPECIES,
    INVALID_CRITICAL_PROPERTIES,
    MISSING_BINARY_INTERACTION,
    INVALID_BINARY_INTERACTION,
    ALPHA_ROOT_INTERVAL_EXIT,
    INVALID_TEMPERATURE,
    OUTSIDE_CERTIFIED_TEMPERATURE_INTERVAL,
    INVALID_DENSITY,
    INVALID_COMPOSITION_SIZE,
    NONFINITE_COMPOSITION,
    NEGATIVE_COMPOSITION,
    COMPOSITION_NOT_NORMALIZED,
    INVALID_IDEAL_CALORIC_STATE,
    INVALID_MOLAR_VOLUME,
    NONFINITE_POTENTIAL,
    NONPOSITIVE_PRESSURE,
    NONPOSITIVE_MOLAR_CV,
    MECHANICAL_INSTABILITY,
    NONPOSITIVE_MOLAR_CP,
    NONPOSITIVE_SOUND_SPEED_SQUARED
};

const char* residualHelmholtzStatusName(PRResidualHelmholtzStatus status);

struct PRIdealMolarCaloricState
{
    // The caller owns the ideal-gas reference and composition terms. These
    // values must all describe the same (T, v, x) state used by the request.
    double internalEnergy = std::numeric_limits<double>::quiet_NaN();
    double entropy = std::numeric_limits<double>::quiet_NaN();
    double cv = std::numeric_limits<double>::quiet_NaN();
};

struct PRResidualHelmholtzRequest
{
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double density = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> moleFractions;
    PRIdealMolarCaloricState ideal;
};

struct PRResidualHelmholtzToken
{
    double residualHelmholtz = 0.0;
    double dResidualHelmholtzDT = 0.0;
    double d2ResidualHelmholtzDT2 = 0.0;
    double dResidualHelmholtzDv = 0.0;
    double d2ResidualHelmholtzDv2 = 0.0;
    double d2ResidualHelmholtzDTDv = 0.0;

    double attraction = 0.0;
    double dAttractionDT = 0.0;
    double d2AttractionDT2 = 0.0;
    double coVolume = 0.0;

    double fade = 0.0;
    double dFadeDT = 0.0;
    double d2FadeDT2 = 0.0;
    double d3FadeDT3 = 0.0;
    bool alphaEvaluated = false;
};

struct PRResidualHelmholtzState
{
    std::vector<double> moleFractions;
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double density = std::numeric_limits<double>::quiet_NaN();
    double molecularWeight = std::numeric_limits<double>::quiet_NaN();
    double molarVolume = std::numeric_limits<double>::quiet_NaN();

    double pressure = std::numeric_limits<double>::quiet_NaN();
    double internalEnergy = std::numeric_limits<double>::quiet_NaN();
    double enthalpy = std::numeric_limits<double>::quiet_NaN();
    double entropy = std::numeric_limits<double>::quiet_NaN();
    double cv = std::numeric_limits<double>::quiet_NaN();
    double cp = std::numeric_limits<double>::quiet_NaN();

    double pressureTemperatureDerivative =
        std::numeric_limits<double>::quiet_NaN();
    double pressureMolarVolumeDerivative =
        std::numeric_limits<double>::quiet_NaN();
    double isothermalPressureDensityDerivative =
        std::numeric_limits<double>::quiet_NaN();
    double soundSpeedSquared = std::numeric_limits<double>::quiet_NaN();

    PRResidualHelmholtzToken token;
};

struct PRResidualHelmholtzResult
{
    PRResidualHelmholtzStatus status =
        PRResidualHelmholtzStatus::INVALID_CONFIGURATION;
    std::optional<PRResidualHelmholtzState> state;
    std::string limitingSpecies;

    bool success() const
    {
        return status == PRResidualHelmholtzStatus::SUCCESS
            && state.has_value();
    }
};

// Candidate high-temperature residual-potential owner. It freezes critical
// properties and every unlike-pair kij at construction, derives all requested
// properties from one residual Helmholtz token, and never repairs caller data.
// It is intentionally disconnected from runtime routing and production entry
// remains compile-time disabled pending physical selection and interval proof.
class PRResidualHelmholtzAuthority
{
private:
    struct RootDerivatives
    {
        double value = std::numeric_limits<double>::quiet_NaN();
        double first = std::numeric_limits<double>::quiet_NaN();
        double second = std::numeric_limits<double>::quiet_NaN();
    };

    struct FadeDerivatives
    {
        double value = 0.0;
        double first = 0.0;
        double second = 0.0;
        double third = 0.0;
    };

    std::vector<std::string> speciesNames_;
    std::vector<CriticalProperties> properties_;
    std::vector<double> binaryInteraction_;
    PRResidualHelmholtzStatus configurationStatus_;
    std::string limitingSpecies_;
    double minimumCertifiedAlphaRoot_;

    static bool finitePositive(double value);
    static double kappa(const CriticalProperties& properties);
    static RootDerivatives alphaRoot(
        const CriticalProperties& properties,
        double temperature);
    static FadeDerivatives fade(double temperature);

    PRResidualHelmholtzStatus validateComposition(
        const std::vector<double>& moleFractions) const;

    PRResidualHelmholtzResult evaluateImpl(
        const PRResidualHelmholtzRequest& request) const;

public:
    static constexpr double universalGasConstant = 8.314462618;
    static constexpr double handoffTemperatureLower = 1000.0;
    static constexpr double handoffTemperatureUpper = 1325.0;
    static constexpr double certifiedTemperatureMaximum = 6000.0;
    static constexpr double compositionTolerance = 1.0e-10;
    static constexpr double minimumAlphaRootMargin = 0.02;
    static constexpr double minimumMolarCv = 1.0e-8;
    static constexpr double minimumSoundSpeedSquared = 1.0e-8;
    static constexpr bool productionCapability = false;

    PRResidualHelmholtzAuthority(
        const SpeciesCriticalDB& criticalData,
        const BinaryInteractionDB& binaryInteraction,
        const std::vector<std::string>& speciesNames);

    bool configured() const
    {
        return configurationStatus_ == PRResidualHelmholtzStatus::SUCCESS;
    }

    PRResidualHelmholtzStatus configurationStatus() const
    {
        return configurationStatus_;
    }

    const std::string& limitingSpecies() const
    {
        return limitingSpecies_;
    }

    double minimumCertifiedAlphaRoot() const
    {
        return minimumCertifiedAlphaRoot_;
    }

    const std::vector<std::string>& speciesNames() const
    {
        return speciesNames_;
    }

    // Default/runtime-facing entry. It cannot execute while the capability is
    // false, even for an otherwise valid diagnostic request.
    PRResidualHelmholtzResult evaluate(
        const PRResidualHelmholtzRequest& request) const;

    // Explicit nonproduction entry used only for isolated evidence. Callers
    // must not publish this state into RFAM/HEM/Fluent runtime owners.
    PRResidualHelmholtzResult evaluateDiagnostic(
        const PRResidualHelmholtzRequest& request) const;
};

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
