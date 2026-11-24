#ifndef PRPhaseCaloric_H
#define PRPhaseCaloric_H

#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "PRMixtureKernel.H"
#include "PengRobinsonEOS.H"

namespace Foam
{
namespace phaseEquilibrium
{

enum class PRPhaseCaloricStatus
{
    SUCCESS,
    INVALID_CONFIGURATION,
    MISSING_CALORIC_DATA,
    MOLECULAR_WEIGHT_MISMATCH,
    MISSING_BINARY_INTERACTION,
    INVALID_PRESSURE,
    INVALID_TEMPERATURE,
    INVALID_COMPOSITION,
    INVALID_PHASE_ROOT,
    UNSUPPORTED_CALORIC_RANGE,
    PHASE_ROOT_FAILURE,
    NONFINITE_PROPERTY,
    EOS_CLOSURE_FAILURE
};

const char* phaseCaloricStatusName(PRPhaseCaloricStatus status);

struct PRPhaseCaloricProperties
{
    PhaseRoot root = PhaseRoot::STABLE;
    std::vector<double> moleFractions;
    std::vector<double> massFractions;
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double molecularWeight = std::numeric_limits<double>::quiet_NaN();
    double compressibility = std::numeric_limits<double>::quiet_NaN();
    double rho = std::numeric_limits<double>::quiet_NaN();
    double h = std::numeric_limits<double>::quiet_NaN();
    double e = std::numeric_limits<double>::quiet_NaN();
    double s = std::numeric_limits<double>::quiet_NaN();
};

struct PRPhaseCaloricResult
{
    PRPhaseCaloricStatus status =
        PRPhaseCaloricStatus::INVALID_CONFIGURATION;
    PhaseEquilibriumStatus phaseStatus =
        PhaseEquilibriumStatus::INVALID_CONFIGURATION;
    std::optional<PRPhaseCaloricProperties> state;
    std::string limitingSpecies;
    double supportedTemperatureMinimum =
        std::numeric_limits<double>::quiet_NaN();
    double supportedTemperatureMaximum =
        std::numeric_limits<double>::quiet_NaN();
    double pressureRelativeResidual =
        std::numeric_limits<double>::infinity();
    double energyIdentityRelativeResidual =
        std::numeric_limits<double>::infinity();

    bool success() const
    {
        return status == PRPhaseCaloricStatus::SUCCESS && state.has_value();
    }
};

struct PRPhaseCaloricTemperatureRange
{
    PRPhaseCaloricStatus status =
        PRPhaseCaloricStatus::INVALID_CONFIGURATION;
    PhaseEquilibriumStatus phaseStatus =
        PhaseEquilibriumStatus::INVALID_CONFIGURATION;
    double minimum = std::numeric_limits<double>::quiet_NaN();
    double maximum = std::numeric_limits<double>::quiet_NaN();
    std::string limitingSpecies;

    bool success() const
    {
        return status == PRPhaseCaloricStatus::SUCCESS
            && phaseStatus == PhaseEquilibriumStatus::SUCCESS
            && std::isfinite(minimum) && std::isfinite(maximum)
            && minimum <= maximum;
    }
};

// Low-level caloric adapter for an explicitly selected PR liquid or vapor
// root. Inputs and phase-equilibrium work use a mole basis; returned CFD
// properties use a mass basis and the existing NASA7 sensible-energy
// reference. Entropy uses the NASA standard state plus ideal mixing/pressure
// and PR departure terms on the same phase composition. This class is a
// phase-property tool, not a stability test or a TP/PH/rho-E flash.
class PRPhaseCaloric
{
private:
    struct TemperatureRange
    {
        double minimum = std::numeric_limits<double>::quiet_NaN();
        double maximum = std::numeric_limits<double>::quiet_NaN();
        double nasaMinimum = std::numeric_limits<double>::quiet_NaN();
        double nasaMaximum = std::numeric_limits<double>::quiet_NaN();
        double lowTemperatureOffset = 0.0;
        double lowTemperatureEntropyOffset = 0.0;
        double handoffMinimum = std::numeric_limits<double>::quiet_NaN();
        double handoffMaximum = std::numeric_limits<double>::quiet_NaN();
        double handoffCpCorrection = 0.0;
        bool hasLowTemperatureTable = false;
        bool hasC1Handoff = false;
    };

    SpeciesCriticalDB criticalData_;
    BinaryInteractionDB kij_;
    std::vector<std::string> speciesNames_;
    PRMixtureKernel kernel_;
    PengRobinsonEOS eos_;
    std::vector<TemperatureRange> caloricRanges_;
    PRPhaseCaloricStatus configurationStatus_;
    PhaseEquilibriumStatus configurationPhaseStatus_;
    std::string configurationSpecies_;

public:
    PRPhaseCaloric(
        const SpeciesCriticalDB& criticalData,
        const BinaryInteractionDB& kij,
        const std::vector<std::string>& speciesNames,
        double compositionTolerance = 1.0e-10);

    bool configured() const
    {
        return configurationStatus_ == PRPhaseCaloricStatus::SUCCESS;
    }

    PRPhaseCaloricStatus configurationStatus() const
    {
        return configurationStatus_;
    }

    PhaseEquilibriumStatus configurationPhaseStatus() const
    {
        return configurationPhaseStatus_;
    }

    const std::string& configurationSpecies() const
    {
        return configurationSpecies_;
    }

    const std::vector<std::string>& speciesNames() const
    {
        return speciesNames_;
    }

    // The basis is irrelevant here: only the active-species set determines
    // the common caloric interval. Fractions must still be normalized.
    PRPhaseCaloricTemperatureRange supportedTemperatureRange(
        const std::vector<double>& fractions) const;

    // Returns the ideal-mixture Cp used by evaluate() on the exact same
    // low-temperature table, C1 handoff correction, and NASA bridge.
    std::optional<double> idealCpMass(
        double temperature,
        const std::vector<double>& massFractions) const;

    PRPhaseCaloricResult evaluate(
        double pressure,
        double temperature,
        const std::vector<double>& moleFractions,
        PhaseRoot root) const;
};

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
