#ifndef PRMixtureKernel_H
#define PRMixtureKernel_H

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "BinaryInteractionDB.H"
#include "SpeciesCriticalDB.H"

namespace Foam
{
namespace phaseEquilibrium
{

enum class PhaseEquilibriumStatus
{
    SUCCESS,
    INVALID_CONFIGURATION,
    UNKNOWN_SPECIES,
    INVALID_CRITICAL_PROPERTIES,
    MISSING_BINARY_INTERACTION,
    INVALID_PRESSURE,
    INVALID_TEMPERATURE,
    INVALID_COMPOSITION_SIZE,
    NONFINITE_COMPOSITION,
    NEGATIVE_COMPOSITION,
    COMPOSITION_NOT_NORMALIZED,
    NONFINITE_PARAMETER,
    NO_PHYSICAL_ROOT,
    INVALID_LOG_ARGUMENT,
    WILSON_K_FAILURE,
    RACHFORD_RICE_NO_BRACKET,
    RACHFORD_RICE_NOT_CONVERGED,
    STABILITY_NOT_CONVERGED,
    FLASH_NOT_CONVERGED,
    PURE_SATURATION_NO_BRACKET,
    PURE_SATURATION_NOT_CONVERGED
};

const char* statusName(PhaseEquilibriumStatus status);

enum class PhaseRoot
{
    LIQUID,
    VAPOR,
    STABLE
};

enum class EquilibriumRegion
{
    UNKNOWN,
    SINGLE_LIQUID,
    SINGLE_VAPOR,
    SINGLE_SUPERCRITICAL,
    TWO_PHASE,
    PURE_COMPONENT_COEXISTENCE
};

struct PRMixtureParameters
{
    PhaseEquilibriumStatus status =
        PhaseEquilibriumStatus::INVALID_CONFIGURATION;
    std::vector<double> moleFractions;
    std::vector<double> speciesA;
    std::vector<double> speciesB;
    std::vector<double> pairA;
    std::vector<double> sumPairA;
    double a = std::numeric_limits<double>::quiet_NaN();
    double b = std::numeric_limits<double>::quiet_NaN();
    double A = std::numeric_limits<double>::quiet_NaN();
    double B = std::numeric_limits<double>::quiet_NaN();

    bool valid() const
    {
        return status == PhaseEquilibriumStatus::SUCCESS;
    }
};

struct PRRootSet
{
    PhaseEquilibriumStatus status =
        PhaseEquilibriumStatus::INVALID_CONFIGURATION;
    PRMixtureParameters mixture;
    std::vector<double> physicalRoots;

    bool valid() const
    {
        return status == PhaseEquilibriumStatus::SUCCESS;
    }
};

struct PRPhaseState
{
    PhaseEquilibriumStatus status =
        PhaseEquilibriumStatus::INVALID_CONFIGURATION;
    PhaseRoot requestedRoot = PhaseRoot::STABLE;
    std::vector<double> moleFractions;
    std::vector<double> lnFugacityCoefficients;
    double compressibility = std::numeric_limits<double>::quiet_NaN();
    double residualGibbsOverRT = std::numeric_limits<double>::quiet_NaN();
    std::size_t physicalRootCount = 0;

    bool valid() const
    {
        return status == PhaseEquilibriumStatus::SUCCESS;
    }
};

// Shared mole-basis Peng-Robinson mixture kernel for phase-equilibrium work.
// Construction freezes critical properties plus kij values/availability, so no
// external database lifetime or later mutation can change a configured flash.
// It never repairs, clips, or silently renormalizes an invalid caller state.
class PRMixtureKernel
{
private:
    std::vector<std::string> speciesNames_;
    std::vector<CriticalProperties> properties_;
    std::vector<double> binaryInteraction_;
    std::vector<bool> hasBinaryInteraction_;
    PhaseEquilibriumStatus configurationStatus_;
    double compositionTolerance_;

    PRPhaseState evaluateAtRoot(
        const PRMixtureParameters& mixture,
        double compressibility,
        PhaseRoot requestedRoot,
        std::size_t rootCount) const;

public:
    static constexpr double universalGasConstant = 8.314462618;

    PRMixtureKernel(
        const SpeciesCriticalDB& criticalData,
        const BinaryInteractionDB& kij,
        const std::vector<std::string>& speciesNames,
        double compositionTolerance = 1.0e-10);

    PhaseEquilibriumStatus configurationStatus() const
    {
        return configurationStatus_;
    }

    bool configured() const
    {
        return configurationStatus_ == PhaseEquilibriumStatus::SUCCESS;
    }

    const std::vector<std::string>& speciesNames() const
    {
        return speciesNames_;
    }

    const std::vector<CriticalProperties>& properties() const
    {
        return properties_;
    }

    double compositionTolerance() const
    {
        return compositionTolerance_;
    }

    PhaseEquilibriumStatus validateMoleFractions(
        const std::vector<double>& moleFractions) const;

    PhaseEquilibriumStatus massToMoleFractions(
        const std::vector<double>& massFractions,
        std::vector<double>& moleFractions) const;

    PhaseEquilibriumStatus moleToMassFractions(
        const std::vector<double>& moleFractions,
        std::vector<double>& massFractions) const;

    double mixtureMolecularWeight(
        const std::vector<double>& moleFractions) const;

    PhaseEquilibriumStatus wilsonKValues(
        double pressure,
        double temperature,
        std::vector<double>& K) const;

    PRMixtureParameters mixtureParameters(
        double pressure,
        double temperature,
        const std::vector<double>& moleFractions) const;

    PRRootSet roots(
        double pressure,
        double temperature,
        const std::vector<double>& moleFractions) const;

    PRPhaseState phaseState(
        double pressure,
        double temperature,
        const std::vector<double>& moleFractions,
        PhaseRoot root) const;
};

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
