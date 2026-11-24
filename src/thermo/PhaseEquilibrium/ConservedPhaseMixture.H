#ifndef ConservedPhaseMixture_H
#define ConservedPhaseMixture_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace Foam
{
namespace phaseEquilibrium
{

enum class ConservedMixtureStatus
{
    SUCCESS,
    INVALID_INPUT,
    INVALID_PHASE_COMPOSITION,
    MATERIAL_BALANCE_FAILED,
    PHASE_ENERGY_IDENTITY_FAILED,
    MIXTURE_ENERGY_IDENTITY_FAILED
};

struct ConservedMixtureTolerances
{
    double compositionAbsolute = 1.0e-10;
    double energyRelative = 1.0e-10;
};

// Intensive properties for one equilibrium phase at the common p and T.
// molecularWeight is kg/mol; rho is kg/m3; h/e are J/kg and s is J/(kg K).
// All must use the same species reference as the CFD conserved energy.
struct ConservedPhaseProperties
{
    double molecularWeight = std::numeric_limits<double>::quiet_NaN();
    double rho = std::numeric_limits<double>::quiet_NaN();
    double h = std::numeric_limits<double>::quiet_NaN();
    double e = std::numeric_limits<double>::quiet_NaN();
    double s = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> massFractions;
};

struct ConservedPhaseMixture
{
    double vaporMoleFraction = std::numeric_limits<double>::quiet_NaN();
    double vaporMassFraction = std::numeric_limits<double>::quiet_NaN();
    double vaporVolumeFraction = std::numeric_limits<double>::quiet_NaN();
    double rho = std::numeric_limits<double>::quiet_NaN();
    double h = std::numeric_limits<double>::quiet_NaN();
    double e = std::numeric_limits<double>::quiet_NaN();
    double s = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> overallMassFractions;
};

struct ConservedMixtureResult
{
    ConservedMixtureStatus status = ConservedMixtureStatus::INVALID_INPUT;
    std::optional<ConservedPhaseMixture> state;
    double materialResidual = std::numeric_limits<double>::infinity();
    double phaseEnergyIdentityResidual =
        std::numeric_limits<double>::infinity();
    double mixtureEnergyIdentityResidual =
        std::numeric_limits<double>::infinity();

    bool success() const
    {
        return status == ConservedMixtureStatus::SUCCESS && state.has_value();
    }
};

namespace detail
{

constexpr double maximumCompositionTolerance = 1.0e-6;
constexpr double maximumEnergyRelativeTolerance = 1.0e-6;

inline bool validTolerance(double value, double maximum)
{
    return std::isfinite(value) && value >= 0.0 && value <= maximum;
}

inline bool validOptionalEntropy(double value)
{
    return std::isfinite(value) || std::isnan(value);
}

inline bool validComposition(
    const std::vector<double>& composition,
    double tolerance)
{
    if (composition.empty()) return false;

    double sum = 0.0;
    for (double value : composition)
    {
        if (!std::isfinite(value) || value < 0.0 || value > 1.0)
            return false;
        sum += value;
    }
    return std::isfinite(sum) && std::abs(sum - 1.0) <= tolerance;
}

inline double relativeResidual(double lhs, double rhs)
{
    const double scale = std::max({std::abs(lhs), std::abs(rhs), 1.0});
    return std::abs(lhs - rhs) / scale;
}

inline double phaseEnergyResidual(
    double pressure,
    const ConservedPhaseProperties& phase)
{
    return relativeResidual(phase.e, phase.h - pressure / phase.rho);
}

} // End namespace detail

// Assemble mass-basis HEM properties from a molar TP-flash split. The expected
// overall composition is an independent input so this function checks, rather
// than assumes, the flash material balance. No output state is exposed unless
// every invariant passes.
inline ConservedMixtureResult assembleConservedPhaseMixture(
    double pressure,
    double vaporMoleFraction,
    const ConservedPhaseProperties& liquid,
    const ConservedPhaseProperties& vapor,
    const std::vector<double>& expectedOverallMassFractions,
    const ConservedMixtureTolerances& tolerances = {})
{
    ConservedMixtureResult result;

    if (!std::isfinite(pressure) || pressure <= 0.0 ||
        !std::isfinite(vaporMoleFraction) ||
        vaporMoleFraction < 0.0 || vaporMoleFraction > 1.0 ||
        !std::isfinite(liquid.molecularWeight) ||
        !std::isfinite(vapor.molecularWeight) ||
        liquid.molecularWeight <= 0.0 || vapor.molecularWeight <= 0.0 ||
        !std::isfinite(liquid.rho) || !std::isfinite(vapor.rho) ||
        liquid.rho <= 0.0 || vapor.rho <= 0.0 ||
        !std::isfinite(liquid.h) || !std::isfinite(liquid.e) ||
        !std::isfinite(vapor.h) || !std::isfinite(vapor.e) ||
        !detail::validOptionalEntropy(liquid.s) ||
        !detail::validOptionalEntropy(vapor.s) ||
        !detail::validTolerance(
            tolerances.compositionAbsolute,
            detail::maximumCompositionTolerance) ||
        !detail::validTolerance(
            tolerances.energyRelative,
            detail::maximumEnergyRelativeTolerance))
        return result;

    if (liquid.massFractions.size() != vapor.massFractions.size() ||
        liquid.massFractions.size() != expectedOverallMassFractions.size() ||
        !detail::validComposition(
            liquid.massFractions, tolerances.compositionAbsolute) ||
        !detail::validComposition(
            vapor.massFractions, tolerances.compositionAbsolute) ||
        !detail::validComposition(
            expectedOverallMassFractions, tolerances.compositionAbsolute))
    {
        result.status = ConservedMixtureStatus::INVALID_PHASE_COMPOSITION;
        return result;
    }

    result.phaseEnergyIdentityResidual = std::max(
        detail::phaseEnergyResidual(pressure, liquid),
        detail::phaseEnergyResidual(pressure, vapor));
    if (!std::isfinite(result.phaseEnergyIdentityResidual) ||
        result.phaseEnergyIdentityResidual > tolerances.energyRelative)
    {
        result.status =
            ConservedMixtureStatus::PHASE_ENERGY_IDENTITY_FAILED;
        return result;
    }

    const double liquidMoles = 1.0 - vaporMoleFraction;
    const double totalMass = liquidMoles * liquid.molecularWeight
        + vaporMoleFraction * vapor.molecularWeight;
    if (!std::isfinite(totalMass) || totalMass <= 0.0) return result;

    const double vaporMassFraction =
        vaporMoleFraction * vapor.molecularWeight / totalMass;
    const double liquidMassFraction = 1.0 - vaporMassFraction;
    const double specificVolume = liquidMassFraction / liquid.rho
        + vaporMassFraction / vapor.rho;
    if (!std::isfinite(vaporMassFraction) ||
        vaporMassFraction < 0.0 || vaporMassFraction > 1.0 ||
        !std::isfinite(specificVolume) || specificVolume <= 0.0)
        return result;

    ConservedPhaseMixture candidate;
    candidate.vaporMoleFraction = vaporMoleFraction;
    candidate.vaporMassFraction = vaporMassFraction;
    candidate.vaporVolumeFraction =
        (vaporMassFraction / vapor.rho) / specificVolume;
    candidate.rho = 1.0 / specificVolume;
    candidate.h = liquidMassFraction * liquid.h
        + vaporMassFraction * vapor.h;
    candidate.e = liquidMassFraction * liquid.e
        + vaporMassFraction * vapor.e;
    if (std::isfinite(liquid.s) && std::isfinite(vapor.s))
        candidate.s = liquidMassFraction * liquid.s
            + vaporMassFraction * vapor.s;
    candidate.overallMassFractions.resize(liquid.massFractions.size(), 0.0);

    result.materialResidual = 0.0;
    for (std::size_t species = 0;
         species < candidate.overallMassFractions.size(); ++species)
    {
        candidate.overallMassFractions[species] =
            liquidMassFraction * liquid.massFractions[species]
            + vaporMassFraction * vapor.massFractions[species];
        result.materialResidual = std::max(
            result.materialResidual,
            std::abs(candidate.overallMassFractions[species]
                - expectedOverallMassFractions[species]));
    }

    if (!std::isfinite(result.materialResidual) ||
        result.materialResidual > tolerances.compositionAbsolute)
    {
        result.status = ConservedMixtureStatus::MATERIAL_BALANCE_FAILED;
        return result;
    }

    result.mixtureEnergyIdentityResidual = detail::relativeResidual(
        candidate.e, candidate.h - pressure / candidate.rho);
    if (!std::isfinite(candidate.vaporVolumeFraction) ||
        candidate.vaporVolumeFraction < 0.0 ||
        candidate.vaporVolumeFraction > 1.0 ||
        !std::isfinite(candidate.rho) || !std::isfinite(candidate.h) ||
        !std::isfinite(candidate.e) ||
        ((std::isfinite(liquid.s) && std::isfinite(vapor.s))
            ? !std::isfinite(candidate.s) : !std::isnan(candidate.s)) ||
        !std::isfinite(result.mixtureEnergyIdentityResidual) ||
        result.mixtureEnergyIdentityResidual > tolerances.energyRelative)
    {
        result.status =
            ConservedMixtureStatus::MIXTURE_ENERGY_IDENTITY_FAILED;
        return result;
    }

    result.status = ConservedMixtureStatus::SUCCESS;
    result.state = std::move(candidate);
    return result;
}

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
