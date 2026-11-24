#include "EquilibriumThermoModel.H"

#include "RealFluidState.H"
#include "SpeciesRegistry.H"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <utility>

namespace Foam
{
namespace thermo
{
constexpr double R_UNIVERSAL = 8.314462618;
namespace
{

constexpr double compositionTolerance = 1.0e-8;
constexpr double closureTolerance = 1.0e-7;
constexpr double universalGasConstant = 8.314462618;

bool normalized(const std::vector<double>& values)
{
    if (values.empty()) return false;
    double sum = 0.0;
    for (double value : values)
    {
        if (!std::isfinite(value) || value < 0.0 || value > 1.0)
            return false;
        sum += value;
    }
    return std::isfinite(sum) && std::abs(sum - 1.0) <= 1.0e-8;
}

bool close(double lhs, double rhs, double tolerance = closureTolerance)
{
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) return false;
    return std::abs(lhs - rhs) <= tolerance
        *std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

bool sameComposition(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs)
{
    if (lhs.size() != rhs.size() || lhs.empty()) return false;
    for (std::size_t i = 0; i < lhs.size(); ++i)
        if (!close(lhs[i], rhs[i], compositionTolerance)) return false;
    return true;
}

bool entropyPayloadValid(bool available, double entropy)
{
    return available ? std::isfinite(entropy) : std::isnan(entropy);
}

bool thermodynamicPayloadClosed(
    const EquilibriumThermodynamicPayload& payload,
    double pressure,
    double temperature,
    double density,
    const EquilibriumAcousticProperties& acoustic)
{
    if (!payload.valid() || payload.closure != acoustic.closure
        || !std::isfinite(pressure) || !std::isfinite(temperature)
        || !std::isfinite(density) || density <= 0.0)
        return false;

    const double densitySquared = density*density;
    const double expectedEnergyDensityDerivative =
        (pressure - temperature*payload.dpdT_rho)/densitySquared;
    const double expectedCp = payload.Cv
        + temperature*payload.dpdT_rho*payload.dpdT_rho
            /(densitySquared*payload.dpdrho_T);
    const double expectedSoundSpeedSquared =
        (payload.Cp/payload.Cv)*payload.dpdrho_T;
    return close(payload.dedT_rho, payload.Cv)
        && close(payload.dedrho_T, expectedEnergyDensityDerivative)
        && close(payload.Cp, expectedCp)
        && close(
            acoustic.pressureDensityDerivative,
            expectedSoundSpeedSquared);
}

bool energyIdentityClosed(
    double internalEnergy,
    double enthalpy,
    double pressure,
    double density)
{
    if (!std::isfinite(internalEnergy) || !std::isfinite(enthalpy)
        || !std::isfinite(pressure) || !std::isfinite(density)
        || density <= 0.0)
        return false;
    const double pressureWork = pressure/density;
    const double residual = internalEnergy - (enthalpy - pressureWork);
    const double scale = std::max({
        std::abs(internalEnergy), std::abs(enthalpy),
        std::abs(pressureWork), 1.0});
    return std::abs(residual) <= closureTolerance*scale;
}

bool phaseValidAtPressure(
    const EquilibriumPhaseProperties& phase,
    double pressure,
    double temperature,
    std::size_t speciesCount)
{
    const bool pv = phase.valid()
        && phase.massFractions.size() == speciesCount
        && phase.moleFractions.size() == speciesCount
        && energyIdentityClosed(phase.e, phase.h, pressure, phase.rho)
        && close(
            phase.compressibility,
            pressure*phase.molecularWeight
                /(phase.rho*universalGasConstant*temperature));
    return pv;
}

bool fraction(double value)
{
    return std::isfinite(value) && value >= 0.0 && value <= 1.0;
}

bool knownAcousticClosure(AcousticClosure closure)
{
    switch (closure)
    {
        case AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE:
        case AcousticClosure::HOMOGENEOUS_EQUILIBRIUM:
            return true;
        case AcousticClosure::UNKNOWN:
            return false;
    }
    return false;
}

bool knownEquilibriumRegion(EquilibriumRegion region)
{
    switch (region)
    {
        case EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE:
        case EquilibriumRegion::SINGLE_LIQUID:
        case EquilibriumRegion::SINGLE_VAPOR:
        case EquilibriumRegion::SUPERCRITICAL:
        case EquilibriumRegion::TWO_PHASE:
        case EquilibriumRegion::CONDENSED_SOLID:
            return true;
        case EquilibriumRegion::UNKNOWN:
            return false;
    }
    return false;
}

bool molecularWeightsForNames(
    const std::vector<std::string>& names,
    std::vector<double>& molecularWeights)
{
    molecularWeights.clear();
    if (names.empty()) return false;

    SpeciesRegistry registry;
    molecularWeights.reserve(names.size());
    for (std::size_t i = 0; i < names.size(); ++i)
    {
        const std::string& name = names[i];
        for (std::size_t previous = 0; previous < i; ++previous)
        {
            if (names[previous] == name)
            {
                molecularWeights.clear();
                return false;
            }
        }
        if (!registry.has(name))
        {
            molecularWeights.clear();
            return false;
        }
        const double molecularWeight =
            registry.molecularWeight(registry.index(name));
        if (!std::isfinite(molecularWeight) || molecularWeight <= 0.0)
        {
            molecularWeights.clear();
            return false;
        }
        molecularWeights.push_back(molecularWeight);
    }
    return true;
}

std::vector<double> molecularWeightsForNames(
    const std::vector<std::string>& names)
{
    std::vector<double> molecularWeights;
    molecularWeightsForNames(names, molecularWeights);
    return molecularWeights;
}

const std::vector<std::string>& defaultLegacySpeciesNames()
{
    static const std::vector<std::string> names = {"O2", "C3H8"};
    return names;
}

const char* legacyStatusName(RealFluidState::Status status)
{
    switch (status)
    {
        case RealFluidState::Status::SUCCESS: return "SUCCESS";
        case RealFluidState::Status::INVERSION_FAILED: return "INVERSION_FAILED";
        case RealFluidState::Status::NEGATIVE_CV: return "NEGATIVE_CV";
        case RealFluidState::Status::NEGATIVE_CP: return "NEGATIVE_CP";
        case RealFluidState::Status::NEGATIVE_SOUNDSPEED:
            return "NEGATIVE_SOUNDSPEED";
        case RealFluidState::Status::NON_PHYSICAL: return "NON_PHYSICAL";
    }
    return "UNKNOWN";
}

EquilibriumThermodynamicPayload makeLegacyThermodynamicPayload(
    const RealFluidState& state)
{
    EquilibriumThermodynamicPayload payload;
    payload.closure = AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE;
    payload.Cp = state.Cp();
    payload.Cv = state.Cv();
    payload.dpdrho_T = state.dpdrho_T();
    payload.dpdT_rho = state.dpdT_rho();
    payload.dedrho_T = state.dedrho_T();
    payload.dedT_rho = state.dedT_rho();
    return payload;
}

bool massToMoleFractions(
    const std::vector<std::string>& names,
    const std::vector<double>& massFractions,
    std::vector<double>& moleFractions,
    double& molecularWeight)
{
    if (names.size() != massFractions.size() || !normalized(massFractions))
        return false;

    SpeciesRegistry registry;
    moleFractions.assign(massFractions.size(), 0.0);
    double totalMolesPerMass = 0.0;
    for (std::size_t i = 0; i < names.size(); ++i)
    {
        if (!registry.has(names[i])) return false;
        const double weight = registry.molecularWeight(registry.index(names[i]));
        if (!std::isfinite(weight) || weight <= 0.0) return false;
        moleFractions[i] = massFractions[i]/weight;
        totalMolesPerMass += moleFractions[i];
    }
    if (!std::isfinite(totalMolesPerMass) || totalMolesPerMass <= 0.0)
        return false;
    for (double& value : moleFractions) value /= totalMolesPerMass;
    molecularWeight = 1.0/totalMolesPerMass;
    return normalized(moleFractions) && std::isfinite(molecularWeight);
}

std::string makeLegacyModelId(const std::vector<std::string>& speciesNames)
{
    const std::string prefix = "legacy-real-fluid-single-phase-v3:";
    if (speciesNames.empty())
        return prefix + "constant-cp1000;species=implicit-o2-c3h8";

    std::string identity = prefix + "nasa7;species="
        + std::to_string(speciesNames.size()) + ';';
    for (const std::string& name : speciesNames)
        identity += std::to_string(name.size()) + ':' + name + ';';
    return identity;
}

} // End anonymous namespace

const char* equilibriumStatusName(EquilibriumStatus status)
{
    switch (status)
    {
        case EquilibriumStatus::SUCCESS: return "SUCCESS";
        case EquilibriumStatus::INVALID_CONFIGURATION:
            return "INVALID_CONFIGURATION";
        case EquilibriumStatus::INVALID_INPUT: return "INVALID_INPUT";
        case EquilibriumStatus::SPECIES_MISMATCH: return "SPECIES_MISMATCH";
        case EquilibriumStatus::PRIMITIVE_INITIALIZATION_UNSUPPORTED:
            return "PRIMITIVE_INITIALIZATION_UNSUPPORTED";
        case EquilibriumStatus::UNSUPPORTED_PROPERTY_RANGE:
            return "UNSUPPORTED_PROPERTY_RANGE";
        case EquilibriumStatus::UNDERDETERMINED_PHASE_FRACTION:
            return "UNDERDETERMINED_PHASE_FRACTION";
        case EquilibriumStatus::NO_EQUILIBRIUM_STATE:
            return "NO_EQUILIBRIUM_STATE";
        case EquilibriumStatus::NOT_BRACKETED: return "NOT_BRACKETED";
        case EquilibriumStatus::NOT_CONVERGED: return "NOT_CONVERGED";
        case EquilibriumStatus::ACOUSTIC_CLOSURE_FAILED:
            return "ACOUSTIC_CLOSURE_FAILED";
        case EquilibriumStatus::NONFINITE_PROPERTY:
            return "NONFINITE_PROPERTY";
        case EquilibriumStatus::RESIDUAL_CONTRACT_FAILED:
            return "RESIDUAL_CONTRACT_FAILED";
    }
    return "UNKNOWN";
}

const char* acousticClosureName(AcousticClosure closure)
{
    switch (closure)
    {
        case AcousticClosure::UNKNOWN: return "UNKNOWN";
        case AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE:
            return "FROZEN_COMPOSITION_SINGLE_PHASE";
        case AcousticClosure::HOMOGENEOUS_EQUILIBRIUM:
            return "HOMOGENEOUS_EQUILIBRIUM";
    }
    return "UNKNOWN";
}

const char* acousticEvidenceName(AcousticEvidence evidence)
{
    switch (evidence)
    {
        case AcousticEvidence::UNVERIFIED: return "UNVERIFIED";
        case AcousticEvidence::ANALYTIC_FROZEN_SINGLE_PHASE:
            return "ANALYTIC_FROZEN_SINGLE_PHASE";
        case AcousticEvidence::CONVERGED_RHOE_ISENTROPE:
            return "CONVERGED_RHOE_ISENTROPE";
        case AcousticEvidence::TABLE_INTERPOLATED:
            return "TABLE_INTERPOLATED";
    }
    return "UNVERIFIED";
}

bool EquilibriumAcousticProperties::valid() const
{
    if (!knownAcousticClosure(closure)
        || !std::isfinite(pressureDensityDerivative)
        || pressureDensityDerivative <= 0.0)
        return false;
    const double speed = std::sqrt(pressureDensityDerivative);
    return std::isfinite(speed) && speed > 0.0;
}

double EquilibriumAcousticProperties::soundSpeed() const
{
    return valid()
        ? std::sqrt(pressureDensityDerivative)
        : std::numeric_limits<double>::quiet_NaN();
}

bool EquilibriumThermodynamicPayload::valid() const
{
    return version == currentVersion && knownAcousticClosure(closure)
        && std::isfinite(Cp) && std::isfinite(Cv)
        && Cp > Cv && Cv > 0.0
        && std::isfinite(dpdrho_T) && dpdrho_T > 0.0
        && std::isfinite(dpdT_rho)
        && std::isfinite(dedrho_T)
        && std::isfinite(dedT_rho) && dedT_rho > 0.0;
}

const char* equilibriumRegionName(EquilibriumRegion region)
{
    switch (region)
    {
        case EquilibriumRegion::UNKNOWN: return "UNKNOWN";
        case EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE:
            return "UNCLASSIFIED_SINGLE_PHASE";
        case EquilibriumRegion::SINGLE_LIQUID: return "SINGLE_LIQUID";
        case EquilibriumRegion::SINGLE_VAPOR: return "SINGLE_VAPOR";
        case EquilibriumRegion::SUPERCRITICAL: return "SUPERCRITICAL";
        case EquilibriumRegion::TWO_PHASE: return "TWO_PHASE";
        case EquilibriumRegion::CONDENSED_SOLID: return "CONDENSED_SOLID";
    }
    return "UNKNOWN";
}

bool EquilibriumPhaseProperties::valid() const
{
    return normalized(massFractions) && normalized(moleFractions)
        && massFractions.size() == moleFractions.size()
        && std::isfinite(molecularWeight) && molecularWeight > 0.0
        && std::isfinite(rho) && rho > 0.0
        && std::isfinite(h) && std::isfinite(e)
        && std::isfinite(compressibility) && compressibility > 0.0
        && entropyPayloadValid(entropyAvailable, s);
}

bool EquilibriumState::valid() const
{
    const bool common = knownEquilibriumRegion(region)
        && normalized(overallMassFractions)
        && std::isfinite(p) && p > 0.0
        && std::isfinite(T) && T > 0.0
        && std::isfinite(rho) && rho > 0.0
        && std::isfinite(h) && std::isfinite(e)
        && acoustic.valid()
        && (!thermodynamics.has_value()
            || thermodynamicPayloadClosed(
                *thermodynamics, p, T, rho, acoustic))
        && entropyPayloadValid(entropyAvailable, s)
        && energyIdentityClosed(e, h, p, rho);
    if (!common) return false;

    if (region == EquilibriumRegion::TWO_PHASE)
    {
        if (acoustic.closure != AcousticClosure::HOMOGENEOUS_EQUILIBRIUM)
            return false;
        if (!liquid.has_value() || !vapor.has_value()
            || singlePhase.has_value()
            || !phaseValidAtPressure(
                *liquid, p, T, overallMassFractions.size())
            || !phaseValidAtPressure(
                *vapor, p, T, overallMassFractions.size())
            || !fraction(vaporMoleFraction)
            || !fraction(vaporMassFraction)
            || !fraction(vaporVolumeFraction))
            return false;

        const double totalMolarMass =
            (1.0 - vaporMoleFraction)*liquid->molecularWeight
            + vaporMoleFraction*vapor->molecularWeight;
        if (!std::isfinite(totalMolarMass) || totalMolarMass <= 0.0)
            return false;
        const double expectedMassFraction =
            vaporMoleFraction*vapor->molecularWeight/totalMolarMass;
        const double liquidMassFraction = 1.0 - expectedMassFraction;
        const double specificVolume = liquidMassFraction/liquid->rho
            + expectedMassFraction/vapor->rho;
        if (!std::isfinite(specificVolume) || specificVolume <= 0.0)
            return false;
        const double expectedVolumeFraction =
            (expectedMassFraction/vapor->rho)/specificVolume;

        if (!close(vaporMassFraction, expectedMassFraction, 1.0e-3)
            || !close(vaporVolumeFraction, expectedVolumeFraction, 1.0e-3)
            || !close(rho, 1.0/specificVolume, 1.0e-3)
            || !close(h, liquidMassFraction*liquid->h
                + vaporMassFraction*vapor->h, 1.0e-3)
            || !close(e, liquidMassFraction*liquid->e
                + vaporMassFraction*vapor->e, 1.0e-3))
            return false;

        for (std::size_t species = 0;
             species < overallMassFractions.size(); ++species)
        {
            const double expected =
                liquidMassFraction*liquid->massFractions[species]
                + vaporMassFraction*vapor->massFractions[species];
            if (!close(
                    overallMassFractions[species], expected,
                    compositionTolerance))
                return false;
        }

        const bool mixtureEntropyAvailable =
            liquid->entropyAvailable && vapor->entropyAvailable;
        return entropyAvailable == mixtureEntropyAvailable
            && (!entropyAvailable
                || close(s, liquidMassFraction*liquid->s
                    + vaporMassFraction*vapor->s));
    }

    if (region == EquilibriumRegion::SINGLE_LIQUID
        || region == EquilibriumRegion::SINGLE_VAPOR)
    {
        // Quality-clamp classification (2026-08-04): the Legacy HEM path
        // commits the saturated endpoint in the liquid/vapour payload with
        // the other phase absent; the PR path represents the same region
        // with the singlePhase payload.  Accept BOTH representations (the
        // PR path predates the clamp classification and uses
        // SINGLE_LIQUID/SINGLE_VAPOR for plain single-phase states).
        const bool wantLiquid = region == EquilibriumRegion::SINGLE_LIQUID;
        if (singlePhase.has_value()
            && !liquid.has_value() && !vapor.has_value())
        {
            // Plain PR single-phase representation: phase fractions are
            // unavailable, not clamped endpoint qualities. Keep the payload
            // canonical so stale HEM fractions cannot survive a region
            // reclassification and leak into downstream diagnostics.
            return std::isnan(vaporMoleFraction)
                && std::isnan(vaporMassFraction)
                && std::isnan(vaporVolumeFraction)
                && phaseValidAtPressure(
                    *singlePhase, p, T, overallMassFractions.size())
                && close(rho, singlePhase->rho)
                && close(h, singlePhase->h)
                && close(e, singlePhase->e)
                && sameComposition(
                    overallMassFractions, singlePhase->massFractions);
        }
        const EquilibriumPhaseProperties* ph = wantLiquid
            ? (liquid.has_value() ? &*liquid : nullptr)
            : (vapor.has_value() ? &*vapor : nullptr);
        if (singlePhase.has_value() || !ph
            || (wantLiquid ? vapor.has_value() : liquid.has_value())
            || !fraction(vaporMassFraction)
            || !fraction(vaporVolumeFraction)
            || (wantLiquid
                ? vaporMassFraction > 1e-6 || vaporVolumeFraction > 1e-6
                : vaporMassFraction < 1.0 - 1e-6
                    || vaporVolumeFraction < 1.0 - 1e-6)
            || !phaseValidAtPressure(
                *ph, p, T, overallMassFractions.size()))
            return false;
        // Quality-clamp classification (2026-08-04): q<=eps commits the
        // saturated endpoint as SINGLE_LIQUID, q>=1-eps as SINGLE_VAPOR.
        // The lever-rule rho carries the q->0/1 residual (e.g. q=6e-9
        // gives rho = rhoL*(1-2e-7)), so the endpoint closure uses a
        // 1e-3 diagnostic tolerance.
        return close(rho, ph->rho, 1.0e-3)
            && close(h, ph->h, 1.0e-3)
            && close(e, ph->e, 1.0e-3)
            && sameComposition(overallMassFractions, ph->massFractions);
    }

    if (region == EquilibriumRegion::CONDENSED_SOLID)
    {
        // Frozen gas + condensed solid: the gas carries the acoustic
        // closure; the condensed phase supplies only volume/mass/energy.
        if (!singlePhase.has_value() || !condensed.has_value()
            || liquid.has_value() || vapor.has_value()
            || !fraction(condensedMassFraction)
            || condensedMassFraction <= 0.0
            || condensedMassFraction >= 1.0
            || !phaseValidAtPressure(
                *singlePhase, p, T, overallMassFractions.size())
            || !phaseValidAtPressure(
                *condensed, p, T, overallMassFractions.size()))
            return false;
        const double gasFraction = 1.0 - condensedMassFraction;
        const double specificVolume = gasFraction/singlePhase->rho
            + condensedMassFraction/condensed->rho;
        if (!std::isfinite(specificVolume) || specificVolume <= 0.0)
            return false;
        if (!close(rho, 1.0/specificVolume)
            || !close(h, gasFraction*singlePhase->h
                + condensedMassFraction*condensed->h)
            || !close(e, gasFraction*singlePhase->e
                + condensedMassFraction*condensed->e))
        for (std::size_t species = 0;
             species < overallMassFractions.size(); ++species)
        {
            const double expected =
                gasFraction*singlePhase->massFractions[species]
                + condensedMassFraction*condensed->massFractions[species];
            if (!close(
                    overallMassFractions[species], expected,
                    compositionTolerance))
                return false;
        }
        return true;
    }

    if (!singlePhase.has_value() || liquid.has_value() || vapor.has_value()
        || !std::isnan(vaporMoleFraction)
        || !std::isnan(vaporMassFraction)
        || !std::isnan(vaporVolumeFraction)
        || !phaseValidAtPressure(
            *singlePhase, p, T, overallMassFractions.size()))
        return false;

    return sameComposition(
            overallMassFractions, singlePhase->massFractions)
        && close(rho, singlePhase->rho)
        && close(h, singlePhase->h)
        && close(e, singlePhase->e)
        && entropyAvailable == singlePhase->entropyAvailable
        && (!entropyAvailable || close(s, singlePhase->s));
}

bool validateEquilibriumSpeciesBasis(
    const EquilibriumState& state,
    const std::vector<double>& orderedSpeciesMolecularWeights,
    const EquilibriumSpeciesBasisTolerances& tolerances)
{
    if (!state.valid()
        || orderedSpeciesMolecularWeights.size()
            != state.overallMassFractions.size()
        || !std::isfinite(tolerances.compositionAbsolute)
        || tolerances.compositionAbsolute < 0.0
        || !std::isfinite(tolerances.molecularWeightRelative)
        || tolerances.molecularWeightRelative < 0.0)
        return false;

    for (double molecularWeight : orderedSpeciesMolecularWeights)
        if (!std::isfinite(molecularWeight) || molecularWeight <= 0.0)
            return false;

    const auto validPhaseBasis = [&](const EquilibriumPhaseProperties& phase)
    {
        if (phase.massFractions.size()
                != orderedSpeciesMolecularWeights.size()
            || phase.moleFractions.size()
                != orderedSpeciesMolecularWeights.size())
            return false;

        double totalMolesPerMass = 0.0;
        for (std::size_t species = 0;
             species < orderedSpeciesMolecularWeights.size(); ++species)
        {
            totalMolesPerMass += phase.massFractions[species]
                /orderedSpeciesMolecularWeights[species];
        }
        if (!std::isfinite(totalMolesPerMass) || totalMolesPerMass <= 0.0)
            return false;

        const double expectedMolecularWeight = 1.0/totalMolesPerMass;
        const double molecularWeightScale = std::max(
            std::abs(expectedMolecularWeight),
            std::abs(phase.molecularWeight));
        if (!std::isfinite(expectedMolecularWeight)
            || molecularWeightScale <= 0.0
            || std::abs(phase.molecularWeight - expectedMolecularWeight)
                > tolerances.molecularWeightRelative*molecularWeightScale)
            return false;

        for (std::size_t species = 0;
             species < orderedSpeciesMolecularWeights.size(); ++species)
        {
            const double expectedMoleFraction =
                (phase.massFractions[species]
                    /orderedSpeciesMolecularWeights[species])
                /totalMolesPerMass;
            if (!std::isfinite(expectedMoleFraction)
                || std::abs(
                    phase.moleFractions[species] - expectedMoleFraction)
                    > tolerances.compositionAbsolute)
                return false;
        }
        return true;
    };

    return (!state.liquid.has_value() || validPhaseBasis(*state.liquid))
        && (!state.vapor.has_value() || validPhaseBasis(*state.vapor))
        && (!state.singlePhase.has_value()
            || validPhaseBasis(*state.singlePhase))
        && (!state.condensed.has_value()
            || validPhaseBasis(*state.condensed));
}

bool EquilibriumThermoModel::stateValidForModel(
    const EquilibriumState& state) const
{
    std::vector<double> molecularWeights;
    return molecularWeightsForNames(speciesNames(), molecularWeights)
        && validateEquilibriumSpeciesBasis(state, molecularWeights);
}

EquilibriumResult EquilibriumThermoModel::initializePTY(
    double,
    double,
    const std::vector<double>&) const
{
    EquilibriumResult result;
    result.status = EquilibriumStatus::PRIMITIVE_INITIALIZATION_UNSUPPORTED;
    result.diagnostics.modelId = modelId();
    result.diagnostics.backendStatus =
        "PRIMITIVE_PTY_REQUIRES_COMPLETE_EQUILIBRIUM_STATE";
    return result;
}

LegacySinglePhaseThermoModel::LegacySinglePhaseThermoModel(
    std::vector<std::string> speciesNames)
:
    speciesNames_(std::move(speciesNames)),
    speciesMolecularWeights_(molecularWeightsForNames(
        speciesNames_.empty() ? defaultLegacySpeciesNames() : speciesNames_)),
    configurationValid_(
        speciesMolecularWeights_.size()
        == (speciesNames_.empty()
            ? defaultLegacySpeciesNames().size() : speciesNames_.size())),
    modelId_(makeLegacyModelId(speciesNames_))
{}

bool LegacySinglePhaseThermoModel::stateValidForModel(
    const EquilibriumState& state) const
{
    return configurationValid_
        && validateEquilibriumSpeciesBasis(
            state, speciesMolecularWeights_)
        && state.thermodynamics.has_value()
        && (state.thermodynamics->closure
                == AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE
            || (state.thermodynamics->closure
                    == AcousticClosure::HOMOGENEOUS_EQUILIBRIUM
                && (state.region == EquilibriumRegion::TWO_PHASE
                    || state.region == EquilibriumRegion::SINGLE_LIQUID
                    || state.region == EquilibriumRegion::SINGLE_VAPOR)));
}

EquilibriumResult LegacySinglePhaseThermoModel::recoverRhoE(
    double density,
    double internalEnergy,
    const std::vector<double>& overallMassFractions,
    const EquilibriumSeed* seed) const
{
#ifdef RFAM_THERMO_DIAGNOSTICS
    std::cerr << "[RECOV-ENTRY] rho=" << density
              << " e=" << internalEnergy << " Y0="
              << (overallMassFractions.empty()
                  ? -1.0 : overallMassFractions[0]) << "\n";
#endif
    EquilibriumResult result;
    result.diagnostics.modelId = modelId_;
    if (!configurationValid_)
    {
        result.status = EquilibriumStatus::INVALID_CONFIGURATION;
        result.diagnostics.backendStatus = "INVALID_CONFIGURATION";
        return result;
    }
    if (!std::isfinite(density) || density <= 0.0
        || !std::isfinite(internalEnergy)
        || !normalized(overallMassFractions))
    {
        result.status = EquilibriumStatus::INVALID_INPUT;
        result.diagnostics.backendStatus = "INVALID_INPUT";
        return result;
    }
    if ((!speciesNames_.empty()
            && speciesNames_.size() != overallMassFractions.size())
        || (speciesNames_.empty() && overallMassFractions.size() != 2))
    {
        result.status = EquilibriumStatus::SPECIES_MISMATCH;
        result.diagnostics.backendStatus = "SPECIES_MISMATCH";
        return result;
    }

    RealFluidState state;
    try
    {
        if (!speciesNames_.empty()) state.setSpecies(speciesNames_);
        if (gergTable_) state.setGERGTable(gergTable_);
        if (saturationTable_) state.setSaturationTable(saturationTable_);
    }
    catch (const std::exception& e)
    {
        // Never swallow the reason: a bare catch(...) here turned every setup
        // failure into INVALID_CONFIGURATION, naming neither the species nor
        // the table, and cost a full rebuild cycle to diagnose.
        result.status = EquilibriumStatus::INVALID_CONFIGURATION;
        result.diagnostics.backendStatus =
            std::string("INVALID_CONFIGURATION: ") + e.what();
        return result;
    }
    const double temperatureSeed = seed && std::isfinite(seed->temperature)
        ? seed->temperature : 300.0;
    const double phaseHintTemperature =
        seed && std::isfinite(seed->previousTemperature)
        ? seed->previousTemperature : temperatureSeed;
    if (seed
        && seed->previousRegion == EquilibriumRegion::TWO_PHASE
        && std::isfinite(seed->previousVaporMassFraction)
        && std::isfinite(phaseHintTemperature))
    {
        state.setTwoPhaseHint(
            phaseHintTemperature, seed->previousVaporMassFraction);
    }
    const RealFluidState::Status status = state.updateRhoEY(
        density, internalEnergy, overallMassFractions, temperatureSeed);
    result.diagnostics.backendStatus = legacyStatusName(status);
    if (status == RealFluidState::Status::SUCCESS
        && state.repairedComposition())
    {
        // 2026-08-07: conserved-field mass roundoff was clipped and
        // renormalised before the EOS closure; audited, not silent.
        result.diagnostics.backendStatus =
            std::string("REPAIRED_ROUNDOFF|") + result.diagnostics.backendStatus;
    }
    else if (status == RealFluidState::Status::SUCCESS
        && state.phaseUnstable())
    {
        // 2026-08-07 cold-region class 3: committed single-phase root on
        // a mechanically unstable branch; audited diagnostic.
        result.diagnostics.backendStatus =
            std::string("PHASE_UNSTABLE|") + result.diagnostics.backendStatus;
    }
    else if (status == RealFluidState::Status::INVERSION_FAILED
        && state.energyNotBracketed())
    {
        // 2026-08-07 cold-region class 2: the isochoric energy bracket
        // failed; report the cold-end gap and equivalent temperature
        // (mK => numerical boundary, 10s K => source-term problem).
        char buf[160];
        std::snprintf(buf, sizeof(buf),
            "ENERGY_NOT_BRACKETED:dT_eq=%.3fK", state.energyNotBracketedDTEq());
        result.diagnostics.backendStatus = buf;
    }
    if (status != RealFluidState::Status::SUCCESS)
    {
        result.status = EquilibriumStatus::NO_EQUILIBRIUM_STATE;
        return result;
    }

    const std::vector<std::string>& names = speciesNames_.empty()
        ? defaultLegacySpeciesNames() : speciesNames_;
    EquilibriumPhaseProperties phase;
    phase.massFractions = overallMassFractions;
    if (!massToMoleFractions(
            names, overallMassFractions, phase.moleFractions,
            phase.molecularWeight))
    {
        result.status = EquilibriumStatus::SPECIES_MISMATCH;
        result.diagnostics.backendStatus = "MOLE_BASIS_CONVERSION_FAILED";
        return result;
    }
    phase.rho = state.rho();
    phase.h = state.h();
    phase.e = state.e();
    phase.compressibility = state.Z();

    EquilibriumState candidate;
    if (state.isTwoPhase())
    {
        // Minimal HEM state (2026-08-04 audit): the saturated (Tsat, psat)
        // is the native p/T.  Liquid/vapour phases are the REAL saturation
        // endpoints from the saturation table (rhoL/rhoV, eL/eV at Tsat);
        // for binary mixtures the phase compositions are the coexisting
        // (x, y) pair (teqp v3 table), NOT the overall composition.
        // Quality clamps are classified: q<=eps -> SINGLE_LIQUID,
        // q>=1-eps -> SINGLE_VAPOR (so chemistry/boundary/acoustic routing
        // sees a single phase), otherwise TWO_PHASE.
        const double q = state.quality();
        // Composition basis: xO2 (liquid mole fraction) from the overall
        // mass fractions; species order is the model's (O2, C3H8) for the
        // 2-species diagnostic models — anything else fails closed.
        double xO2 = 1.0;
        if (overallMassFractions.size() >= 2)
        {
            // Phase split on the O2/C3H8 sub-composition only; other
            // species (H2O/CO2 LOX near-field dilution, residues after
            // chemistry) are passive diluents carried in the overall
            // mass fractions.  The HEM endpoints are taken from the
            // nearest v3 composition level (levelRows snapping).
            double nO2 = overallMassFractions[0] / 0.0319988;
            double nC3 = overallMassFractions[1] / 0.04409562;
            if ((nO2 + nC3) <= 0.0)
            {
                result.status = EquilibriumStatus::NONFINITE_PROPERTY;
                result.diagnostics.backendStatus =
                    "TWO_PHASE_UNSUPPORTED_SPECIES_BASIS";
                return result;
            }
            xO2 = nO2 / (nO2 + nC3);
        }
        double rhoV = 0.0, rhoL = 0.0, eV = 0.0, eL = 0.0;
        bool haveEndpoints = saturationTable_
            && saturationTable_->saturationEndpoints(
                state.T(), xO2, rhoV, rhoL, &eV, &eL);
        if (!haveEndpoints)
        {
            result.status = EquilibriumStatus::NONFINITE_PROPERTY;
            result.diagnostics.backendStatus = "TWO_PHASE_ENDPOINTS_UNAVAILABLE";
            return result;
        }
        // Coexisting vapour composition (teqp v3); NaN -> fall back to x
        // (pure fluids and uncovered mixture rows).
        // Near-azeotropic composition (2026-08-31): the minimal HEM does
        // NOT redistribute mass between phases; the 0.72-0.78 band sits at
        // the low-temperature O2/C3H8 azeotrope where y ~ x, so the vapor
        // phase keeps the liquid composition (xO2).  A composition-separated
        // flash would make the overall composition inconsistent with the
        // two-phase lever rule in EquilibriumState::valid().  The teqp
        // vaporComposition y is deliberately NOT used here.
        double yO2 = xO2;
        auto makePhase = [&](double compO2, double rho, double e)
            -> EquilibriumPhaseProperties
        {
            // Full 9-species phase composition: the O2/C3H8 sub-split is
            // (compO2, 1-compO2); passive diluents (H2O/CO2 residues)
            // keep their overall mass-fraction ratios within both phases.
            EquilibriumPhaseProperties ph;
            std::vector<double> y9(overallMassFractions.size(), 0.0);
            std::vector<double> x9(overallMassFractions.size(), 0.0);
            double diluentMass = 0.0;
            for (std::size_t si = 2; si < overallMassFractions.size(); ++si)
            {
                y9[si] = overallMassFractions[si];
                diluentMass += y9[si];
            }
            y9[0] = compO2 * (1.0 - diluentMass);
            y9[1] = (1.0 - compO2) * (1.0 - diluentMass);
            double molesPerMass = 0.0;
            for (std::size_t si = 0; si < y9.size(); ++si)
                molesPerMass += y9[si] / speciesMolecularWeights_[si];
            for (std::size_t si = 0; si < y9.size(); ++si)
                x9[si] = (y9[si] / speciesMolecularWeights_[si])
                    / molesPerMass;
            ph.massFractions = std::move(y9);
            ph.moleFractions = std::move(x9);
            ph.molecularWeight = 1.0 / molesPerMass;
            ph.rho = rho;
            ph.h = e + state.p() / rho;
            ph.e = e;
            ph.compressibility = state.p() * ph.molecularWeight
                / (rho * R_UNIVERSAL * state.T());
            return ph;
        };
        const EquilibriumPhaseProperties liquidPhase =
            makePhase(xO2, rhoL, eL);
        const EquilibriumPhaseProperties vaporPhase =
            makePhase(yO2, rhoV, eV);
        candidate.vaporMassFraction = q;
        candidate.vaporMoleFraction =
            q*liquidPhase.molecularWeight
            /(q*liquidPhase.molecularWeight
                + (1.0 - q)*vaporPhase.molecularWeight);
        candidate.liquid = liquidPhase;
        candidate.singlePhase.reset();
        // Acoustics are a DIAGNOSTIC PLACEHOLDER (self-consistent
        // ideal-gas-like closure, NOT a converged derivative); evidence is
        // UNVERIFIED and ModelCapabilities::phaseAwareSoundSpeed=false
        // keeps it out of production fluxes/CFL.  SINGLE_LIQUID/VAPOR use
        // the frozen-single-phase closure so single-phase consumers accept
        // the state (region carries the phase identity).
        candidate.acoustic.pressureDensityDerivative = state.a()*state.a();
        // NOTE (2026-08-04 audit): SINGLE_LIQUID/VAPOR keep the
        // HOMOGENEOUS_EQUILIBRIUM + UNVERIFIED acoustic closure so the
        // placeholder sound speed is NOT admitted to production fluxes
        // (AcousticStateAdmission rejects it); the region itself already
        // routes chemistry/boundaries to the correct single phase.
        candidate.acoustic.closure =
            AcousticClosure::HOMOGENEOUS_EQUILIBRIUM;
        candidate.acoustic.evidence =
            AcousticEvidence::UNVERIFIED;
        if (q <= 1e-8)
        {
            candidate.region = EquilibriumRegion::SINGLE_LIQUID;
            candidate.vaporVolumeFraction = 0.0;
            candidate.vapor.reset();
        }
        else if (q >= 1.0 - 1e-8)
        {
            candidate.region = EquilibriumRegion::SINGLE_VAPOR;
            candidate.vaporVolumeFraction = 1.0;
            candidate.vapor = vaporPhase;
            candidate.liquid.reset();
        }
        else
        {
            candidate.region = EquilibriumRegion::TWO_PHASE;
            const double specificVolume = (1.0 - q) / rhoL + q / rhoV;
            candidate.vaporVolumeFraction =
                (q / rhoV) / specificVolume;
            candidate.vapor = vaporPhase;
        }
    }
    else if (state.isCondensed())
    {
        // Frozen gas + condensed solid: singlePhase carries the gas
        // sub-composition; condensed carries the solid species with the
        // solid density.  Acoustic closure is the frozen gas (diagnostic
        // placeholder mixture estimate), consistent with the HEM route.
        const std::vector<double>& gasY = state.gasMassFractions();
        const std::vector<double>& solY =
            state.condensedMassFractions();
        const double cmf = state.condensedMassFraction();
        EquilibriumPhaseProperties gasPhase;
        gasPhase.massFractions = gasY;
        double molesPerMass = 0.0;
        for (std::size_t si = 0; si < gasY.size(); ++si)
            molesPerMass += gasY[si] / speciesMolecularWeights_[si];
        std::vector<double> x9(gasY.size(), 0.0);
        for (std::size_t si = 0; si < gasY.size(); ++si)
            x9[si] = (gasY[si] / speciesMolecularWeights_[si])
                / molesPerMass;
        gasPhase.moleFractions = std::move(x9);
        gasPhase.molecularWeight = 1.0 / molesPerMass;
        gasPhase.rho = state.gasRho();
        gasPhase.h = state.gasE() + state.p() / state.gasRho();
        gasPhase.e = state.gasE();
        gasPhase.compressibility = state.p() * gasPhase.molecularWeight
            / (state.gasRho() * R_UNIVERSAL * state.T());
        candidate.singlePhase = std::move(gasPhase);
        EquilibriumPhaseProperties solidPhase;
        // Phase-internal normalised composition (sum=1): the condensed
        // species renormalised to the solid fraction; valid() composes
        // overall = (1-cmf)*gasY + cmf*solidY.
        std::vector<double> solYnorm = solY;
        double solSum = 0.0;
        for (double sy : solYnorm) solSum += sy;
        if (solSum > 0.0)
            for (double& sy : solYnorm) sy /= solSum;
        solidPhase.massFractions = solYnorm;
        double sMolesPerMass = 0.0;
        for (std::size_t si = 0; si < solYnorm.size(); ++si)
            sMolesPerMass += solYnorm[si] / speciesMolecularWeights_[si];
        std::vector<double> sx9(solYnorm.size(), 0.0);
        for (std::size_t si = 0; si < solYnorm.size(); ++si)
            sx9[si] = (solYnorm[si] / speciesMolecularWeights_[si])
                / sMolesPerMass;
        solidPhase.moleFractions = std::move(sx9);
        solidPhase.molecularWeight = 1.0 / sMolesPerMass;
        const double solidRho = 1.0 / (
            (solYnorm.size() > 8 ? solYnorm[8] / 1562.0 : 0.0)
            + (solYnorm.size() > 6 ? solYnorm[6] / 917.0 : 0.0));
        solidPhase.rho = solidRho;
        // Per-kg-solid sensible energy: e_total = (1-cmf)*eGas + cmf*eSolid.
        solidPhase.e = (state.e() - state.gasE() * (1.0 - cmf)) / cmf;
        solidPhase.h = solidPhase.e + state.p() / solidRho;
        solidPhase.compressibility = state.p() * solidPhase.molecularWeight
            / (solidRho * R_UNIVERSAL * state.T());
        candidate.condensed = std::move(solidPhase);
        candidate.condensedMassFraction = cmf;
        candidate.region = EquilibriumRegion::CONDENSED_SOLID;
        candidate.vaporMassFraction = 1.0 - cmf;
        candidate.vaporMoleFraction = 1.0 - cmf;
        candidate.vaporVolumeFraction = state.gasRho() > 0.0
            ? ((1.0 - cmf) / state.gasRho())
                / ((1.0 - cmf) / state.gasRho() + cmf / solidRho)
            : std::numeric_limits<double>::quiet_NaN();
        candidate.acoustic.closure =
            AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE;
        candidate.acoustic.evidence =
            AcousticEvidence::ANALYTIC_FROZEN_SINGLE_PHASE;
        candidate.acoustic.pressureDensityDerivative = state.a()*state.a();
        candidate.liquid.reset();
        candidate.vapor.reset();
    }
    else
    {
        candidate.region = EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE;
        candidate.acoustic.closure =
            AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE;
        candidate.acoustic.evidence =
            AcousticEvidence::ANALYTIC_FROZEN_SINGLE_PHASE;
        candidate.acoustic.pressureDensityDerivative = state.a()*state.a();
        candidate.singlePhase = std::move(phase);
    }
    candidate.overallMassFractions = overallMassFractions;
    candidate.p = state.p();
    candidate.T = state.T();
    candidate.rho = state.rho();
    candidate.h = state.h();
    candidate.e = state.e();
    candidate.thermodynamics =
        makeLegacyThermodynamicPayload(state);
    if (candidate.thermodynamics)
        candidate.thermodynamics->closure = candidate.acoustic.closure;

    result.diagnostics.densityResidual = candidate.rho - density;
    result.diagnostics.energyResidual = candidate.e - internalEnergy;
    result.diagnostics.energyIdentityResidual =
        candidate.e - (candidate.h - candidate.p/candidate.rho);
    if (!stateValidForModel(candidate))
    {
        if(std::getenv("RFAM_VALID_DIAG"))
        std::cerr << "[VALID-DIAG] stateValidForModel failed: region="
                  << static_cast<int>(candidate.region)
                  << " closure="
                  << static_cast<int>(candidate.acoustic.closure)
                  << " sp=" << candidate.singlePhase.has_value()
                  << " cd=" << candidate.condensed.has_value()
                  << " liq=" << candidate.liquid.has_value()
                  << " vap=" << candidate.vapor.has_value()
                  << " cmf=" << candidate.condensedMassFraction
                  << " rho=" << candidate.rho
                  << " T=" << candidate.T
                  << " p=" << candidate.p
                  << " e=" << candidate.e
                  << " h=" << candidate.h
                  << " thermo=" << candidate.thermodynamics.has_value()
                  << " vmf=" << candidate.vaporMassFraction
                  << " cfg=" << configurationValid_
                  << " nY=" << candidate.overallMassFractions.size()
                  << " nMW=" << speciesMolecularWeights_.size()
                  << " Cp=" << (candidate.thermodynamics
                      ? candidate.thermodynamics->Cp : -1.0)
                  << " Cv=" << (candidate.thermodynamics
                      ? candidate.thermodynamics->Cv : -1.0)
                  << " dpdrho=" << (candidate.thermodynamics
                      ? candidate.thermodynamics->dpdrho_T : -1.0)
                  << " dpdT=" << (candidate.thermodynamics
                      ? candidate.thermodynamics->dpdT_rho : -1.0)
                  << " a2=" << candidate.acoustic.pressureDensityDerivative
                  << "\n";
        result.status = EquilibriumStatus::NONFINITE_PROPERTY;
        result.diagnostics.backendStatus = "INVALID_COMMITTED_STATE";
        return result;
    }

    result.status = EquilibriumStatus::SUCCESS;
    result.state = std::move(candidate);
    return result;
}

EquilibriumResult LegacySinglePhaseThermoModel::initializePTY(
    double pressure,
    double temperature,
    const std::vector<double>& overallMassFractions) const
{
    EquilibriumResult result;
    result.diagnostics.modelId = modelId_;
    if (!configurationValid_)
    {
        result.status = EquilibriumStatus::INVALID_CONFIGURATION;
        result.diagnostics.backendStatus = "INVALID_CONFIGURATION";
        return result;
    }
    if (!std::isfinite(pressure) || pressure <= 0.0
        || !std::isfinite(temperature) || temperature <= 0.0
        || !normalized(overallMassFractions))
    {
        result.status = EquilibriumStatus::INVALID_INPUT;
        result.diagnostics.backendStatus = "INVALID_INPUT";
        return result;
    }
    if ((!speciesNames_.empty()
            && speciesNames_.size() != overallMassFractions.size())
        || (speciesNames_.empty() && overallMassFractions.size() != 2))
    {
        result.status = EquilibriumStatus::SPECIES_MISMATCH;
        result.diagnostics.backendStatus = "SPECIES_MISMATCH";
        return result;
    }

    RealFluidState state;
    try
    {
        if (!speciesNames_.empty()) state.setSpecies(speciesNames_);
        if (gergTable_) state.setGERGTable(gergTable_);
        if (saturationTable_) state.setSaturationTable(saturationTable_);
    }
    catch (const std::exception& e)
    {
        // Never swallow the reason: a bare catch(...) here turned every setup
        // failure into INVALID_CONFIGURATION, naming neither the species nor
        // the table, and cost a full rebuild cycle to diagnose.
        result.status = EquilibriumStatus::INVALID_CONFIGURATION;
        result.diagnostics.backendStatus =
            std::string("INVALID_CONFIGURATION: ") + e.what();
        return result;
    }
    const RealFluidState::Status status =
        state.updatePTY(pressure, temperature, overallMassFractions);
    result.diagnostics.backendStatus = legacyStatusName(status);
    if (status != RealFluidState::Status::SUCCESS)
    {
        result.status = EquilibriumStatus::NO_EQUILIBRIUM_STATE;
        return result;
    }

    const std::vector<std::string>& names = speciesNames_.empty()
        ? defaultLegacySpeciesNames() : speciesNames_;
    EquilibriumPhaseProperties phase;
    phase.massFractions = overallMassFractions;
    if (!massToMoleFractions(
            names, overallMassFractions, phase.moleFractions,
            phase.molecularWeight))
    {
        result.status = EquilibriumStatus::SPECIES_MISMATCH;
        result.diagnostics.backendStatus = "MOLE_BASIS_CONVERSION_FAILED";
        return result;
    }
    phase.rho = state.rho();
    phase.h = state.h();
    phase.e = state.e();
    phase.compressibility = state.Z();

    EquilibriumState candidate;
    if (state.isTwoPhase())
    {
        // Minimal HEM state (2026-08-04 audit): the saturated (Tsat, psat)
        // is the native p/T.  Liquid/vapour phases are the REAL saturation
        // endpoints from the saturation table (rhoL/rhoV, eL/eV at Tsat);
        // for binary mixtures the phase compositions are the coexisting
        // (x, y) pair (teqp v3 table), NOT the overall composition.
        // Quality clamps are classified: q<=eps -> SINGLE_LIQUID,
        // q>=1-eps -> SINGLE_VAPOR (so chemistry/boundary/acoustic routing
        // sees a single phase), otherwise TWO_PHASE.
        const double q = state.quality();
        // Composition basis: xO2 (liquid mole fraction) from the overall
        // mass fractions; species order is the model's (O2, C3H8) for the
        // 2-species diagnostic models — anything else fails closed.
        double xO2 = 1.0;
        if (overallMassFractions.size() >= 2)
        {
            // Phase split on the O2/C3H8 sub-composition only; other
            // species (H2O/CO2 LOX near-field dilution, residues after
            // chemistry) are passive diluents carried in the overall
            // mass fractions.  The HEM endpoints are taken from the
            // nearest v3 composition level (levelRows snapping).
            double nO2 = overallMassFractions[0] / 0.0319988;
            double nC3 = overallMassFractions[1] / 0.04409562;
            if ((nO2 + nC3) <= 0.0)
            {
                result.status = EquilibriumStatus::NONFINITE_PROPERTY;
                result.diagnostics.backendStatus =
                    "TWO_PHASE_UNSUPPORTED_SPECIES_BASIS";
                return result;
            }
            xO2 = nO2 / (nO2 + nC3);
        }
        double rhoV = 0.0, rhoL = 0.0, eV = 0.0, eL = 0.0;
        bool haveEndpoints = saturationTable_
            && saturationTable_->saturationEndpoints(
                state.T(), xO2, rhoV, rhoL, &eV, &eL);
        if (!haveEndpoints)
        {
            result.status = EquilibriumStatus::NONFINITE_PROPERTY;
            result.diagnostics.backendStatus = "TWO_PHASE_ENDPOINTS_UNAVAILABLE";
            return result;
        }
        // Coexisting vapour composition (teqp v3); NaN -> fall back to x
        // (pure fluids and uncovered mixture rows).
        // Near-azeotropic composition (2026-08-31): the minimal HEM does
        // NOT redistribute mass between phases; the 0.72-0.78 band sits at
        // the low-temperature O2/C3H8 azeotrope where y ~ x, so the vapor
        // phase keeps the liquid composition (xO2).  A composition-separated
        // flash would make the overall composition inconsistent with the
        // two-phase lever rule in EquilibriumState::valid().  The teqp
        // vaporComposition y is deliberately NOT used here.
        double yO2 = xO2;
        auto makePhase = [&](double compO2, double rho, double e)
            -> EquilibriumPhaseProperties
        {
            // Full 9-species phase composition: the O2/C3H8 sub-split is
            // (compO2, 1-compO2); passive diluents (H2O/CO2 residues)
            // keep their overall mass-fraction ratios within both phases.
            EquilibriumPhaseProperties ph;
            std::vector<double> y9(overallMassFractions.size(), 0.0);
            std::vector<double> x9(overallMassFractions.size(), 0.0);
            double diluentMass = 0.0;
            for (std::size_t si = 2; si < overallMassFractions.size(); ++si)
            {
                y9[si] = overallMassFractions[si];
                diluentMass += y9[si];
            }
            y9[0] = compO2 * (1.0 - diluentMass);
            y9[1] = (1.0 - compO2) * (1.0 - diluentMass);
            double molesPerMass = 0.0;
            for (std::size_t si = 0; si < y9.size(); ++si)
                molesPerMass += y9[si] / speciesMolecularWeights_[si];
            for (std::size_t si = 0; si < y9.size(); ++si)
                x9[si] = (y9[si] / speciesMolecularWeights_[si])
                    / molesPerMass;
            ph.massFractions = std::move(y9);
            ph.moleFractions = std::move(x9);
            ph.molecularWeight = 1.0 / molesPerMass;
            ph.rho = rho;
            ph.h = e + state.p() / rho;
            ph.e = e;
            ph.compressibility = state.p() * ph.molecularWeight
                / (rho * R_UNIVERSAL * state.T());
            return ph;
        };
        const EquilibriumPhaseProperties liquidPhase =
            makePhase(xO2, rhoL, eL);
        const EquilibriumPhaseProperties vaporPhase =
            makePhase(yO2, rhoV, eV);
        candidate.vaporMassFraction = q;
        candidate.vaporMoleFraction =
            q*liquidPhase.molecularWeight
            /(q*liquidPhase.molecularWeight
                + (1.0 - q)*vaporPhase.molecularWeight);
        candidate.liquid = liquidPhase;
        candidate.singlePhase.reset();
        // Acoustics are a DIAGNOSTIC PLACEHOLDER (self-consistent
        // ideal-gas-like closure, NOT a converged derivative); evidence is
        // UNVERIFIED and ModelCapabilities::phaseAwareSoundSpeed=false
        // keeps it out of production fluxes/CFL.  SINGLE_LIQUID/VAPOR use
        // the frozen-single-phase closure so single-phase consumers accept
        // the state (region carries the phase identity).
        candidate.acoustic.pressureDensityDerivative = state.a()*state.a();
        // NOTE (2026-08-04 audit): SINGLE_LIQUID/VAPOR keep the
        // HOMOGENEOUS_EQUILIBRIUM + UNVERIFIED acoustic closure so the
        // placeholder sound speed is NOT admitted to production fluxes
        // (AcousticStateAdmission rejects it); the region itself already
        // routes chemistry/boundaries to the correct single phase.
        candidate.acoustic.closure =
            AcousticClosure::HOMOGENEOUS_EQUILIBRIUM;
        candidate.acoustic.evidence =
            AcousticEvidence::UNVERIFIED;
        if (q <= 1e-8)
        {
            candidate.region = EquilibriumRegion::SINGLE_LIQUID;
            candidate.vaporVolumeFraction = 0.0;
            candidate.vapor.reset();
        }
        else if (q >= 1.0 - 1e-8)
        {
            candidate.region = EquilibriumRegion::SINGLE_VAPOR;
            candidate.vaporVolumeFraction = 1.0;
            candidate.vapor = vaporPhase;
            candidate.liquid.reset();
        }
        else
        {
            candidate.region = EquilibriumRegion::TWO_PHASE;
            const double specificVolume = (1.0 - q) / rhoL + q / rhoV;
            candidate.vaporVolumeFraction =
                (q / rhoV) / specificVolume;
            candidate.vapor = vaporPhase;
        }
    }
    else
    {
        candidate.region = EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE;
        candidate.acoustic.closure =
            AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE;
        candidate.acoustic.evidence =
            AcousticEvidence::ANALYTIC_FROZEN_SINGLE_PHASE;
        candidate.acoustic.pressureDensityDerivative = state.a()*state.a();
        candidate.singlePhase = std::move(phase);
    }
    candidate.overallMassFractions = overallMassFractions;
    candidate.p = state.p();
    candidate.T = state.T();
    candidate.rho = state.rho();
    candidate.h = state.h();
    candidate.e = state.e();
    candidate.thermodynamics =
        makeLegacyThermodynamicPayload(state);
    if (candidate.thermodynamics)
        candidate.thermodynamics->closure = candidate.acoustic.closure;
    result.diagnostics.densityResidual = 0.0;
    result.diagnostics.energyResidual = 0.0;
    result.diagnostics.energyIdentityResidual =
        candidate.e - (candidate.h - candidate.p/candidate.rho);
    if (!stateValidForModel(candidate))
    {
        if(std::getenv("RFAM_VALID_DIAG"))
        std::cerr << "[VALID-DIAG] stateValidForModel failed: region="
                  << static_cast<int>(candidate.region)
                  << " rho=" << candidate.rho
                  << " T=" << candidate.T
                  << " p=" << candidate.p
                  << " e=" << candidate.e
                  << " h=" << candidate.h
                  << " thermo=" << candidate.thermodynamics.has_value()
                  << " vmf=" << candidate.vaporMassFraction
                  << " cfg=" << configurationValid_
                  << " nY=" << candidate.overallMassFractions.size()
                  << " nMW=" << speciesMolecularWeights_.size()
                  << " Cp=" << (candidate.thermodynamics
                      ? candidate.thermodynamics->Cp : -1.0)
                  << " Cv=" << (candidate.thermodynamics
                      ? candidate.thermodynamics->Cv : -1.0)
                  << " dpdrho=" << (candidate.thermodynamics
                      ? candidate.thermodynamics->dpdrho_T : -1.0)
                  << " dpdT=" << (candidate.thermodynamics
                      ? candidate.thermodynamics->dpdT_rho : -1.0)
                  << " a2=" << candidate.acoustic.pressureDensityDerivative
                  << "\n";
        result.status = EquilibriumStatus::NONFINITE_PROPERTY;
        result.diagnostics.backendStatus = "INVALID_COMMITTED_STATE";
        return result;
    }

    result.status = EquilibriumStatus::SUCCESS;
    result.state = std::move(candidate);
    return result;
}

std::shared_ptr<const EquilibriumThermoModel>
makeLegacySinglePhaseThermoModel(
    const std::vector<std::string>& speciesNames)
{
    return std::make_shared<const LegacySinglePhaseThermoModel>(speciesNames);
}

} // End namespace thermo
} // End namespace Foam
