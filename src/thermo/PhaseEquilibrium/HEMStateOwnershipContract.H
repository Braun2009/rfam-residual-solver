#ifndef HEMStateOwnershipContract_H
#define HEMStateOwnershipContract_H

#include "BinaryInteractionDB.H"
#include "ConservedVars.H"
#include "EquilibriumThermoModel.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace Foam
{
namespace phaseEquilibrium
{

using thermo::EquilibriumResult;
using thermo::EquilibriumSeed;
using thermo::EquilibriumState;
using thermo::EquilibriumThermoModel;

// This is a narrow transaction contract for future HEM Box B/C harnesses.
// It is deliberately not a transport, phase-change, or production-readiness
// implementation.  The only authoritative state is rho-e-Y; phase fractions
// and phase compositions enter only through the immutable equilibrium cache.
enum class HEMStateOwnershipStatus
{
    SUCCESS,
    INVALID_CONFIGURATION,
    INVALID_COMMITTED_STATE,
    NOT_TWO_PHASE,
    INVALID_VAPOR_SOURCE,
    SOURCE_MASS_NOT_CONSERVATIVE,
    VAPOR_SPECIES_DEPLETED,
    CANDIDATE_SPECIES_INVALID,
    CANDIDATE_ENERGY_INVALID,
    REFLASH_FAILED,
    REFLASH_RESIDUAL_MISMATCH
};

inline const char* hemStateOwnershipStatusName(
    HEMStateOwnershipStatus status)
{
    switch (status)
    {
        case HEMStateOwnershipStatus::SUCCESS: return "SUCCESS";
        case HEMStateOwnershipStatus::INVALID_CONFIGURATION:
            return "INVALID_CONFIGURATION";
        case HEMStateOwnershipStatus::INVALID_COMMITTED_STATE:
            return "INVALID_COMMITTED_STATE";
        case HEMStateOwnershipStatus::NOT_TWO_PHASE:
            return "NOT_TWO_PHASE";
        case HEMStateOwnershipStatus::INVALID_VAPOR_SOURCE:
            return "INVALID_VAPOR_SOURCE";
        case HEMStateOwnershipStatus::SOURCE_MASS_NOT_CONSERVATIVE:
            return "SOURCE_MASS_NOT_CONSERVATIVE";
        case HEMStateOwnershipStatus::VAPOR_SPECIES_DEPLETED:
            return "VAPOR_SPECIES_DEPLETED";
        case HEMStateOwnershipStatus::CANDIDATE_SPECIES_INVALID:
            return "CANDIDATE_SPECIES_INVALID";
        case HEMStateOwnershipStatus::CANDIDATE_ENERGY_INVALID:
            return "CANDIDATE_ENERGY_INVALID";
        case HEMStateOwnershipStatus::REFLASH_FAILED:
            return "REFLASH_FAILED";
        case HEMStateOwnershipStatus::REFLASH_RESIDUAL_MISMATCH:
            return "REFLASH_RESIDUAL_MISMATCH";
    }
    return "INVALID_CONFIGURATION";
}

struct HEMStateOwnershipTolerances
{
    static constexpr double maximumConcentrationAbsolute = 1.0e-6;
    static constexpr double maximumEnergyAbsolute = 1.0;
    static constexpr double maximumSourceMassAbsolute = 1.0e-9;
    static constexpr double maximumRelative = 1.0e-6;

    double concentrationAbsolute = 1.0e-10;
    double concentrationRelative = 1.0e-8;
    double energyAbsolute = 1.0e-5;
    double energyRelative = 1.0e-8;
    double sourceMassAbsolute = 1.0e-12;
    double sourceMassRelative = 1.0e-10;

    bool valid() const
    {
        return std::isfinite(concentrationAbsolute)
            && std::isfinite(concentrationRelative)
            && std::isfinite(energyAbsolute)
            && std::isfinite(energyRelative)
            && std::isfinite(sourceMassAbsolute)
            && std::isfinite(sourceMassRelative)
            && concentrationAbsolute >= 0.0
            && concentrationAbsolute <= maximumConcentrationAbsolute
            && concentrationRelative >= 0.0
            && concentrationRelative <= maximumRelative
            && energyAbsolute >= 0.0
            && energyAbsolute <= maximumEnergyAbsolute
            && energyRelative >= 0.0
            && energyRelative <= maximumRelative
            && sourceMassAbsolute >= 0.0
            && sourceMassAbsolute <= maximumSourceMassAbsolute
            && sourceMassRelative >= 0.0
            && sourceMassRelative <= maximumRelative;
    }
};

// Each entry is an integrated source increment per unit vapor volume.  It is
// not a transported phase state: alpha_v is read from the committed flash and
// converts the increment to the mixture-volume conserved update.
struct HEMVaporConservativeIncrement
{
    std::vector<double> deltaRhoYPerVaporVolume;
    double deltaRhoEPerVaporVolume =
        std::numeric_limits<double>::quiet_NaN();
};

struct HEMStateOwnershipResult
{
    HEMStateOwnershipStatus status =
        HEMStateOwnershipStatus::INVALID_CONFIGURATION;
    std::optional<ConservedVars> candidate;
    std::optional<EquilibriumState> reflashed;

    bool success() const
    {
        return status == HEMStateOwnershipStatus::SUCCESS
            && candidate.has_value() && reflashed.has_value()
            && reflashed->valid();
    }
};

class HEMRhoEStateOwnershipTransaction
{
    const EquilibriumThermoModel& model_;
    HEMStateOwnershipTolerances tolerances_;

    static bool close(
        double lhs, double rhs, double absolute, double relative)
    {
        return std::isfinite(lhs) && std::isfinite(rhs)
            && std::abs(lhs - rhs) <= absolute + relative
                * std::max(std::abs(lhs), std::abs(rhs));
    }

    static double inventoryRoundoff(
        double available, double increment)
    {
        return 64.0*std::numeric_limits<double>::epsilon()
            *std::max({std::abs(available), std::abs(increment), 1.0});
    }

    static bool normalizedCandidate(
        const ConservedVars& candidate, std::vector<double>& massFractions,
        const HEMStateOwnershipTolerances& tolerances)
    {
        if (!std::isfinite(candidate.rho) || candidate.rho <= 0.0
            || candidate.rhoY.empty())
            return false;

        massFractions.assign(candidate.rhoY.size(), 0.0);
        double sum = 0.0;
        for (std::size_t species = 0;
             species < candidate.rhoY.size(); ++species)
        {
            const double rhoY = candidate.rhoY[species];
            if (!std::isfinite(rhoY) || rhoY < 0.0) return false;
            massFractions[species] = rhoY/candidate.rho;
            if (!std::isfinite(massFractions[species])) return false;
            sum += massFractions[species];
        }
        return close(
            sum, 1.0, tolerances.concentrationAbsolute,
            tolerances.concentrationRelative);
    }

    static bool internalEnergy(
        const ConservedVars& state, double& result)
    {
        if (!std::isfinite(state.rho) || state.rho <= 0.0
            || !std::isfinite(state.rhoE) || state.rhoU.size() != 3)
            return false;
        double kinetic = 0.0;
        for (double momentum : state.rhoU)
        {
            if (!std::isfinite(momentum)) return false;
            const double velocity = momentum/state.rho;
            kinetic += 0.5*velocity*velocity;
        }
        result = state.rhoE/state.rho - kinetic;
        return std::isfinite(result);
    }

    bool committedStateValid(
        const ConservedVars& committed,
        const EquilibriumState& equilibrium) const
    {
        if (!equilibrium.valid() || !model_.stateValidForModel(equilibrium)
            || committed.rhoY.size()
                != equilibrium.overallMassFractions.size())
            return false;

        if (!close(
                committed.rho, equilibrium.rho,
                tolerances_.concentrationAbsolute,
                tolerances_.concentrationRelative))
            return false;

        double committedEnergy = 0.0;
        if (!internalEnergy(committed, committedEnergy)
            || !close(
                committedEnergy, equilibrium.e,
                tolerances_.energyAbsolute, tolerances_.energyRelative))
            return false;

        for (std::size_t species = 0;
             species < committed.rhoY.size(); ++species)
        {
            if (!close(
                    committed.rhoY[species],
                    committed.rho*equilibrium.overallMassFractions[species],
                    tolerances_.concentrationAbsolute,
                    tolerances_.concentrationRelative))
                return false;
        }
        return true;
    }

public:
    explicit HEMRhoEStateOwnershipTransaction(
        const EquilibriumThermoModel& model,
        const HEMStateOwnershipTolerances& tolerances = {})
        : model_(model), tolerances_(tolerances)
    {
    }

    // Does not mutate committed.  The caller may publish candidate/reflashed
    // only after success, preserving full-step rollback semantics.
    HEMStateOwnershipResult stageVaporConservativeUpdate(
        const ConservedVars& committed,
        const EquilibriumState& equilibrium,
        const HEMVaporConservativeIncrement& vaporIncrement) const
    {
        HEMStateOwnershipResult result;
        if (!tolerances_.valid())
        {
            result.status = HEMStateOwnershipStatus::INVALID_CONFIGURATION;
            return result;
        }
        if (equilibrium.region != thermo::EquilibriumRegion::TWO_PHASE)
        {
            result.status = HEMStateOwnershipStatus::NOT_TWO_PHASE;
            return result;
        }
        if (!equilibrium.vapor.has_value()
            || !std::isfinite(equilibrium.vaporVolumeFraction)
            || equilibrium.vaporVolumeFraction <= 0.0
            || equilibrium.vaporVolumeFraction > 1.0
            || !committedStateValid(committed, equilibrium))
        {
            result.status = HEMStateOwnershipStatus::INVALID_COMMITTED_STATE;
            return result;
        }
        if (vaporIncrement.deltaRhoYPerVaporVolume.size()
                != committed.rhoY.size()
            || !std::isfinite(vaporIncrement.deltaRhoEPerVaporVolume))
        {
            result.status = HEMStateOwnershipStatus::INVALID_VAPOR_SOURCE;
            return result;
        }

        double sourceMass = 0.0;
        for (double increment : vaporIncrement.deltaRhoYPerVaporVolume)
        {
            if (!std::isfinite(increment))
            {
                result.status = HEMStateOwnershipStatus::INVALID_VAPOR_SOURCE;
                return result;
            }
            sourceMass += increment;
        }
        if (!close(
                sourceMass, 0.0, tolerances_.sourceMassAbsolute,
                tolerances_.sourceMassRelative))
        {
            result.status = HEMStateOwnershipStatus::SOURCE_MASS_NOT_CONSERVATIVE;
            return result;
        }

        const thermo::EquilibriumPhaseProperties& vapor = *equilibrium.vapor;
        if (vapor.massFractions.size() != committed.rhoY.size()
            || !std::isfinite(vapor.rho) || vapor.rho <= 0.0)
        {
            result.status = HEMStateOwnershipStatus::INVALID_COMMITTED_STATE;
            return result;
        }
        for (std::size_t species = 0;
             species < vaporIncrement.deltaRhoYPerVaporVolume.size();
             ++species)
        {
            const double available =
                vapor.rho*vapor.massFractions[species];
            const double increment =
                vaporIncrement.deltaRhoYPerVaporVolume[species];
            if (!std::isfinite(available) || available < 0.0
                || available + increment
                    < -inventoryRoundoff(available, increment))
            {
                result.status = HEMStateOwnershipStatus::VAPOR_SPECIES_DEPLETED;
                return result;
            }
        }

        ConservedVars candidate = committed;
        const double vaporVolumeFraction = equilibrium.vaporVolumeFraction;
        for (std::size_t species = 0;
             species < candidate.rhoY.size(); ++species)
        {
            candidate.rhoY[species] += vaporVolumeFraction
                * vaporIncrement.deltaRhoYPerVaporVolume[species];
        }
        candidate.rhoE += vaporVolumeFraction
            * vaporIncrement.deltaRhoEPerVaporVolume;

        std::vector<double> candidateY;
        if (!normalizedCandidate(candidate, candidateY, tolerances_))
        {
            result.status = HEMStateOwnershipStatus::CANDIDATE_SPECIES_INVALID;
            return result;
        }

        double candidateEnergy = 0.0;
        if (!internalEnergy(candidate, candidateEnergy))
        {
            result.status = HEMStateOwnershipStatus::CANDIDATE_ENERGY_INVALID;
            return result;
        }

        EquilibriumSeed seed;
        seed.pressure = equilibrium.p;
        seed.temperature = equilibrium.T;
        const EquilibriumResult recovered = model_.recoverRhoE(
            candidate.rho, candidateEnergy, candidateY, &seed);
        if (!recovered.success() || !model_.stateValidForModel(*recovered.state))
        {
            result.status = HEMStateOwnershipStatus::REFLASH_FAILED;
            return result;
        }

        const EquilibriumState& reflashed = *recovered.state;
        if (!close(
                reflashed.rho, candidate.rho,
                tolerances_.concentrationAbsolute,
                tolerances_.concentrationRelative)
            || !close(
                reflashed.e, candidateEnergy, tolerances_.energyAbsolute,
                tolerances_.energyRelative)
            || reflashed.overallMassFractions.size() != candidateY.size())
        {
            result.status = HEMStateOwnershipStatus::REFLASH_RESIDUAL_MISMATCH;
            return result;
        }
        for (std::size_t species = 0; species < candidateY.size(); ++species)
        {
            if (!close(
                    reflashed.overallMassFractions[species], candidateY[species],
                    tolerances_.concentrationAbsolute,
                    tolerances_.concentrationRelative))
            {
                result.status = HEMStateOwnershipStatus::REFLASH_RESIDUAL_MISMATCH;
                return result;
            }
        }

        result.status = HEMStateOwnershipStatus::SUCCESS;
        result.candidate = std::move(candidate);
        result.reflashed = reflashed;
        return result;
    }
};

enum class HEMBoxPreconditionStatus
{
    ELIGIBLE_FOR_NONPRODUCTION_BOX,
    OWNERSHIP_TRANSACTION_UNVERIFIED,
    O2_C3H8_KIJ_UNVALIDATED,
    PHASE_TOPOLOGY_UNVALIDATED,
    BOX_B_UNVALIDATED,
    R2_PHYSICS_UNVERIFIED
};

inline const char* hemBoxPreconditionStatusName(
    HEMBoxPreconditionStatus status)
{
    switch (status)
    {
        case HEMBoxPreconditionStatus::ELIGIBLE_FOR_NONPRODUCTION_BOX:
            return "ELIGIBLE_FOR_NONPRODUCTION_BOX";
        case HEMBoxPreconditionStatus::OWNERSHIP_TRANSACTION_UNVERIFIED:
            return "OWNERSHIP_TRANSACTION_UNVERIFIED";
        case HEMBoxPreconditionStatus::O2_C3H8_KIJ_UNVALIDATED:
            return "O2_C3H8_KIJ_UNVALIDATED";
        case HEMBoxPreconditionStatus::PHASE_TOPOLOGY_UNVALIDATED:
            return "PHASE_TOPOLOGY_UNVALIDATED";
        case HEMBoxPreconditionStatus::BOX_B_UNVALIDATED:
            return "BOX_B_UNVALIDATED";
        case HEMBoxPreconditionStatus::R2_PHYSICS_UNVERIFIED:
            return "R2_PHYSICS_UNVERIFIED";
    }
    return "OWNERSHIP_TRANSACTION_UNVERIFIED";
}

enum class HEMPhaseTopologyEvidence
{
    UNVALIDATED,
    TEST_ONLY,
    VALIDATED_MAXIMUM_TWO_PHASES
};

// A converged two-phase PR flash does not prove that a third equilibrium
// phase is absent. Production eligibility therefore requires independent
// phase-topology evidence covering both the thermodynamic and composition
// envelope. This record is only an input gate; it does not implement VLLE.
struct HEMPhaseTopologyReadiness
{
    HEMPhaseTopologyEvidence evidence =
        HEMPhaseTopologyEvidence::UNVALIDATED;
    std::string firstSpecies;
    std::string secondSpecies;
    BinaryInteractionApplicability applicability;
    double minimumO2MoleFraction = std::numeric_limits<double>::quiet_NaN();
    double maximumO2MoleFraction = std::numeric_limits<double>::quiet_NaN();
    std::size_t maximumEquilibriumPhaseCount = 0;
    bool vlleObservedOrPredicted = false;
    std::string sourceId;
    std::string datasetId;

    bool structurallyValid() const
    {
        return !firstSpecies.empty() && !secondSpecies.empty()
            && firstSpecies < secondSpecies
            && applicability.valid()
            && std::isfinite(minimumO2MoleFraction)
            && std::isfinite(maximumO2MoleFraction)
            && minimumO2MoleFraction >= 0.0
            && maximumO2MoleFraction <= 1.0
            && minimumO2MoleFraction <= maximumO2MoleFraction
            && maximumEquilibriumPhaseCount > 0
            && !sourceId.empty() && !datasetId.empty();
    }

    bool productionReady(
        const std::string& requiredFirstSpecies,
        const std::string& requiredSecondSpecies,
        const BinaryInteractionApplicability& requiredApplicability,
        double requiredMinimumO2MoleFraction,
        double requiredMaximumO2MoleFraction) const
    {
        return structurallyValid()
            && firstSpecies == requiredFirstSpecies
            && secondSpecies == requiredSecondSpecies
            && requiredApplicability.valid()
            && std::isfinite(requiredMinimumO2MoleFraction)
            && std::isfinite(requiredMaximumO2MoleFraction)
            && requiredMinimumO2MoleFraction >= 0.0
            && requiredMaximumO2MoleFraction <= 1.0
            && requiredMinimumO2MoleFraction
                <= requiredMaximumO2MoleFraction
            && evidence
                == HEMPhaseTopologyEvidence::VALIDATED_MAXIMUM_TWO_PHASES
            && maximumEquilibriumPhaseCount <= 2
            && !vlleObservedOrPredicted
            && applicability.covers(requiredApplicability)
            && minimumO2MoleFraction <= requiredMinimumO2MoleFraction
            && maximumO2MoleFraction >= requiredMaximumO2MoleFraction;
    }
};

struct HEMBoxBPrerequisites
{
    bool rhoEYOwnershipTransactionValidated = false;
    BinaryInteractionReadinessAudit o2C3H8KijReadiness;
    BinaryInteractionApplicability requiredApplicability;
    double requiredMinimumO2MoleFraction =
        std::numeric_limits<double>::quiet_NaN();
    double requiredMaximumO2MoleFraction =
        std::numeric_limits<double>::quiet_NaN();
    HEMPhaseTopologyReadiness phaseTopology;
};

inline bool exactBinaryInteractionAuditReady(
    const BinaryInteractionReadinessAudit& audit,
    const std::string& requiredFirstSpecies,
    const std::string& requiredSecondSpecies,
    const BinaryInteractionApplicability& requiredApplicability)
{
    return audit.productionReady()
        && audit.pairs.size() == 1
        && audit.pairs[0].first == requiredFirstSpecies
        && audit.pairs[0].second == requiredSecondSpecies
        && audit.requiredEnvelope.covers(requiredApplicability);
}

inline HEMBoxPreconditionStatus evaluateHEMBoxBPrerequisites(
    const HEMBoxBPrerequisites& prerequisites)
{
    const std::string firstSpecies = "C3H8";
    const std::string secondSpecies = "O2";
    if (!prerequisites.rhoEYOwnershipTransactionValidated)
        return HEMBoxPreconditionStatus::OWNERSHIP_TRANSACTION_UNVERIFIED;
    if (!exactBinaryInteractionAuditReady(
            prerequisites.o2C3H8KijReadiness,
            firstSpecies,
            secondSpecies,
            prerequisites.requiredApplicability))
        return HEMBoxPreconditionStatus::O2_C3H8_KIJ_UNVALIDATED;
    if (!prerequisites.phaseTopology.productionReady(
            firstSpecies,
            secondSpecies,
            prerequisites.requiredApplicability,
            prerequisites.requiredMinimumO2MoleFraction,
            prerequisites.requiredMaximumO2MoleFraction))
        return HEMBoxPreconditionStatus::PHASE_TOPOLOGY_UNVALIDATED;
    return HEMBoxPreconditionStatus::ELIGIBLE_FOR_NONPRODUCTION_BOX;
}

struct HEMBoxCPrerequisites
{
    HEMBoxPreconditionStatus boxBPrerequisites =
        HEMBoxPreconditionStatus::BOX_B_UNVALIDATED;
    bool boxBConservativeMixingValidated = false;
    bool r2PhysicsValidated = false;
};

inline HEMBoxPreconditionStatus evaluateHEMBoxCPrerequisites(
    const HEMBoxCPrerequisites& prerequisites)
{
    if (prerequisites.boxBPrerequisites
            != HEMBoxPreconditionStatus::ELIGIBLE_FOR_NONPRODUCTION_BOX
        || !prerequisites.boxBConservativeMixingValidated)
        return HEMBoxPreconditionStatus::BOX_B_UNVALIDATED;
    if (!prerequisites.r2PhysicsValidated)
        return HEMBoxPreconditionStatus::R2_PHYSICS_UNVERIFIED;
    return HEMBoxPreconditionStatus::ELIGIBLE_FOR_NONPRODUCTION_BOX;
}

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
