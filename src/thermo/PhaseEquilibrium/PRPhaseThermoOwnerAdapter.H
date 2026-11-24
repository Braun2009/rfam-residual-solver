#ifndef PRPhaseThermoOwnerAdapter_H
#define PRPhaseThermoOwnerAdapter_H

#include "CompleteCellState.H"
#include "PRPhaseDisplacementAuthority.H"
#include "PRPhaseRhoERecovery.H"
#include "StateContract.H"

#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace Foam
{
namespace phaseEquilibrium
{

enum class PRPhaseThermoOwnerStatus
{
    SUCCESS,
    INVALID_CONFIGURATION,
    INVALID_PARENT_ADMISSION,
    INVALID_PARENT_STATE_CONTRACT,
    INVALID_PARENT_PHASE_ADMISSION,
    TRANSACTION_ALREADY_ACTIVE,
    NO_ACTIVE_TRANSACTION,
    PARENT_BINDING_DRIFT,
    STALE_TOKEN,
    INVALID_REQUEST,
    DISPLACEMENT_NOT_AUTHORIZED,
    DISPLACEMENT_AUTHORITY_FAILED,
    SEQUENCE_OVERFLOW,
    KERNEL_FAILURE,
    RECOVERY_CONTRACT_FAILED
};

const char* phaseThermoOwnerStatusName(PRPhaseThermoOwnerStatus status);

struct PRPhaseThermoParentAdmission
{
    std::string parentOwnerId;
    std::uint64_t processRuntimeEpoch = 0;
    std::uint64_t parentGeneration = 0;
    PRPhaseIdentity phaseIdentity = PRPhaseIdentity::VAPOR;
    PRPhaseRootIdentity rootIdentity =
        PRPhaseRootIdentity::PR_VAPOR_ROOT;
};

// All three references are parent-owned and immutable for the duration of one
// call. The adapter retains fingerprints and scalar/vector copies only; it
// never stores a pointer to a Domain or publishes a complete state.
struct PRPhaseThermoParentView
{
    const ConservativeStateContract* stateContract = nullptr;
    const CompleteCellState* committed = nullptr;
    const thermo::EquilibriumThermoModel* model = nullptr;
    PRPhaseThermoParentAdmission admission;

    PRPhaseThermoParentView(
        const ConservativeStateContract& parentStateContract,
        const CompleteCellState& parentCommitted,
        const thermo::EquilibriumThermoModel& parentModel,
        const PRPhaseThermoParentAdmission& parentAdmission)
    :
        stateContract(&parentStateContract),
        committed(&parentCommitted),
        model(&parentModel),
        admission(parentAdmission)
    {}
};

struct PRPhaseThermoLiveToken
{
    std::string parentOwnerId;
    std::string transactionId;
    std::uint64_t processRuntimeEpoch = 0;
    std::uint64_t parentGeneration = 0;
    std::uint64_t transactionSequence = 0;
    std::string parentBindingDigest;
    PRPhaseRootIdentity rootIdentity =
        PRPhaseRootIdentity::PR_VAPOR_ROOT;
    std::string tokenDigest;
};

struct PRPhaseThermoOwnerRequest
{
    double phaseDensity = std::numeric_limits<double>::quiet_NaN();
    double phaseSpecificInternalEnergy =
        std::numeric_limits<double>::quiet_NaN();
    std::vector<double> phaseMassFractions;
    std::optional<double> pressureSeed;
    std::optional<double> temperatureSeed;
};

struct PRPhaseThermoOwnerResiduals
{
    double pressureRelative = std::numeric_limits<double>::infinity();
    double temperatureRelative = std::numeric_limits<double>::infinity();
    double compressibilityRelative =
        std::numeric_limits<double>::infinity();
    double densityRelative = std::numeric_limits<double>::infinity();
    double specificInternalEnergyRelative =
        std::numeric_limits<double>::infinity();
    double specificEnthalpyRelative =
        std::numeric_limits<double>::infinity();
    double massFractionL1 = std::numeric_limits<double>::infinity();
    double eosRoot = std::numeric_limits<double>::infinity();
    double departureCaloric = std::numeric_limits<double>::infinity();
    double cvMinusDedTRho = std::numeric_limits<double>::infinity();
    double rootContinuity = std::numeric_limits<double>::infinity();

    bool accepted() const;
};

struct PRPhaseThermoOwnerDiagnostics
{
    bool parentStateContractValidated = false;
    bool parentFingerprintsComputed = false;
    bool inputTokenAdmitted = false;
    bool kernelInvoked = false;
    bool resultTokenIssued = false;
    bool displacementAuthorityInvoked = false;
    bool requestCanonicalizedToAdmittedState = false;
    PRPhaseRhoEStatus kernelStatus =
        PRPhaseRhoEStatus::INVALID_CONFIGURATION;
    PRPhaseDisplacementStatus displacementStatus =
        PRPhaseDisplacementStatus::INVALID_CONFIGURATION;
    std::optional<PRPhaseDisplacementDiagnostics> displacementDiagnostics;
    std::string committedParentQFingerprint;
    std::string parentPhasePayloadFingerprint;
    std::string parentBindingDigest;
};

struct PRPhaseThermoOwnerResult
{
    PRPhaseThermoOwnerStatus status =
        PRPhaseThermoOwnerStatus::INVALID_CONFIGURATION;
    std::optional<PRPhaseCaloricProperties> phaseState;
    std::optional<PRPhaseRhoEDerivatives> derivatives;
    std::optional<PRPhaseThermoLiveToken> acceptedToken;
    PRPhaseThermoOwnerResiduals residuals;
    PRPhaseThermoOwnerDiagnostics diagnostics;

    bool success() const
    {
        return status == PRPhaseThermoOwnerStatus::SUCCESS
            && phaseState.has_value() && derivatives.has_value()
            && acceptedToken.has_value()
            && (diagnostics.displacementAuthorityInvoked
                ? diagnostics.displacementStatus
                    == PRPhaseDisplacementStatus::SUCCESS
                : residuals.accepted());
    }
};

// A numerical solver may evaluate the same trial state repeatedly and in an
// arbitrary order while constructing a Jacobian.  Probe results therefore
// carry no publishable token and never advance or close the live transaction.
struct PRPhaseThermoProbeResult
{
    PRPhaseThermoOwnerStatus status =
        PRPhaseThermoOwnerStatus::INVALID_CONFIGURATION;
    std::optional<PRPhaseCaloricProperties> phaseState;
    std::optional<PRPhaseRhoEDerivatives> derivatives;
    PRPhaseThermoOwnerResiduals residuals;
    PRPhaseThermoOwnerDiagnostics diagnostics;

    bool success() const
    {
        return status == PRPhaseThermoOwnerStatus::SUCCESS
            && phaseState.has_value() && derivatives.has_value()
            && (diagnostics.displacementAuthorityInvoked
                ? diagnostics.displacementStatus
                    == PRPhaseDisplacementStatus::SUCCESS
                : residuals.accepted());
    }
};

// Token-bound adapter between a complete immutable HEM parent and the
// trusted-caller phase rho-e recovery/displacement kernels. Numerical probes
// are pure; only open() and recover() can issue or replace a live token.
class PRPhaseThermoOwnerAdapter
{
private:
    struct ValidatedParent
    {
        const thermo::EquilibriumPhaseProperties* phase = nullptr;
    };

    struct ParentBinding
    {
        std::string committedParentQFingerprint;
        std::string parentPhasePayloadFingerprint;
        std::string parentBindingDigest;
        std::string displacementAuthorityFingerprint;
        std::string parentOwnerId;
        std::string modelId;
        std::uint64_t processRuntimeEpoch = 0;
        std::uint64_t parentGeneration = 0;
        PRPhaseIdentity phaseIdentity = PRPhaseIdentity::VAPOR;
        PRPhaseRootIdentity rootIdentity =
            PRPhaseRootIdentity::PR_VAPOR_ROOT;
    };

    struct ActiveTransaction
    {
        ParentBinding binding;
        std::string transactionId;
        std::uint64_t transactionSequence = 0;
        std::string currentTokenDigest;
        double admittedDensity = std::numeric_limits<double>::quiet_NaN();
        double admittedSpecificInternalEnergy =
            std::numeric_limits<double>::quiet_NaN();
        double admittedTemperature = std::numeric_limits<double>::quiet_NaN();
        std::vector<double> admittedMassFractions;
    };

    const PRPhaseRhoERecovery& recovery_;
    const PRPhaseDisplacementAuthority* displacementAuthority_ = nullptr;
    const thermo::EquilibriumThermoModel* modelAuthority_ = nullptr;
    std::string orderedSpeciesFingerprint_;
    std::string optionsFingerprint_;
    std::string displacementAuthorityFingerprint_;
    bool configurationValid_ = false;
    std::unordered_map<std::string, std::uint64_t> transactionOrdinals_;
    std::optional<ActiveTransaction> active_;

    PRPhaseThermoOwnerStatus validateParent(
        const PRPhaseThermoParentView& parent,
        ValidatedParent& validated,
        PRPhaseThermoOwnerDiagnostics& diagnostics) const;

    ParentBinding fingerprintParent(
        const PRPhaseThermoParentView& parent,
        const ValidatedParent& validated) const;

    PRPhaseThermoOwnerResult recoverBound(
        const PRPhaseThermoParentView& parent,
        const ValidatedParent& validated,
        const ParentBinding& binding,
        const PRPhaseThermoOwnerRequest& request,
        const std::string& transactionId,
        std::uint64_t nextSequence,
        PRPhaseThermoOwnerDiagnostics diagnostics) const;

    PRPhaseThermoOwnerResult recoverDisplacedBound(
        const PRPhaseThermoParentView& parent,
        const ParentBinding& binding,
        const PRPhaseThermoOwnerRequest& request,
        const std::string& transactionId,
        std::uint64_t nextSequence,
        PRPhaseThermoOwnerDiagnostics diagnostics) const;

    PRPhaseThermoLiveToken issueToken(
        const ParentBinding& binding,
        const std::string& transactionId,
        std::uint64_t sequence) const;

    PRPhaseThermoOwnerResult reject(
        PRPhaseThermoOwnerStatus status,
        PRPhaseThermoOwnerDiagnostics diagnostics,
        bool closeTransaction);

    static bool sameBinding(
        const ParentBinding& lhs,
        const ParentBinding& rhs);

    static bool sameToken(
        const PRPhaseThermoLiveToken& token,
        const ActiveTransaction& active);

public:
    explicit PRPhaseThermoOwnerAdapter(
        const PRPhaseRhoERecovery& recovery);
    PRPhaseThermoOwnerAdapter(
        const PRPhaseRhoERecovery& recovery,
        const PRPhaseDisplacementAuthority& displacementAuthority,
        const thermo::EquilibriumThermoModel& modelAuthority);
    PRPhaseThermoOwnerAdapter(PRPhaseRhoERecovery&&) = delete;
    PRPhaseThermoOwnerAdapter(const PRPhaseThermoOwnerAdapter&) = delete;
    PRPhaseThermoOwnerAdapter& operator=(
        const PRPhaseThermoOwnerAdapter&) = delete;
    PRPhaseThermoOwnerAdapter(PRPhaseThermoOwnerAdapter&&) = delete;
    PRPhaseThermoOwnerAdapter& operator=(
        PRPhaseThermoOwnerAdapter&&) = delete;

    PRPhaseThermoOwnerResult open(
        const PRPhaseThermoParentView& parent,
        double temperatureSeed);

    PRPhaseThermoOwnerResult recover(
        const PRPhaseThermoParentView& parent,
        const PRPhaseThermoOwnerRequest& request,
        const PRPhaseThermoLiveToken& token);

    PRPhaseThermoProbeResult probe(
        const PRPhaseThermoParentView& parent,
        const PRPhaseThermoOwnerRequest& request,
        const PRPhaseThermoLiveToken& token) const;

    void closeTransaction() noexcept { active_.reset(); }
    void abort() noexcept { closeTransaction(); }
    bool transactionActive() const noexcept { return active_.has_value(); }
    std::uint64_t currentSequence() const noexcept
    {
        return active_ ? active_->transactionSequence : 0U;
    }
};

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
