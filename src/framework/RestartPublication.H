#ifndef RestartPublication_H
#define RestartPublication_H

#include "RestartClockAdmission.H"

#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace Foam
{

enum class RestartPublicationStatus : std::uint8_t
{
    SUCCESS = 0,
    INVALID_REQUEST,
    INVALID_STORE_RESULT,
    INVALID_BASE_ROOT,
    LEDGER_REJECTED,
    LEDGER_CONFLICT,
    MANIFEST_REJECTED,
    CLOCK_REJECTED,
    OWNER_STAGER_MISSING,
    OWNER_STAGER_UNEXPECTED,
    OWNER_STAGE_REJECTED,
    OWNER_STAGE_EXCEPTION,
    OWNER_BINDING_MISMATCH,
    OWNER_LEDGER_MISMATCH,
    CROSS_VALIDATOR_MISSING,
    CROSS_VALIDATION_REJECTED,
    CROSS_VALIDATION_EXCEPTION,
    STALE_PUBLICATION_BASE,
    PUBLICATION_ALREADY_CONSUMED,
    INTERNAL_ERROR
};

const char* restartPublicationStatusName(
    RestartPublicationStatus status) noexcept;

enum class RestartGenerationLedgerDisposition : std::uint8_t
{
    ACCEPTED = 0,
    REJECTED
};

struct RestartGenerationLedgerEntry
{
    RestartGenerationLedgerDisposition disposition =
        RestartGenerationLedgerDisposition::REJECTED;
    std::string generationId;
    std::uint64_t ordinal = 0;
    std::string generationChecksumSha256;
    RestartPublicationStatus status = RestartPublicationStatus::INVALID_REQUEST;
    std::optional<RestartOwner> owner;
    std::string detail;

    bool accepted() const noexcept
    {
        return disposition == RestartGenerationLedgerDisposition::ACCEPTED
            && status == RestartPublicationStatus::SUCCESS;
    }

    bool rejected() const noexcept
    {
        return disposition == RestartGenerationLedgerDisposition::REJECTED
            && status != RestartPublicationStatus::SUCCESS;
    }
};

bool validRestartGenerationLedgerEntry(
    const RestartGenerationLedgerEntry& entry) noexcept;

bool validRestartAcceptedGenerationLedger(
    const std::vector<RestartGenerationLedgerEntry>& entries) noexcept;

// Diagnostic-only. Entries in this ledger are never generation-selection
// authority and are never copied into the published live-state root.
class RestartRejectedGenerationLedger
{
    std::vector<RestartGenerationLedgerEntry> entries_;

public:
    bool record(const RestartGenerationLedgerEntry& entry) noexcept;

    const std::vector<RestartGenerationLedgerEntry>& entries() const noexcept
    {
        return entries_;
    }
};

struct RestartDomainStepLedger
{
    std::uint64_t acceptedSteps = 0;
    std::uint64_t rejectedSteps = 0;

    bool consistentWith(
        const RestartGenerationState& generation) const noexcept
    {
        return acceptedSteps == generation.acceptedIteration;
    }
};

// Owner implementations decode into immutable detached objects. No callback
// receives a live Domain, SolverBase, or owner object. The fluid owner must
// supply the Domain step ledger through this base; no other owner may supply
// one. The payload checksum remains the authority for the decoded values.
class RestartDetachedOwnerState
{
    std::optional<RestartDomainStepLedger> domainStepLedger_;

protected:
    explicit RestartDetachedOwnerState(
        std::optional<RestartDomainStepLedger> domainStepLedger = std::nullopt)
        : domainStepLedger_(std::move(domainStepLedger))
    {}

public:
    virtual ~RestartDetachedOwnerState() = default;

    const std::optional<RestartDomainStepLedger>& domainStepLedger()
        const noexcept
    {
        return domainStepLedger_;
    }
};

struct RestartOwnerStageInput
{
    const RestartManifest& manifest;
    const RestartOwnerRecord& owner;
    const RestartPayloadCandidate& payload;
};

struct RestartOwnerStageResult
{
    RestartOwner owner = RestartOwner::FLUID;
    std::string generationId;
    std::string identitySha256;
    std::string payloadSha256;
    std::shared_ptr<const RestartDetachedOwnerState> state;
    std::string detail;

    bool success() const noexcept
    {
        return state != nullptr;
    }
};

using RestartOwnerStageFunction = RestartOwnerStageResult (*)(
    const RestartOwnerStageInput& input,
    void* context);

struct RestartOwnerStager
{
    RestartOwner owner = RestartOwner::FLUID;
    RestartOwnerStageFunction stage = nullptr;
    void* context = nullptr;
};

struct RestartPublishedOwner
{
    RestartOwnerRecord record;
    std::string stagedGenerationId;
    std::string stagedIdentitySha256;
    std::string stagedPayloadSha256;
    std::shared_ptr<const RestartDetachedOwnerState> state;

    bool valid() const noexcept;
};

class RestartPublicationRoot
{
    RestartManifest manifest_;
    RestartClockPublication clock_;
    std::array<RestartPublishedOwner,restartOwnerCount> owners_{};
    std::vector<RestartGenerationLedgerEntry> acceptedLedger_;

    RestartPublicationRoot() = default;

    friend class RestartPublicationStager;

public:
    const RestartManifest& manifest() const noexcept
    {
        return manifest_;
    }

    const RestartClockPublication& clock() const noexcept
    {
        return clock_;
    }

    const std::array<RestartPublishedOwner,restartOwnerCount>& owners()
        const noexcept
    {
        return owners_;
    }

    const std::vector<RestartGenerationLedgerEntry>& acceptedLedger()
        const noexcept
    {
        return acceptedLedger_;
    }

    bool valid() const noexcept;
};

struct RestartCrossValidationResult
{
    bool accepted = false;
    std::string detail;
};

using RestartCrossValidationFunction = RestartCrossValidationResult (*)(
    const RestartPublicationRoot& root,
    void* context);

struct RestartPublicationRequest
{
    RestartGenerationLedgerEntry ledgerEntry;
    std::array<RestartOwnerStager,restartOwnerCount> ownerStagers{};
    RestartCrossValidationFunction crossValidate = nullptr;
    void* crossValidationContext = nullptr;
};

struct RestartPublicationCommitResult
{
    RestartPublicationStatus status = RestartPublicationStatus::INTERNAL_ERROR;
    std::string detail;
    std::optional<RestartGenerationLedgerEntry> rejectedLedgerEntry;

    bool success() const noexcept
    {
        return status == RestartPublicationStatus::SUCCESS;
    }
};

class RestartStagedPublication
{
    std::shared_ptr<const RestartPublicationRoot> expectedBase_;
    std::shared_ptr<const RestartPublicationRoot> candidate_;
    RestartGenerationLedgerEntry ledgerEntry_;
    bool consumed_ = false;

    RestartStagedPublication(
        std::shared_ptr<const RestartPublicationRoot> expectedBase,
        std::shared_ptr<const RestartPublicationRoot> candidate,
        RestartGenerationLedgerEntry ledgerEntry);

    friend class RestartPublicationStager;

public:
    RestartStagedPublication(const RestartStagedPublication&) = delete;
    RestartStagedPublication& operator=(const RestartStagedPublication&) = delete;

    const RestartPublicationRoot& candidate() const noexcept
    {
        return *candidate_;
    }

    RestartPublicationCommitResult publish(
        std::shared_ptr<const RestartPublicationRoot>& liveRoot) noexcept;
};

struct RestartPublicationStageResult
{
    RestartPublicationStatus status = RestartPublicationStatus::INTERNAL_ERROR;
    std::string detail;
    std::optional<RestartGenerationLedgerEntry> rejectedLedgerEntry;
    std::unique_ptr<RestartStagedPublication> publication;

    bool success() const noexcept
    {
        return status == RestartPublicationStatus::SUCCESS
            && publication != nullptr && !rejectedLedgerEntry.has_value();
    }
};

// Candidate-only coordinator for an offline immutable-root contract. It is not
// wired to SectorCheckpoint, SolverBase, the live Domain, generation selection,
// or a quiescence mechanism, and it grants no restart or production capability.
class RestartPublicationStager
{
public:
    static RestartPublicationStageResult stage(
        const std::shared_ptr<const RestartPublicationRoot>& currentRoot,
        const RestartGenerationStoreLoadResult& load,
        const RestartManifestExpectation& expected,
        const RestartCheckpointClockObservation& checkpoint,
        const RestartPublicationRequest& request) noexcept;
};

} // End namespace Foam

#endif
