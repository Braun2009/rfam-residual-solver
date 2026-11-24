#ifndef SectorRestartPublicationStagers_H
#define SectorRestartPublicationStagers_H

#include "RestartPublication.H"
#include "SectorRestartOwnerCodec.H"

#include <filesystem>
#include <memory>

namespace Foam
{

class SectorTurbulenceOperator;

class SectorRestartDetachedFluidState final
    : public RestartDetachedOwnerState
{
    std::shared_ptr<const Domain> domain_;
    RestartOwnerRecord record_;

public:
    SectorRestartDetachedFluidState(
        std::unique_ptr<Domain> domain,
        RestartOwnerRecord record);

    const Domain& domain() const noexcept
    {
        return *domain_;
    }

    const RestartOwnerRecord& record() const noexcept
    {
        return record_;
    }
};

class SectorRestartDetachedTurbulenceState final
    : public RestartDetachedOwnerState
{
    std::shared_ptr<const SectorRestartTurbulenceState> fields_;
    RestartOwnerRecord record_;

public:
    SectorRestartDetachedTurbulenceState(
        std::unique_ptr<SectorRestartTurbulenceState> fields,
        RestartOwnerRecord record);

    const SectorRestartTurbulenceState& fields() const noexcept
    {
        return *fields_;
    }

    const RestartOwnerRecord& record() const noexcept
    {
        return record_;
    }
};

struct SectorRestartFluidStageContext
{
    const Domain* targetTemplate = nullptr;
};

struct SectorRestartTurbulenceStageContext
{
    const Domain* targetDomain = nullptr;
    const SectorTurbulenceOperator* targetTurbulence = nullptr;
};

struct SectorRestartPublicationCrossValidationContext
{
    const std::string* explicitGenerationId = nullptr;
    const RestartManifestExpectation* expected = nullptr;
    const Domain* targetTemplate = nullptr;
    const SectorTurbulenceOperator* targetTurbulence = nullptr;
};

// The caller supplies an exact generation ID and its already accepted ledger
// entry. This request deliberately does not discover or select "latest".
struct SectorRestartExplicitGenerationStageRequest
{
    std::filesystem::path storeRoot;
    std::string explicitGenerationId;
    RestartGenerationLedgerEntry selectedLedgerEntry;
    RestartManifestExpectation expected;
    RestartCheckpointClockObservation checkpoint;
    const Domain* targetTemplate = nullptr;
    const SectorTurbulenceOperator* targetTurbulence = nullptr;
};

struct SectorRestartExplicitGenerationStageResult
{
    RestartGenerationStoreStatus loadStatus =
        RestartGenerationStoreStatus::INVALID_INPUT;
    bool loadAttempted = false;
    std::string loadDetail;
    RestartPublicationStageResult staging;

    bool success() const noexcept
    {
        return loadAttempted
            && loadStatus == RestartGenerationStoreStatus::SUCCESS
            && staging.success();
    }
};

RestartOwnerStageResult stageSectorRestartFluidOwner(
    const RestartOwnerStageInput& input,
    void* context);

RestartOwnerStageResult stageSectorRestartTurbulenceOwner(
    const RestartOwnerStageInput& input,
    void* context);

RestartCrossValidationResult validateSectorRestartPublicationRoot(
    const RestartPublicationRoot& root,
    void* context);

// Candidate-only: performs an actual stable generation-store load followed by
// detached fluid/turbulence staging. It never publishes a live root.
SectorRestartExplicitGenerationStageResult
stageSectorRestartExplicitGeneration(
    const std::shared_ptr<const RestartPublicationRoot>& currentRoot,
    const SectorRestartExplicitGenerationStageRequest& request) noexcept;

} // End namespace Foam

#endif
