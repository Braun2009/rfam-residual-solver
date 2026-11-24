#include "RestartPublication.H"

#include <exception>
#include <utility>

namespace Foam
{
namespace
{

constexpr std::size_t maximumLedgerDetailBytes = 4096U;

bool validOwner(RestartOwner owner) noexcept
{
    return static_cast<std::size_t>(owner) < restartOwnerCount;
}

bool validPublicationStatus(RestartPublicationStatus status) noexcept
{
    switch(status)
    {
        case RestartPublicationStatus::SUCCESS:
        case RestartPublicationStatus::INVALID_REQUEST:
        case RestartPublicationStatus::INVALID_STORE_RESULT:
        case RestartPublicationStatus::INVALID_BASE_ROOT:
        case RestartPublicationStatus::LEDGER_REJECTED:
        case RestartPublicationStatus::LEDGER_CONFLICT:
        case RestartPublicationStatus::MANIFEST_REJECTED:
        case RestartPublicationStatus::CLOCK_REJECTED:
        case RestartPublicationStatus::OWNER_STAGER_MISSING:
        case RestartPublicationStatus::OWNER_STAGER_UNEXPECTED:
        case RestartPublicationStatus::OWNER_STAGE_REJECTED:
        case RestartPublicationStatus::OWNER_STAGE_EXCEPTION:
        case RestartPublicationStatus::OWNER_BINDING_MISMATCH:
        case RestartPublicationStatus::OWNER_LEDGER_MISMATCH:
        case RestartPublicationStatus::CROSS_VALIDATOR_MISSING:
        case RestartPublicationStatus::CROSS_VALIDATION_REJECTED:
        case RestartPublicationStatus::CROSS_VALIDATION_EXCEPTION:
        case RestartPublicationStatus::STALE_PUBLICATION_BASE:
        case RestartPublicationStatus::PUBLICATION_ALREADY_CONSUMED:
        case RestartPublicationStatus::INTERNAL_ERROR:
            return true;
    }
    return false;
}

std::string boundedDetail(
    const std::string& detail,
    const char* fallback)
{
    const std::string& source = detail.empty() ? std::string(fallback) : detail;
    if(source.size() <= maximumLedgerDetailBytes) return source;
    return source.substr(0U,maximumLedgerDetailBytes);
}

std::optional<RestartGenerationLedgerEntry> rejectedEntry(
    const RestartGenerationLedgerEntry& source,
    RestartPublicationStatus status,
    const std::optional<RestartOwner>& owner,
    const std::string& detail)
{
    if(!RestartGenerationStore::validGenerationId(source.generationId)
       ||!validRestartSha256(source.generationChecksumSha256)
       ||status==RestartPublicationStatus::SUCCESS)
        return std::nullopt;

    RestartGenerationLedgerEntry rejected;
    rejected.disposition=RestartGenerationLedgerDisposition::REJECTED;
    rejected.generationId=source.generationId;
    rejected.ordinal=source.ordinal;
    rejected.generationChecksumSha256=source.generationChecksumSha256;
    rejected.status=status;
    rejected.owner=owner;
    rejected.detail=boundedDetail(detail,"restart publication rejected");
    if(!validRestartGenerationLedgerEntry(rejected)) return std::nullopt;
    return rejected;
}

RestartPublicationStageResult stageFailure(
    RestartPublicationStatus status,
    const std::string& detail,
    const RestartGenerationLedgerEntry* ledgerEntry = nullptr,
    const std::optional<RestartOwner>& owner = std::nullopt)
{
    RestartPublicationStageResult result;
    result.status=status;
    result.detail=detail;
    if(ledgerEntry)
        result.rejectedLedgerEntry=
            rejectedEntry(*ledgerEntry,status,owner,detail);
    return result;
}

RestartPublicationCommitResult commitFailure(
    RestartPublicationStatus status,
    const std::string& detail,
    const RestartGenerationLedgerEntry* ledgerEntry = nullptr)
{
    RestartPublicationCommitResult result;
    result.status=status;
    result.detail=detail;
    if(ledgerEntry)
        result.rejectedLedgerEntry=
            rejectedEntry(*ledgerEntry,status,std::nullopt,detail);
    return result;
}

const RestartPayloadCandidate* findPayload(
    const RestartCandidatePackage& candidate,
    RestartOwner owner) noexcept
{
    for(const RestartPayloadCandidate& payload:candidate.payloads)
        if(payload.owner==owner) return &payload;
    return nullptr;
}

bool sameOwnerRecord(
    const RestartOwnerRecord& left,
    const RestartOwnerRecord& right) noexcept
{
    return left.owner==right.owner
        &&left.disposition==right.disposition
        &&left.publicationOrder==right.publicationOrder
        &&left.generationId==right.generationId
        &&left.identitySha256==right.identitySha256
        &&left.payloadBytes==right.payloadBytes
        &&left.payloadSha256==right.payloadSha256;
}

bool sameClock(
    const RestartClockPublication& clock,
    const RestartGenerationState& generation) noexcept
{
    if(clock.physicalTime!=generation.physicalTime
       ||clock.acceptedIteration
            !=static_cast<int>(generation.acceptedIteration)
       ||clock.previousAcceptedDt.has_value()
            !=generation.previousAcceptedDtPresent)
        return false;
    return !clock.previousAcceptedDt.has_value()
        ||*clock.previousAcceptedDt==generation.previousAcceptedDt;
}

bool ledgerAllows(
    const std::vector<RestartGenerationLedgerEntry>& accepted,
    const RestartGenerationLedgerEntry& proposed) noexcept
{
    if(!validRestartAcceptedGenerationLedger(accepted)||!proposed.accepted())
        return false;
    if(!accepted.empty()&&proposed.ordinal<=accepted.back().ordinal)
        return false;
    for(const RestartGenerationLedgerEntry& entry:accepted)
        if(entry.generationId==proposed.generationId
           ||entry.ordinal==proposed.ordinal)
            return false;
    return true;
}

} // namespace

const char* restartPublicationStatusName(
    RestartPublicationStatus status) noexcept
{
    switch(status)
    {
        case RestartPublicationStatus::SUCCESS: return "SUCCESS";
        case RestartPublicationStatus::INVALID_REQUEST:
            return "INVALID_REQUEST";
        case RestartPublicationStatus::INVALID_STORE_RESULT:
            return "INVALID_STORE_RESULT";
        case RestartPublicationStatus::INVALID_BASE_ROOT:
            return "INVALID_BASE_ROOT";
        case RestartPublicationStatus::LEDGER_REJECTED:
            return "LEDGER_REJECTED";
        case RestartPublicationStatus::LEDGER_CONFLICT:
            return "LEDGER_CONFLICT";
        case RestartPublicationStatus::MANIFEST_REJECTED:
            return "MANIFEST_REJECTED";
        case RestartPublicationStatus::CLOCK_REJECTED:
            return "CLOCK_REJECTED";
        case RestartPublicationStatus::OWNER_STAGER_MISSING:
            return "OWNER_STAGER_MISSING";
        case RestartPublicationStatus::OWNER_STAGER_UNEXPECTED:
            return "OWNER_STAGER_UNEXPECTED";
        case RestartPublicationStatus::OWNER_STAGE_REJECTED:
            return "OWNER_STAGE_REJECTED";
        case RestartPublicationStatus::OWNER_STAGE_EXCEPTION:
            return "OWNER_STAGE_EXCEPTION";
        case RestartPublicationStatus::OWNER_BINDING_MISMATCH:
            return "OWNER_BINDING_MISMATCH";
        case RestartPublicationStatus::OWNER_LEDGER_MISMATCH:
            return "OWNER_LEDGER_MISMATCH";
        case RestartPublicationStatus::CROSS_VALIDATOR_MISSING:
            return "CROSS_VALIDATOR_MISSING";
        case RestartPublicationStatus::CROSS_VALIDATION_REJECTED:
            return "CROSS_VALIDATION_REJECTED";
        case RestartPublicationStatus::CROSS_VALIDATION_EXCEPTION:
            return "CROSS_VALIDATION_EXCEPTION";
        case RestartPublicationStatus::STALE_PUBLICATION_BASE:
            return "STALE_PUBLICATION_BASE";
        case RestartPublicationStatus::PUBLICATION_ALREADY_CONSUMED:
            return "PUBLICATION_ALREADY_CONSUMED";
        case RestartPublicationStatus::INTERNAL_ERROR:
            return "INTERNAL_ERROR";
    }
    return "INTERNAL_ERROR";
}

bool validRestartGenerationLedgerEntry(
    const RestartGenerationLedgerEntry& entry) noexcept
{
    if(!RestartGenerationStore::validGenerationId(entry.generationId)
       ||!validRestartSha256(entry.generationChecksumSha256)
       ||!validPublicationStatus(entry.status)
       ||(entry.owner.has_value()&&!validOwner(*entry.owner)))
        return false;
    if(entry.disposition==RestartGenerationLedgerDisposition::ACCEPTED)
        return entry.status==RestartPublicationStatus::SUCCESS
            &&!entry.owner.has_value()&&entry.detail.empty();
    if(entry.disposition==RestartGenerationLedgerDisposition::REJECTED)
        return entry.status!=RestartPublicationStatus::SUCCESS
            &&!entry.detail.empty()
            &&entry.detail.size()<=maximumLedgerDetailBytes;
    return false;
}

bool validRestartAcceptedGenerationLedger(
    const std::vector<RestartGenerationLedgerEntry>& entries) noexcept
{
    for(std::size_t i=0;i<entries.size();++i)
    {
        if(!validRestartGenerationLedgerEntry(entries[i])
           ||!entries[i].accepted())
            return false;
        if(i>0U&&entries[i].ordinal<=entries[i-1U].ordinal) return false;
        for(std::size_t j=0;j<i;++j)
            if(entries[j].generationId==entries[i].generationId
               ||entries[j].ordinal==entries[i].ordinal)
                return false;
    }
    return true;
}

bool RestartRejectedGenerationLedger::record(
    const RestartGenerationLedgerEntry& entry) noexcept
{
    if(!entry.rejected()||!validRestartGenerationLedgerEntry(entry))
        return false;
    try
    {
        entries_.push_back(entry);
        return true;
    }
    catch(...)
    {
        return false;
    }
}

bool RestartPublishedOwner::valid() const noexcept
{
    if(!validOwner(record.owner)
       ||!RestartGenerationStore::validGenerationId(record.generationId)
       ||!validRestartSha256(record.identitySha256)
       ||!validRestartSha256(record.payloadSha256)
       ||stagedGenerationId!=record.generationId
       ||stagedIdentitySha256!=record.identitySha256
       ||stagedPayloadSha256!=record.payloadSha256)
        return false;
    if(record.disposition==RestartOwnerDisposition::STATE_PAYLOAD)
    {
        if(record.payloadBytes==0U||state==nullptr) return false;
        const bool isFluid=record.owner==RestartOwner::FLUID;
        return state->domainStepLedger().has_value()==isFluid;
    }
    if(record.disposition==RestartOwnerDisposition::DISABLED
       ||record.disposition==RestartOwnerDisposition::IDENTITY_ONLY)
        return record.payloadBytes==0U&&state==nullptr;
    return false;
}

bool RestartPublicationRoot::valid() const noexcept
{
    try
    {
        if(!RestartGenerationStore::validGenerationId(manifest_.generationId)
           ||!validRestartSha256(manifest_.generationChecksumSha256)
           ||manifest_.generationChecksumSha256
                !=restartManifestChecksum(manifest_)
           ||manifest_.owners.size()!=restartOwnerCount
           ||!clock_.valid()
           ||!sameClock(clock_,manifest_.generationState)
           ||!validRestartAcceptedGenerationLedger(acceptedLedger_)
           ||acceptedLedger_.empty())
            return false;

        const RestartGenerationLedgerEntry& latest=acceptedLedger_.back();
        if(latest.generationId!=manifest_.generationId
           ||latest.generationChecksumSha256
                !=manifest_.generationChecksumSha256)
            return false;

        const auto& roster=restartOwnerRoster();
        for(std::size_t i=0;i<owners_.size();++i)
            if(owners_[i].record.owner!=roster[i]
               ||owners_[i].record.publicationOrder
                    !=static_cast<std::uint32_t>(i)
               ||!sameOwnerRecord(owners_[i].record,manifest_.owners[i])
               ||!owners_[i].valid())
                return false;
        if(owners_[static_cast<std::size_t>(RestartOwner::FLUID)]
                .record.disposition!=RestartOwnerDisposition::STATE_PAYLOAD)
            return false;
        const std::optional<RestartDomainStepLedger>& domainLedger=
            owners_[static_cast<std::size_t>(RestartOwner::FLUID)]
                .state->domainStepLedger();
        if(!domainLedger.has_value()
           ||!domainLedger->consistentWith(manifest_.generationState))
            return false;
        return true;
    }
    catch(...)
    {
        return false;
    }
}

RestartStagedPublication::RestartStagedPublication(
    std::shared_ptr<const RestartPublicationRoot> expectedBase,
    std::shared_ptr<const RestartPublicationRoot> candidate,
    RestartGenerationLedgerEntry ledgerEntry)
    : expectedBase_(std::move(expectedBase)),candidate_(std::move(candidate)),
      ledgerEntry_(std::move(ledgerEntry))
{}

RestartPublicationCommitResult RestartStagedPublication::publish(
    std::shared_ptr<const RestartPublicationRoot>& liveRoot) noexcept
{
    try
    {
        if(consumed_)
            return commitFailure(
                RestartPublicationStatus::PUBLICATION_ALREADY_CONSUMED,
                "staged publication was already consumed",&ledgerEntry_);
        if(!candidate_||!candidate_->valid())
            return commitFailure(RestartPublicationStatus::INTERNAL_ERROR,
                "staged publication root is invalid",&ledgerEntry_);
        if(liveRoot.get()!=expectedBase_.get())
            return commitFailure(
                RestartPublicationStatus::STALE_PUBLICATION_BASE,
                "live publication root changed after staging",&ledgerEntry_);

        RestartPublicationCommitResult result;
        result.status=RestartPublicationStatus::SUCCESS;
        result.detail="detached owner root and accepted ledger published";

        // Build every potentially allocating return value before the
        // no-throw linearization point. A failure return must never coexist
        // with an already changed live root.
        std::shared_ptr<const RestartPublicationRoot> published=candidate_;
        liveRoot.swap(candidate_);
        candidate_.swap(published);
        consumed_=true;
        return result;
    }
    catch(...)
    {
        return RestartPublicationCommitResult{};
    }
}

RestartPublicationStageResult RestartPublicationStager::stage(
    const std::shared_ptr<const RestartPublicationRoot>& currentRoot,
    const RestartGenerationStoreLoadResult& load,
    const RestartManifestExpectation& expected,
    const RestartCheckpointClockObservation& checkpoint,
    const RestartPublicationRequest& request) noexcept
{
    try
    {
        if(!validRestartGenerationLedgerEntry(request.ledgerEntry))
            return stageFailure(RestartPublicationStatus::INVALID_REQUEST,
                "invalid generation ledger entry");
        if(!request.ledgerEntry.accepted())
            return stageFailure(RestartPublicationStatus::LEDGER_REJECTED,
                "rejected generations are not selectable",
                &request.ledgerEntry);
        if(request.crossValidate==nullptr)
            return stageFailure(
                RestartPublicationStatus::CROSS_VALIDATOR_MISSING,
                "a complete-root cross validator is required",
                &request.ledgerEntry);

        const auto& roster=restartOwnerRoster();
        for(std::size_t i=0;i<request.ownerStagers.size();++i)
            if(request.ownerStagers[i].owner!=roster[i])
                return stageFailure(RestartPublicationStatus::INVALID_REQUEST,
                    "owner stager roster is not canonical",
                    &request.ledgerEntry,
                    validOwner(request.ownerStagers[i].owner)
                        ?std::optional<RestartOwner>(request.ownerStagers[i].owner)
                        :std::nullopt);

        if(currentRoot&&!currentRoot->valid())
            return stageFailure(RestartPublicationStatus::INVALID_BASE_ROOT,
                "current publication root is invalid",
                &request.ledgerEntry);
        const std::vector<RestartGenerationLedgerEntry> emptyLedger;
        const auto& accepted=currentRoot
            ?currentRoot->acceptedLedger():emptyLedger;
        if(!ledgerAllows(accepted,request.ledgerEntry))
            return stageFailure(RestartPublicationStatus::LEDGER_CONFLICT,
                "accepted generation or ordinal conflicts with the live ledger",
                &request.ledgerEntry);

        if(!load.success()||!load.candidate.has_value())
            return stageFailure(
                RestartPublicationStatus::INVALID_STORE_RESULT,
                "generation store did not expose a committed candidate",
                &request.ledgerEntry);
        const RestartCandidatePackage& candidate=*load.candidate;
        if(candidate.manifest.generationId!=request.ledgerEntry.generationId
           ||candidate.manifest.generationChecksumSha256
                !=request.ledgerEntry.generationChecksumSha256)
            return stageFailure(RestartPublicationStatus::LEDGER_REJECTED,
                "ledger identity does not match the loaded generation",
                &request.ledgerEntry);

        const RestartManifestValidationResult manifestValidation=
            RestartManifestValidator::validate(candidate,expected);
        if(!manifestValidation.success())
            return stageFailure(RestartPublicationStatus::MANIFEST_REJECTED,
                std::string("manifest revalidation failed: ")
                    +restartManifestStatusName(manifestValidation.status)
                    +": "+manifestValidation.detail,
                &request.ledgerEntry);

        if(candidate.manifest.owners[
                static_cast<std::size_t>(RestartOwner::FLUID)].disposition
            !=RestartOwnerDisposition::STATE_PAYLOAD)
            return stageFailure(RestartPublicationStatus::MANIFEST_REJECTED,
                "composite publication requires a fluid Domain payload",
                &request.ledgerEntry,RestartOwner::FLUID);

        const RestartClockAdmissionResult clock=RestartClockAdmission::stage(
            load,checkpoint,expected.requirePreviousAcceptedDt);
        if(!clock.success())
            return stageFailure(RestartPublicationStatus::CLOCK_REJECTED,
                std::string("clock admission failed: ")
                    +restartClockAdmissionStatusName(clock.status),
                &request.ledgerEntry);

        std::shared_ptr<RestartPublicationRoot> root(
            new RestartPublicationRoot());
        root->manifest_=candidate.manifest;
        root->clock_=*clock.publication;
        root->acceptedLedger_=accepted;
        root->acceptedLedger_.push_back(request.ledgerEntry);

        for(std::size_t i=0;i<roster.size();++i)
        {
            const RestartOwnerRecord& owner=candidate.manifest.owners[i];
            const RestartOwnerStager& stager=request.ownerStagers[i];
            RestartPublishedOwner& published=root->owners_[i];
            published.record=owner;
            published.stagedGenerationId=owner.generationId;
            published.stagedIdentitySha256=owner.identitySha256;
            published.stagedPayloadSha256=owner.payloadSha256;

            if(owner.disposition!=RestartOwnerDisposition::STATE_PAYLOAD)
            {
                if(stager.stage!=nullptr)
                    return stageFailure(
                        RestartPublicationStatus::OWNER_STAGER_UNEXPECTED,
                        "non-state owner supplied a state stager",
                        &request.ledgerEntry,owner.owner);
                continue;
            }
            if(stager.stage==nullptr)
                return stageFailure(
                    RestartPublicationStatus::OWNER_STAGER_MISSING,
                    "state-payload owner has no detached stager",
                    &request.ledgerEntry,owner.owner);
            const RestartPayloadCandidate* payload=
                findPayload(candidate,owner.owner);
            if(!payload)
                return stageFailure(
                    RestartPublicationStatus::MANIFEST_REJECTED,
                    "validated package lost a state-owner payload",
                    &request.ledgerEntry,owner.owner);

            RestartOwnerStageResult staged;
            try
            {
                const RestartOwnerStageInput input{
                    candidate.manifest,owner,*payload};
                staged=stager.stage(input,stager.context);
            }
            catch(const std::exception& error)
            {
                return stageFailure(
                    RestartPublicationStatus::OWNER_STAGE_EXCEPTION,
                    std::string("owner stager threw: ")+error.what(),
                    &request.ledgerEntry,owner.owner);
            }
            catch(...)
            {
                return stageFailure(
                    RestartPublicationStatus::OWNER_STAGE_EXCEPTION,
                    "owner stager threw a non-standard exception",
                    &request.ledgerEntry,owner.owner);
            }
            if(!staged.success())
                return stageFailure(
                    RestartPublicationStatus::OWNER_STAGE_REJECTED,
                    boundedDetail(staged.detail,"owner staging rejected payload"),
                    &request.ledgerEntry,owner.owner);
            if(staged.owner!=owner.owner
               ||staged.generationId!=owner.generationId
               ||staged.identitySha256!=owner.identitySha256
               ||staged.payloadSha256!=owner.payloadSha256)
                return stageFailure(
                    RestartPublicationStatus::OWNER_BINDING_MISMATCH,
                    "detached owner state is not bound to its manifest record",
                    &request.ledgerEntry,owner.owner);
            const bool isFluid=owner.owner==RestartOwner::FLUID;
            const std::optional<RestartDomainStepLedger>& domainLedger=
                staged.state->domainStepLedger();
            if(domainLedger.has_value()!=isFluid
               ||(isFluid&&!domainLedger->consistentWith(
                    candidate.manifest.generationState)))
                return stageFailure(
                    RestartPublicationStatus::OWNER_LEDGER_MISMATCH,
                    isFluid
                        ?"fluid Domain accepted-step ledger does not match the coordinator clock"
                        :"non-fluid owner supplied a Domain step ledger",
                    &request.ledgerEntry,owner.owner);
            published.state=std::move(staged.state);
        }

        if(!root->valid())
            return stageFailure(RestartPublicationStatus::INTERNAL_ERROR,
                "assembled publication root failed invariant validation",
                &request.ledgerEntry);

        RestartCrossValidationResult cross;
        try
        {
            cross=request.crossValidate(
                *root,request.crossValidationContext);
        }
        catch(const std::exception& error)
        {
            return stageFailure(
                RestartPublicationStatus::CROSS_VALIDATION_EXCEPTION,
                std::string("cross validator threw: ")+error.what(),
                &request.ledgerEntry);
        }
        catch(...)
        {
            return stageFailure(
                RestartPublicationStatus::CROSS_VALIDATION_EXCEPTION,
                "cross validator threw a non-standard exception",
                &request.ledgerEntry);
        }
        if(!cross.accepted)
            return stageFailure(
                RestartPublicationStatus::CROSS_VALIDATION_REJECTED,
                boundedDetail(cross.detail,
                    "complete publication root was rejected"),
                &request.ledgerEntry);

        std::shared_ptr<const RestartPublicationRoot> immutableRoot=root;
        RestartPublicationStageResult result;
        result.status=RestartPublicationStatus::SUCCESS;
        result.detail="all state owners and clock staged in one detached root";
        result.publication=std::unique_ptr<RestartStagedPublication>(
            new RestartStagedPublication(
                currentRoot,std::move(immutableRoot),request.ledgerEntry));
        return result;
    }
    catch(const std::exception& error)
    {
        try
        {
            return stageFailure(RestartPublicationStatus::INTERNAL_ERROR,
                std::string("unexpected staging exception: ")+error.what(),
                &request.ledgerEntry);
        }
        catch(...)
        {
            return RestartPublicationStageResult{};
        }
    }
    catch(...)
    {
        return RestartPublicationStageResult{};
    }
}

} // End namespace Foam
