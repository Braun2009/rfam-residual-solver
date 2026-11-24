#include "RestartLiveRoot.H"

#include <cmath>
#include <limits>
#include <mutex>
#include <utility>

namespace Foam
{
namespace
{

bool clockBoundToManifest(const RestartPublicationRoot& root) noexcept
{
    const RestartGenerationState& generation=
        root.manifest().generationState;
    const RestartClockPublication& clock=root.clock();
    if(!clock.valid()
       ||!std::isfinite(generation.physicalTime)
       ||generation.physicalTime<0.0
       ||std::signbit(generation.physicalTime)
       ||generation.acceptedIteration>
            static_cast<std::uint64_t>(std::numeric_limits<int>::max())
       ||clock.physicalTime!=generation.physicalTime
       ||clock.acceptedIteration
            !=static_cast<int>(generation.acceptedIteration)
       ||clock.previousAcceptedDt.has_value()
            !=generation.previousAcceptedDtPresent)
        return false;
    if(generation.previousAcceptedDtPresent)
        return std::isfinite(generation.previousAcceptedDt)
            &&generation.previousAcceptedDt>0.0
            &&*clock.previousAcceptedDt==generation.previousAcceptedDt;
    return generation.previousAcceptedDt==0.0
        &&!std::signbit(generation.previousAcceptedDt);
}

RestartLiveRootPublicationResult failure(
    RestartLiveRootStatus status,
    const std::string& detail)
{
    RestartLiveRootPublicationResult result;
    result.status=status;
    result.detail=detail;
    return result;
}

} // namespace

const char* restartLiveRootStatusName(
    RestartLiveRootStatus status) noexcept
{
    switch(status)
    {
        case RestartLiveRootStatus::SUCCESS: return "SUCCESS";
        case RestartLiveRootStatus::INVALID_REQUEST:
            return "INVALID_REQUEST";
        case RestartLiveRootStatus::READERS_OR_WRITER_ACTIVE:
            return "READERS_OR_WRITER_ACTIVE";
        case RestartLiveRootStatus::INVALID_LIVE_ROOT:
            return "INVALID_LIVE_ROOT";
        case RestartLiveRootStatus::QUIESCENCE_REJECTED:
            return "QUIESCENCE_REJECTED";
        case RestartLiveRootStatus::QUIESCENCE_EXCEPTION:
            return "QUIESCENCE_EXCEPTION";
        case RestartLiveRootStatus::PREPUBLICATION_REJECTED:
            return "PREPUBLICATION_REJECTED";
        case RestartLiveRootStatus::PREPUBLICATION_EXCEPTION:
            return "PREPUBLICATION_EXCEPTION";
        case RestartLiveRootStatus::STAGED_PUBLICATION_REJECTED:
            return "STAGED_PUBLICATION_REJECTED";
        case RestartLiveRootStatus::TEST_FAULT_ROLLED_BACK:
            return "TEST_FAULT_ROLLED_BACK";
        case RestartLiveRootStatus::INTERNAL_ERROR:
            return "INTERNAL_ERROR";
    }
    return "INTERNAL_ERROR";
}

RestartLiveRootReadLease::RestartLiveRootReadLease(
    std::shared_lock<std::shared_mutex>&& lock,
    std::shared_ptr<const RestartPublicationRoot> root) noexcept
    : lock_(std::move(lock)),root_(std::move(root))
{}

RestartLiveRootReadLease RestartLiveRoot::read() const noexcept
{
    try
    {
        std::shared_lock<std::shared_mutex> lock(mutex_);
        return RestartLiveRootReadLease(std::move(lock),root_);
    }
    catch(...)
    {
        return RestartLiveRootReadLease{};
    }
}

std::shared_ptr<const RestartPublicationRoot>
RestartLiveRoot::stagingBase() const noexcept
{
    try
    {
        const std::shared_lock<std::shared_mutex> lock(mutex_);
        return root_;
    }
    catch(...)
    {
        return {};
    }
}

RestartLiveRootPublicationResult RestartLiveRoot::tryPublish(
    RestartStagedPublication& publication,
    const RestartLiveRootPublicationRequest& request) noexcept
{
    try
    {
        if(!request.observeQuiescence||!request.prepublishValidate
           ||request.fault>RestartLiveRootFaultPoint::AFTER_ROOT_SWAP)
            return failure(RestartLiveRootStatus::INVALID_REQUEST,
                "quiescence and prepublication validators are required");

        std::unique_lock<std::shared_mutex> lock(mutex_,std::try_to_lock);
        if(!lock.owns_lock())
            return failure(
                RestartLiveRootStatus::READERS_OR_WRITER_ACTIVE,
                "a routed reader or another publisher is active");
        if(root_&&!root_->valid())
            return failure(RestartLiveRootStatus::INVALID_LIVE_ROOT,
                "current routed root is invalid");

        RestartRuntimeQuiescenceResult quiescence;
        try
        {
            quiescence=request.observeQuiescence(
                request.quiescenceContext);
        }
        catch(const std::exception& error)
        {
            return failure(RestartLiveRootStatus::QUIESCENCE_EXCEPTION,
                std::string("quiescence probe threw: ")+error.what());
        }
        catch(...)
        {
            return failure(RestartLiveRootStatus::QUIESCENCE_EXCEPTION,
                "quiescence probe threw a non-standard exception");
        }
        if(!quiescence.accepted)
            return failure(RestartLiveRootStatus::QUIESCENCE_REJECTED,
                quiescence.detail.empty()
                    ?"runtime quiescence was not established"
                    :quiescence.detail);

        const RestartPublicationRoot& candidate=publication.candidate();
        if(!candidate.valid()||!clockBoundToManifest(candidate))
            return failure(RestartLiveRootStatus::PREPUBLICATION_REJECTED,
                "candidate clock is not exactly bound to its manifest");
        RestartCrossValidationResult validation;
        try
        {
            validation=request.prepublishValidate(
                candidate,request.prepublishValidationContext);
        }
        catch(const std::exception& error)
        {
            return failure(
                RestartLiveRootStatus::PREPUBLICATION_EXCEPTION,
                std::string("prepublication validator threw: ")
                    +error.what());
        }
        catch(...)
        {
            return failure(
                RestartLiveRootStatus::PREPUBLICATION_EXCEPTION,
                "prepublication validator threw a non-standard exception");
        }
        if(!validation.accepted)
            return failure(RestartLiveRootStatus::PREPUBLICATION_REJECTED,
                validation.detail.empty()
                    ?"candidate root failed prepublication validation"
                    :validation.detail);

        std::shared_ptr<const RestartPublicationRoot> next=root_;
        RestartPublicationCommitResult stagedCommit=
            publication.publish(next);
        if(!stagedCommit.success())
        {
            RestartLiveRootPublicationResult result=failure(
                RestartLiveRootStatus::STAGED_PUBLICATION_REJECTED,
                stagedCommit.detail);
            result.stagedCommit=std::move(stagedCommit);
            return result;
        }
        if(!next||!next->valid()||!clockBoundToManifest(*next))
            return failure(RestartLiveRootStatus::INTERNAL_ERROR,
                "staged publication returned an invalid root");

        if(request.fault==RestartLiveRootFaultPoint::AFTER_ROOT_SWAP)
        {
            RestartLiveRootPublicationResult result=failure(
                RestartLiveRootStatus::TEST_FAULT_ROLLED_BACK,
                "injected post-swap fault rolled back under exclusive lock");
            result.stagedCommit=std::move(stagedCommit);
            result.rollbackCompleted=true;

            root_.swap(next);
            root_.swap(next);
            return result;
        }

        RestartLiveRootPublicationResult result;
        result.status=RestartLiveRootStatus::SUCCESS;
        result.detail="one immutable root published at quiescence";
        result.stagedCommit=std::move(stagedCommit);

        // Everything that can allocate is complete before this one no-throw
        // root swap. Routed readers cannot enter until the unique lock leaves.
        root_.swap(next);
        return result;
    }
    catch(const std::exception& error)
    {
        try
        {
            return failure(RestartLiveRootStatus::INTERNAL_ERROR,
                std::string("unexpected live-root exception: ")
                    +error.what());
        }
        catch(...)
        {
            return RestartLiveRootPublicationResult{};
        }
    }
    catch(...)
    {
        return RestartLiveRootPublicationResult{};
    }
}

} // End namespace Foam
