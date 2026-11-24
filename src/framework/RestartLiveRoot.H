#ifndef RestartLiveRoot_H
#define RestartLiveRoot_H

#include "RestartPublication.H"

#include <memory>
#include <shared_mutex>
#include <string>

namespace Foam
{

enum class RestartLiveRootStatus : std::uint8_t
{
    SUCCESS = 0,
    INVALID_REQUEST,
    READERS_OR_WRITER_ACTIVE,
    INVALID_LIVE_ROOT,
    QUIESCENCE_REJECTED,
    QUIESCENCE_EXCEPTION,
    PREPUBLICATION_REJECTED,
    PREPUBLICATION_EXCEPTION,
    STAGED_PUBLICATION_REJECTED,
    TEST_FAULT_ROLLED_BACK,
    INTERNAL_ERROR
};

const char* restartLiveRootStatusName(
    RestartLiveRootStatus status) noexcept;

struct RestartRuntimeQuiescenceResult
{
    bool accepted = false;
    std::string detail;
};

using RestartRuntimeQuiescenceFunction =
    RestartRuntimeQuiescenceResult (*)(void* context);

// Test-only fault point. The swap-back occurs while the exclusive live-root
// lock is retained, so no routed reader can observe the candidate root.
enum class RestartLiveRootFaultPoint : std::uint8_t
{
    NONE = 0,
    AFTER_ROOT_SWAP
};

struct RestartLiveRootPublicationRequest
{
    RestartRuntimeQuiescenceFunction observeQuiescence = nullptr;
    void* quiescenceContext = nullptr;
    RestartCrossValidationFunction prepublishValidate = nullptr;
    void* prepublishValidationContext = nullptr;
    RestartLiveRootFaultPoint fault = RestartLiveRootFaultPoint::NONE;
};

struct RestartLiveRootPublicationResult
{
    RestartLiveRootStatus status = RestartLiveRootStatus::INTERNAL_ERROR;
    std::string detail;
    RestartPublicationCommitResult stagedCommit;
    bool rollbackCompleted = false;

    bool success() const noexcept
    {
        return status == RestartLiveRootStatus::SUCCESS
            && stagedCommit.success() && !rollbackCompleted;
    }
};

class RestartLiveRoot;

// A read lease is the only runtime-facing access to the current immutable
// root. Its shared lock excludes publication until every routed reader leaves.
class RestartLiveRootReadLease
{
    std::shared_lock<std::shared_mutex> lock_;
    std::shared_ptr<const RestartPublicationRoot> root_;

    RestartLiveRootReadLease(
        std::shared_lock<std::shared_mutex>&& lock,
        std::shared_ptr<const RestartPublicationRoot> root) noexcept;

    friend class RestartLiveRoot;

public:
    RestartLiveRootReadLease() = default;
    RestartLiveRootReadLease(RestartLiveRootReadLease&&) noexcept = default;
    RestartLiveRootReadLease& operator=(
        RestartLiveRootReadLease&&) noexcept = default;
    RestartLiveRootReadLease(const RestartLiveRootReadLease&) = delete;
    RestartLiveRootReadLease& operator=(
        const RestartLiveRootReadLease&) = delete;

    bool available() const noexcept { return root_ != nullptr; }
    const RestartPublicationRoot* get() const noexcept { return root_.get(); }
    const RestartPublicationRoot& operator*() const { return *root_; }
    const RestartPublicationRoot* operator->() const noexcept
    { return root_.get(); }
};

// Candidate-only single-process indirection. Staging occurs outside this
// class. tryPublish() acquires exclusive routed-reader quiescence, rechecks the
// complete candidate, and performs one noexcept shared_ptr swap.
class RestartLiveRoot
{
    mutable std::shared_mutex mutex_;
    std::shared_ptr<const RestartPublicationRoot> root_;

public:
    RestartLiveRoot() = default;
    RestartLiveRoot(const RestartLiveRoot&) = delete;
    RestartLiveRoot& operator=(const RestartLiveRoot&) = delete;

    RestartLiveRootReadLease read() const noexcept;

    // A short-lived snapshot for detached staging. It does not retain a read
    // lease; stale-base admission is rechecked by the staged publication.
    std::shared_ptr<const RestartPublicationRoot> stagingBase() const noexcept;

    RestartLiveRootPublicationResult tryPublish(
        RestartStagedPublication& publication,
        const RestartLiveRootPublicationRequest& request) noexcept;
};

} // End namespace Foam

#endif
