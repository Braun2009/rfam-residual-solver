#ifndef RestartGenerationStore_H
#define RestartGenerationStore_H

#include "RestartManifestCodec.H"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>

namespace Foam
{

enum class RestartGenerationStoreStatus : std::uint8_t
{
    SUCCESS = 0,
    INVALID_INPUT,
    UNSAFE_GENERATION_ID,
    STORE_ROOT_UNAVAILABLE,
    GENERATION_CONFLICT,
    INCOMPLETE_GENERATION,
    CODEC_REJECTED,
    STAGE_OPEN_FAILED,
    STAGE_WRITE_FAILED,
    STAGE_FLUSH_FAILED,
    STAGE_CLOSE_FAILED,
    FILE_SYNC_FAILED,
    STAGE_READBACK_FAILED,
    ATOMIC_RENAME_FAILED,
    DIRECTORY_SYNC_FAILED,
    NOT_COMMITTED,
    MANIFEST_READ_FAILED,
    MANIFEST_REJECTED,
    GENERATION_MISMATCH,
    OWNER_FILE_MISSING,
    OWNER_FILE_UNEXPECTED,
    OWNER_GENERATION_MISMATCH,
    OWNER_SIZE_MISMATCH,
    OWNER_CHECKSUM_MISMATCH,
    OWNER_CONTENT_MISMATCH,
    FILE_LAYOUT_REJECTED,
    PACKAGE_REJECTED,
    STORE_ROOT_CHANGED,
    COMMIT_STATE_UNCERTAIN,
    INTERNAL_ERROR
};

const char* restartGenerationStoreStatusName(
    RestartGenerationStoreStatus status);

enum class RestartGenerationStoreFaultPoint : std::uint8_t
{
    NONE = 0,
    PAYLOAD_OPEN,
    PAYLOAD_WRITE,
    PAYLOAD_FLUSH,
    PAYLOAD_CLOSE,
    PAYLOAD_FILE_SYNC,
    PAYLOAD_READBACK,
    PAYLOAD_RENAME,
    PAYLOAD_DIRECTORY_SYNC,
    MANIFEST_OPEN,
    MANIFEST_WRITE,
    MANIFEST_FLUSH,
    MANIFEST_CLOSE,
    MANIFEST_FILE_SYNC,
    MANIFEST_READBACK,
    MANIFEST_RENAME,
    MANIFEST_DIRECTORY_SYNC,
    MANIFEST_ROLLBACK_REMOVE,
    MANIFEST_ROLLBACK_DIRECTORY_SYNC
};

// Runs immediately before the final retained root/generation identity check,
// including after a durably synchronized manifest-marker rollback.
using RestartGenerationStoreFinalIdentityHook = void (*)(void*) noexcept;

// Fault injection is a test-only boundary. Payload fault points apply to the
// selected owner; manifest fault points ignore owner.
struct RestartGenerationStoreFault
{
    RestartGenerationStoreFaultPoint point =
        RestartGenerationStoreFaultPoint::NONE;
    RestartOwner owner = RestartOwner::FLUID;
    RestartGenerationStoreFinalIdentityHook beforeFinalIdentityCheck = nullptr;
    void* hookContext = nullptr;
};

struct RestartGenerationStorePublishResult
{
    RestartGenerationStoreStatus status =
        RestartGenerationStoreStatus::INVALID_INPUT;
    std::string detail;
    std::filesystem::path manifestPath;

    bool success() const
    {
        return status == RestartGenerationStoreStatus::SUCCESS
            && !manifestPath.empty();
    }
};

struct RestartGenerationStoreLoadResult
{
    RestartGenerationStoreStatus status =
        RestartGenerationStoreStatus::INVALID_INPUT;
    std::string detail;
    std::optional<RestartCandidatePackage> candidate;

    bool success() const
    {
        return status == RestartGenerationStoreStatus::SUCCESS
            && candidate.has_value();
    }
};

// Candidate-only, single-process generation storage. This class does not
// select a generation, restore an owner, mutate Domain, or authorize a run.
class RestartGenerationStore
{
public:
    static constexpr const char* manifestFileName = "manifest.rfamrm2";

    static bool validGenerationId(
        const std::string& generationId) noexcept;

    static std::filesystem::path generationDirectory(
        const std::filesystem::path& storeRoot,
        const std::string& generationId);

    static std::string ownerPayloadFileName(RestartOwner owner);

    static RestartGenerationStorePublishResult publish(
        const std::filesystem::path& storeRoot,
        const RestartCandidatePackage& candidate,
        const RestartManifestExpectation& expected,
        const RestartGenerationStoreFault& fault =
            RestartGenerationStoreFault{}) noexcept;

    static RestartGenerationStoreLoadResult load(
        const std::filesystem::path& storeRoot,
        const std::string& generationId,
        const RestartManifestExpectation& expected,
        const RestartGenerationStoreFault& fault =
            RestartGenerationStoreFault{}) noexcept;
};

} // End namespace Foam

#endif
