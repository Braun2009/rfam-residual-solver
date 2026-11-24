#ifndef RestartManifest_H
#define RestartManifest_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace Foam
{

enum class RestartOwner : std::uint8_t
{
    FLUID = 0,
    TURBULENCE,
    CHT,
    CONDENSED_CARBON,
    CHEMISTRY,
    TOPOLOGY_ENVELOPE,
    TRANSPORT,
    COUNT
};

constexpr std::size_t restartOwnerCount =
    static_cast<std::size_t>(RestartOwner::COUNT);

const char* restartOwnerName(RestartOwner owner);
const std::array<RestartOwner,restartOwnerCount>& restartOwnerRoster();

enum class RestartOwnerDisposition : std::uint8_t
{
    DISABLED = 0,
    IDENTITY_ONLY,
    STATE_PAYLOAD
};

const char* restartOwnerDispositionName(RestartOwnerDisposition disposition);

enum class RestartExecutionScope : std::uint8_t
{
    SINGLE_PROCESS = 0,
    MPI_PARTITIONED
};

enum class RestartByteOrder : std::uint8_t
{
    LITTLE_ENDIAN_IEEE754 = 0,
    BIG_ENDIAN_IEEE754
};

RestartByteOrder nativeRestartByteOrder();

enum class RestartPhaseAssertion : std::uint8_t
{
    UNKNOWN = 0,
    SINGLE_LIQUID,
    SINGLE_VAPOR,
    SUPERCRITICAL,
    TWO_PHASE,
    THREE_PHASE_REQUIRED
};

const char* restartPhaseAssertionName(RestartPhaseAssertion assertion);

struct RestartIdentitySet
{
    std::string configurationSha256;
    std::string buildSha256;
    std::string meshSha256;
    std::string thermoModelId;
    std::string kijManifestSha256;
    std::string mechanismSha256;
    std::string chemistryIntegratorSha256;
    std::string transportSha256;
    std::string topologyEnvelopeSha256;
};

struct RestartOwnerRecord
{
    RestartOwner owner = RestartOwner::FLUID;
    RestartOwnerDisposition disposition = RestartOwnerDisposition::DISABLED;
    std::uint32_t publicationOrder = 0;
    std::string generationId;
    std::string identitySha256;
    std::uint64_t payloadBytes = 0;
    std::string payloadSha256;
};

struct RestartGenerationState
{
    double physicalTime = std::numeric_limits<double>::quiet_NaN();
    std::uint64_t acceptedIteration = 0;
    bool previousAcceptedDtPresent = false;
    double previousAcceptedDt = 0.0;
};

// The fingerprint covers the complete canonical rho-e-Y recovery, including
// region, quality, phase presence/composition/caloric data, mixture caloric
// data, acoustic closure, and thermodynamic derivatives. The serialized
// assertion is never an alternative authority to rho-e-Y recovery.
struct RestartRecoverySeed
{
    std::uint64_t cell = 0;
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double temperature = std::numeric_limits<double>::quiet_NaN();
    RestartPhaseAssertion expectedPhase = RestartPhaseAssertion::UNKNOWN;
    std::string phaseIdentitySha256;
};

struct RestartManifest
{
    static constexpr std::uint32_t currentVersion = 2;
    static constexpr const char* schemaId = "rfam.restart-manifest/v2";
    static constexpr const char* floatingPointFormatId =
        "ieee754-binary64";
    static constexpr const char* singleProcessDecompositionId =
        "single-process-global-domain/v1";

    std::string schema = schemaId;
    std::uint32_t version = currentVersion;
    std::string generationId;
    std::string generationChecksumSha256;
    RestartGenerationState generationState;
    RestartExecutionScope executionScope =
        RestartExecutionScope::SINGLE_PROCESS;
    std::uint32_t rankCount = 1;
    std::uint32_t rankIndex = 0;
    bool collectiveCommit = false;
    std::string decompositionId = singleProcessDecompositionId;
    RestartByteOrder byteOrder = RestartByteOrder::LITTLE_ENDIAN_IEEE754;
    std::string floatingPointFormat = floatingPointFormatId;
    bool trailingDataPermitted = false;
    std::uint64_t cellCount = 0;
    RestartIdentitySet identities;
    std::vector<RestartOwnerRecord> owners;
    std::vector<RestartRecoverySeed> recoverySeeds;
};

struct RestartPayloadCandidate
{
    RestartOwner owner = RestartOwner::FLUID;
    std::string generationId;
    std::string bytes;
    std::uint64_t consumedBytes = 0;
};

// A strict decoder supplies manifestTrailingBytes and each consumedBytes value.
// The validator rejects any undecoded suffix before it exposes publicationOrder.
struct RestartCandidatePackage
{
    RestartManifest manifest;
    std::uint64_t manifestTrailingBytes = 0;
    std::vector<RestartPayloadCandidate> payloads;
};

struct RestartOwnerExpectation
{
    RestartOwner owner = RestartOwner::FLUID;
    RestartOwnerDisposition disposition = RestartOwnerDisposition::DISABLED;
    std::string identitySha256;
};

struct RestartManifestExpectation
{
    RestartIdentitySet identities;
    RestartByteOrder byteOrder = RestartByteOrder::LITTLE_ENDIAN_IEEE754;
    std::string floatingPointFormat = RestartManifest::floatingPointFormatId;
    std::string decompositionId =
        RestartManifest::singleProcessDecompositionId;
    std::uint64_t cellCount = 0;
    bool requirePreviousAcceptedDt = false;
    std::uint64_t maximumAcceptedIteration =
        static_cast<std::uint64_t>(std::numeric_limits<int>::max());
    std::array<RestartOwnerExpectation,restartOwnerCount> owners{};
};

enum class RestartManifestStatus
{
    SUCCESS,
    INVALID_EXPECTATION,
    INVALID_SCHEMA,
    UNSUPPORTED_VERSION,
    INVALID_GENERATION,
    INVALID_GENERATION_STATE,
    UNSUPPORTED_DISTRIBUTED_RESTART,
    UNSUPPORTED_ENCODING,
    TRAILING_DATA,
    IDENTITY_MISMATCH,
    INVALID_OWNER_ROSTER,
    OWNER_DISPOSITION_MISMATCH,
    OWNER_GENERATION_MISMATCH,
    OWNER_IDENTITY_MISMATCH,
    PAYLOAD_MISSING,
    PAYLOAD_UNEXPECTED,
    PAYLOAD_GENERATION_MISMATCH,
    PAYLOAD_SIZE_MISMATCH,
    PAYLOAD_CHECKSUM_MISMATCH,
    INVALID_RECOVERY_SEED,
    UNSUPPORTED_PHASE_TOPOLOGY,
    GENERATION_CHECKSUM_MISMATCH
};

const char* restartManifestStatusName(RestartManifestStatus status);

struct RestartManifestValidationResult
{
    RestartManifestStatus status = RestartManifestStatus::INVALID_SCHEMA;
    std::string detail;
    std::vector<RestartOwner> publicationOrder;
    std::optional<RestartGenerationState> generationState;

    bool success() const
    {
        return status == RestartManifestStatus::SUCCESS
            && publicationOrder.size() == restartOwnerCount
            && generationState.has_value();
    }
};

bool validRestartSha256(const std::string& value);
std::string restartManifestChecksum(const RestartManifest& manifest);

class RestartManifestValidator
{
public:
    static RestartManifestValidationResult validate(
        const RestartCandidatePackage& candidate,
        const RestartManifestExpectation& expected);
};

} // End namespace Foam

#endif
