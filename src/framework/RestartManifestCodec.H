#ifndef RestartManifestCodec_H
#define RestartManifestCodec_H

#include "RestartManifest.H"

#include <cstdint>
#include <optional>
#include <string>

namespace Foam
{

enum class RestartCodecStatus : std::uint8_t
{
    SUCCESS = 0,
    INVALID_INPUT,
    INVALID_MAGIC,
    UNSUPPORTED_CODEC_VERSION,
    TRUNCATED_INPUT,
    TRAILING_DATA,
    DUPLICATE_FIELD,
    UNKNOWN_FIELD,
    MISSING_FIELD,
    NONCANONICAL_FIELD_ORDER,
    TYPE_MISMATCH,
    RANGE_ERROR,
    RESOURCE_LIMIT,
    MALFORMED_RECORD,
    MANIFEST_REJECTED
};

const char* restartCodecStatusName(RestartCodecStatus status);

struct RestartCodecEncodeResult
{
    RestartCodecStatus status = RestartCodecStatus::INVALID_INPUT;
    std::string detail;
    std::string bytes;

    bool success() const
    {
        return status == RestartCodecStatus::SUCCESS && !bytes.empty();
    }
};

struct RestartCodecDecodeResult
{
    RestartCodecStatus status = RestartCodecStatus::INVALID_INPUT;
    std::string detail;
    std::optional<RestartCandidatePackage> candidate;

    bool success() const
    {
        return status == RestartCodecStatus::SUCCESS
            && candidate.has_value();
    }
};

// Deterministic binary codec for the candidate-only restart package. The
// encoded bytes contain no publication instruction and perform no file I/O.
class RestartManifestCodec
{
public:
    static constexpr std::uint32_t diskFormatVersion = 1U;

    static RestartCodecEncodeResult encode(
        const RestartCandidatePackage& candidate,
        const RestartManifestExpectation& expected) noexcept;

    static RestartCodecDecodeResult decode(
        const std::string& bytes,
        const RestartManifestExpectation& expected) noexcept;
};

} // End namespace Foam

#endif
