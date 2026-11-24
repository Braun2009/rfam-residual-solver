#include "RestartManifest.H"

#include "ThermoDataProvenance.H"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Foam
{
namespace
{

using thermoProvenance::appendDouble;
using thermoProvenance::appendSize;
using thermoProvenance::appendString;

RestartManifestValidationResult failure(
    RestartManifestStatus status,
    const std::string& detail)
{
    RestartManifestValidationResult result;
    result.status = status;
    result.detail = detail;
    return result;
}

bool validToken(const std::string& value,std::size_t maximumLength)
{
    if(value.empty()||value.size()>maximumLength) return false;
    return std::all_of(value.begin(),value.end(),[](unsigned char c)
    {
        return (c>='a'&&c<='z')||(c>='A'&&c<='Z')
            ||(c>='0'&&c<='9')||c=='-'||c=='_'||c=='.'||c==':';
    });
}

bool validIdentityText(const std::string& value)
{
    if(value.empty()||value.size()>4096U) return false;
    return std::all_of(value.begin(),value.end(),[](unsigned char c)
    { return c>=0x21U&&c<=0x7eU; });
}

bool sameIdentities(
    const RestartIdentitySet& actual,
    const RestartIdentitySet& expected)
{
    return actual.configurationSha256==expected.configurationSha256
        &&actual.buildSha256==expected.buildSha256
        &&actual.meshSha256==expected.meshSha256
        &&actual.thermoModelId==expected.thermoModelId
        &&actual.kijManifestSha256==expected.kijManifestSha256
        &&actual.mechanismSha256==expected.mechanismSha256
        &&actual.chemistryIntegratorSha256
            ==expected.chemistryIntegratorSha256
        &&actual.transportSha256==expected.transportSha256
        &&actual.topologyEnvelopeSha256
            ==expected.topologyEnvelopeSha256;
}

bool validIdentities(const RestartIdentitySet& identities)
{
    return validRestartSha256(identities.configurationSha256)
        &&validRestartSha256(identities.buildSha256)
        &&validRestartSha256(identities.meshSha256)
        &&validIdentityText(identities.thermoModelId)
        &&validRestartSha256(identities.kijManifestSha256)
        &&validRestartSha256(identities.mechanismSha256)
        &&validRestartSha256(identities.chemistryIntegratorSha256)
        &&validRestartSha256(identities.transportSha256)
        &&validRestartSha256(identities.topologyEnvelopeSha256);
}

bool validDisposition(RestartOwnerDisposition disposition)
{
    return disposition==RestartOwnerDisposition::DISABLED
        ||disposition==RestartOwnerDisposition::IDENTITY_ONLY
        ||disposition==RestartOwnerDisposition::STATE_PAYLOAD;
}

bool validPhaseAssertion(RestartPhaseAssertion assertion)
{
    return assertion==RestartPhaseAssertion::SINGLE_LIQUID
        ||assertion==RestartPhaseAssertion::SINGLE_VAPOR
        ||assertion==RestartPhaseAssertion::SUPERCRITICAL
        ||assertion==RestartPhaseAssertion::TWO_PHASE;
}

void appendIdentitySet(
    std::string& material,
    const RestartIdentitySet& identities)
{
    appendString(material,"identity.configuration",identities.configurationSha256);
    appendString(material,"identity.build",identities.buildSha256);
    appendString(material,"identity.mesh",identities.meshSha256);
    appendString(material,"identity.thermo",identities.thermoModelId);
    appendString(material,"identity.kij",identities.kijManifestSha256);
    appendString(material,"identity.mechanism",identities.mechanismSha256);
    appendString(material,"identity.chemistry-integrator",
        identities.chemistryIntegratorSha256);
    appendString(material,"identity.transport",identities.transportSha256);
    appendString(material,"identity.topology",identities.topologyEnvelopeSha256);
}

void appendUnsigned64(
    std::string& material,
    const std::string& key,
    std::uint64_t value)
{
    material += key+'='+std::to_string(value)+';';
}

} // namespace

const char* restartOwnerName(RestartOwner owner)
{
    switch(owner)
    {
        case RestartOwner::FLUID: return "fluid";
        case RestartOwner::TURBULENCE: return "turbulence";
        case RestartOwner::CHT: return "cht";
        case RestartOwner::CONDENSED_CARBON: return "condensed-carbon";
        case RestartOwner::CHEMISTRY: return "chemistry";
        case RestartOwner::TOPOLOGY_ENVELOPE: return "topology-envelope";
        case RestartOwner::TRANSPORT: return "transport";
        case RestartOwner::COUNT: return "invalid-owner";
    }
    return "invalid-owner";
}

const std::array<RestartOwner,restartOwnerCount>& restartOwnerRoster()
{
    static const std::array<RestartOwner,restartOwnerCount> roster = {
        RestartOwner::FLUID,
        RestartOwner::TURBULENCE,
        RestartOwner::CHT,
        RestartOwner::CONDENSED_CARBON,
        RestartOwner::CHEMISTRY,
        RestartOwner::TOPOLOGY_ENVELOPE,
        RestartOwner::TRANSPORT
    };
    return roster;
}

const char* restartOwnerDispositionName(RestartOwnerDisposition disposition)
{
    switch(disposition)
    {
        case RestartOwnerDisposition::DISABLED: return "disabled";
        case RestartOwnerDisposition::IDENTITY_ONLY: return "identity-only";
        case RestartOwnerDisposition::STATE_PAYLOAD: return "state-payload";
    }
    return "invalid-disposition";
}

RestartByteOrder nativeRestartByteOrder()
{
    const std::uint16_t value=1U;
    const unsigned char* bytes=
        reinterpret_cast<const unsigned char*>(&value);
    return bytes[0]==1U
        ?RestartByteOrder::LITTLE_ENDIAN_IEEE754
        :RestartByteOrder::BIG_ENDIAN_IEEE754;
}

const char* restartPhaseAssertionName(RestartPhaseAssertion assertion)
{
    switch(assertion)
    {
        case RestartPhaseAssertion::UNKNOWN: return "unknown";
        case RestartPhaseAssertion::SINGLE_LIQUID: return "single-liquid";
        case RestartPhaseAssertion::SINGLE_VAPOR: return "single-vapor";
        case RestartPhaseAssertion::SUPERCRITICAL: return "supercritical";
        case RestartPhaseAssertion::TWO_PHASE: return "two-phase";
        case RestartPhaseAssertion::THREE_PHASE_REQUIRED:
            return "three-phase-required";
    }
    return "unknown";
}

const char* restartManifestStatusName(RestartManifestStatus status)
{
    switch(status)
    {
        case RestartManifestStatus::SUCCESS: return "SUCCESS";
        case RestartManifestStatus::INVALID_EXPECTATION:
            return "INVALID_EXPECTATION";
        case RestartManifestStatus::INVALID_SCHEMA: return "INVALID_SCHEMA";
        case RestartManifestStatus::UNSUPPORTED_VERSION:
            return "UNSUPPORTED_VERSION";
        case RestartManifestStatus::INVALID_GENERATION:
            return "INVALID_GENERATION";
        case RestartManifestStatus::INVALID_GENERATION_STATE:
            return "INVALID_GENERATION_STATE";
        case RestartManifestStatus::UNSUPPORTED_DISTRIBUTED_RESTART:
            return "UNSUPPORTED_DISTRIBUTED_RESTART";
        case RestartManifestStatus::UNSUPPORTED_ENCODING:
            return "UNSUPPORTED_ENCODING";
        case RestartManifestStatus::TRAILING_DATA: return "TRAILING_DATA";
        case RestartManifestStatus::IDENTITY_MISMATCH:
            return "IDENTITY_MISMATCH";
        case RestartManifestStatus::INVALID_OWNER_ROSTER:
            return "INVALID_OWNER_ROSTER";
        case RestartManifestStatus::OWNER_DISPOSITION_MISMATCH:
            return "OWNER_DISPOSITION_MISMATCH";
        case RestartManifestStatus::OWNER_GENERATION_MISMATCH:
            return "OWNER_GENERATION_MISMATCH";
        case RestartManifestStatus::OWNER_IDENTITY_MISMATCH:
            return "OWNER_IDENTITY_MISMATCH";
        case RestartManifestStatus::PAYLOAD_MISSING: return "PAYLOAD_MISSING";
        case RestartManifestStatus::PAYLOAD_UNEXPECTED:
            return "PAYLOAD_UNEXPECTED";
        case RestartManifestStatus::PAYLOAD_GENERATION_MISMATCH:
            return "PAYLOAD_GENERATION_MISMATCH";
        case RestartManifestStatus::PAYLOAD_SIZE_MISMATCH:
            return "PAYLOAD_SIZE_MISMATCH";
        case RestartManifestStatus::PAYLOAD_CHECKSUM_MISMATCH:
            return "PAYLOAD_CHECKSUM_MISMATCH";
        case RestartManifestStatus::INVALID_RECOVERY_SEED:
            return "INVALID_RECOVERY_SEED";
        case RestartManifestStatus::UNSUPPORTED_PHASE_TOPOLOGY:
            return "UNSUPPORTED_PHASE_TOPOLOGY";
        case RestartManifestStatus::GENERATION_CHECKSUM_MISMATCH:
            return "GENERATION_CHECKSUM_MISMATCH";
    }
    return "INVALID_SCHEMA";
}

bool validRestartSha256(const std::string& value)
{
    return value.size()==64U&&std::all_of(
        value.begin(),value.end(),[](unsigned char c)
        { return (c>='0'&&c<='9')||(c>='a'&&c<='f'); });
}

std::string restartManifestChecksum(const RestartManifest& manifest)
{
    std::string material;
    appendString(material,"schema",manifest.schema);
    appendSize(material,"version",manifest.version);
    appendString(material,"generation",manifest.generationId);
    appendDouble(material,"generation.physical-time",
        manifest.generationState.physicalTime);
    appendUnsigned64(material,"generation.accepted-iteration",
        manifest.generationState.acceptedIteration);
    appendSize(material,"generation.previous-dt-present",
        manifest.generationState.previousAcceptedDtPresent?1U:0U);
    appendDouble(material,"generation.previous-dt",
        manifest.generationState.previousAcceptedDt);
    appendSize(material,"execution-scope",
        static_cast<std::size_t>(manifest.executionScope));
    appendSize(material,"rank-count",manifest.rankCount);
    appendSize(material,"rank-index",manifest.rankIndex);
    appendSize(material,"collective-commit",manifest.collectiveCommit?1U:0U);
    appendString(material,"decomposition",manifest.decompositionId);
    appendSize(material,"byte-order",
        static_cast<std::size_t>(manifest.byteOrder));
    appendString(material,"floating-point",manifest.floatingPointFormat);
    appendSize(material,"trailing-permitted",
        manifest.trailingDataPermitted?1U:0U);
    appendUnsigned64(material,"cell-count",manifest.cellCount);
    appendIdentitySet(material,manifest.identities);
    appendSize(material,"owner-count",manifest.owners.size());
    for(std::size_t i=0;i<manifest.owners.size();++i)
    {
        const RestartOwnerRecord& owner=manifest.owners[i];
        const std::string prefix="owner["+std::to_string(i)+"].";
        appendSize(material,prefix+"kind",static_cast<std::size_t>(owner.owner));
        appendSize(material,prefix+"disposition",
            static_cast<std::size_t>(owner.disposition));
        appendSize(material,prefix+"publication-order",owner.publicationOrder);
        appendString(material,prefix+"generation",owner.generationId);
        appendString(material,prefix+"identity",owner.identitySha256);
        appendUnsigned64(material,prefix+"payload-bytes",owner.payloadBytes);
        appendString(material,prefix+"payload-sha256",owner.payloadSha256);
    }
    appendSize(material,"seed-count",manifest.recoverySeeds.size());
    for(std::size_t i=0;i<manifest.recoverySeeds.size();++i)
    {
        const RestartRecoverySeed& seed=manifest.recoverySeeds[i];
        const std::string prefix="seed["+std::to_string(i)+"].";
        appendUnsigned64(material,prefix+"cell",seed.cell);
        appendDouble(material,prefix+"pressure",seed.pressure);
        appendDouble(material,prefix+"temperature",seed.temperature);
        appendSize(material,prefix+"phase",
            static_cast<std::size_t>(seed.expectedPhase));
        appendString(material,prefix+"phase-identity",
            seed.phaseIdentitySha256);
    }
    return thermoProvenance::sha256Hex(material);
}

RestartManifestValidationResult RestartManifestValidator::validate(
    const RestartCandidatePackage& candidate,
    const RestartManifestExpectation& expected)
{
    const RestartManifest& manifest=candidate.manifest;
    const auto& roster=restartOwnerRoster();

    if(!validIdentities(expected.identities)
       ||expected.byteOrder!=nativeRestartByteOrder()
       ||expected.floatingPointFormat!=RestartManifest::floatingPointFormatId
       ||expected.decompositionId
            !=RestartManifest::singleProcessDecompositionId
       ||expected.cellCount==0U
       ||expected.maximumAcceptedIteration
            >static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
        return failure(RestartManifestStatus::INVALID_EXPECTATION,
            "invalid restart expectation");
    for(std::size_t i=0;i<roster.size();++i)
    {
        const RestartOwnerExpectation& owner=expected.owners[i];
        if(owner.owner!=roster[i]||!validDisposition(owner.disposition)
           ||!validRestartSha256(owner.identitySha256))
            return failure(RestartManifestStatus::INVALID_EXPECTATION,
                "invalid expected owner roster");
    }
    if(expected.owners[static_cast<std::size_t>(RestartOwner::FLUID)]
            .disposition!=RestartOwnerDisposition::STATE_PAYLOAD)
        return failure(RestartManifestStatus::INVALID_EXPECTATION,
            "fluid state payload is required for restart");

    if(manifest.schema!=RestartManifest::schemaId)
        return failure(RestartManifestStatus::INVALID_SCHEMA,
            "restart schema mismatch");
    if(manifest.version!=RestartManifest::currentVersion)
        return failure(RestartManifestStatus::UNSUPPORTED_VERSION,
            "restart version is unsupported");
    if(!validToken(manifest.generationId,128U)
       ||!validRestartSha256(manifest.generationChecksumSha256))
        return failure(RestartManifestStatus::INVALID_GENERATION,
            "invalid restart generation identity");
    const RestartGenerationState& generation=manifest.generationState;
    if(!std::isfinite(generation.physicalTime)
       ||generation.physicalTime<0.0
       ||std::signbit(generation.physicalTime)
       ||generation.acceptedIteration>expected.maximumAcceptedIteration
       ||generation.previousAcceptedDtPresent
            !=expected.requirePreviousAcceptedDt
       ||(generation.previousAcceptedDtPresent
            &&(!std::isfinite(generation.previousAcceptedDt)
               ||generation.previousAcceptedDt<=0.0))
       ||(!generation.previousAcceptedDtPresent
            &&(generation.previousAcceptedDt!=0.0
               ||std::signbit(generation.previousAcceptedDt))))
        return failure(RestartManifestStatus::INVALID_GENERATION_STATE,
            "invalid or unsupported restart clock/step state");
    if(manifest.executionScope!=RestartExecutionScope::SINGLE_PROCESS
       ||manifest.rankCount!=1U||manifest.rankIndex!=0U
       ||manifest.collectiveCommit
       ||manifest.decompositionId
            !=RestartManifest::singleProcessDecompositionId
       ||manifest.decompositionId!=expected.decompositionId)
        return failure(
            RestartManifestStatus::UNSUPPORTED_DISTRIBUTED_RESTART,
            "MPI, decomposition, and collective commit are unsupported");
    if(manifest.byteOrder!=expected.byteOrder
       ||manifest.byteOrder!=nativeRestartByteOrder()
       ||manifest.floatingPointFormat!=expected.floatingPointFormat
       ||!std::numeric_limits<double>::is_iec559||sizeof(double)!=8U)
        return failure(RestartManifestStatus::UNSUPPORTED_ENCODING,
            "restart numeric encoding is unsupported");
    if(manifest.trailingDataPermitted||candidate.manifestTrailingBytes!=0U)
        return failure(RestartManifestStatus::TRAILING_DATA,
            "manifest trailing data is forbidden");
    if(manifest.cellCount!=expected.cellCount
       ||!validIdentities(manifest.identities)
       ||!sameIdentities(manifest.identities,expected.identities))
        return failure(RestartManifestStatus::IDENTITY_MISMATCH,
            "restart configuration identity mismatch");

    if(manifest.owners.size()!=roster.size())
        return failure(RestartManifestStatus::INVALID_OWNER_ROSTER,
            "restart owner count mismatch");
    for(std::size_t i=0;i<roster.size();++i)
    {
        const RestartOwnerRecord& owner=manifest.owners[i];
        const RestartOwnerExpectation& ownerExpected=expected.owners[i];
        if(owner.owner!=roster[i]
           ||owner.publicationOrder!=static_cast<std::uint32_t>(i))
            return failure(RestartManifestStatus::INVALID_OWNER_ROSTER,
                "restart owner order mismatch");
        if(!validDisposition(owner.disposition)
           ||owner.disposition!=ownerExpected.disposition)
            return failure(
                RestartManifestStatus::OWNER_DISPOSITION_MISMATCH,
                std::string("owner disposition mismatch: ")
                    +restartOwnerName(owner.owner));
        if(owner.generationId!=manifest.generationId)
            return failure(RestartManifestStatus::OWNER_GENERATION_MISMATCH,
                std::string("owner generation mismatch: ")
                    +restartOwnerName(owner.owner));
        if(!validRestartSha256(owner.identitySha256)
           ||owner.identitySha256!=ownerExpected.identitySha256)
            return failure(RestartManifestStatus::OWNER_IDENTITY_MISMATCH,
                std::string("owner identity mismatch: ")
                    +restartOwnerName(owner.owner));
        if(!validRestartSha256(owner.payloadSha256))
            return failure(RestartManifestStatus::PAYLOAD_CHECKSUM_MISMATCH,
                std::string("invalid payload checksum: ")
                    +restartOwnerName(owner.owner));
        if(owner.disposition==RestartOwnerDisposition::STATE_PAYLOAD)
        {
            if(owner.payloadBytes==0U)
                return failure(RestartManifestStatus::PAYLOAD_MISSING,
                    std::string("empty state payload: ")
                        +restartOwnerName(owner.owner));
        }
        else if(owner.payloadBytes!=0U
                ||owner.payloadSha256
                    !=thermoProvenance::sha256Hex(std::string()))
            return failure(RestartManifestStatus::PAYLOAD_UNEXPECTED,
                std::string("non-state owner has payload: ")
                    +restartOwnerName(owner.owner));
    }

    std::array<const RestartPayloadCandidate*,restartOwnerCount> payloads{};
    for(const RestartPayloadCandidate& payload:candidate.payloads)
    {
        const std::size_t ownerIndex=static_cast<std::size_t>(payload.owner);
        if(ownerIndex>=restartOwnerCount||payloads[ownerIndex]!=nullptr)
            return failure(RestartManifestStatus::PAYLOAD_UNEXPECTED,
                "unknown or duplicate owner payload");
        payloads[ownerIndex]=&payload;
    }
    for(std::size_t i=0;i<roster.size();++i)
    {
        const RestartOwnerRecord& owner=manifest.owners[i];
        const RestartPayloadCandidate* payload=payloads[i];
        if(owner.disposition!=RestartOwnerDisposition::STATE_PAYLOAD)
        {
            if(payload!=nullptr)
                return failure(RestartManifestStatus::PAYLOAD_UNEXPECTED,
                    std::string("unexpected payload: ")
                        +restartOwnerName(owner.owner));
            continue;
        }
        if(payload==nullptr)
            return failure(RestartManifestStatus::PAYLOAD_MISSING,
                std::string("missing payload: ")
                    +restartOwnerName(owner.owner));
        if(payload->generationId!=manifest.generationId)
            return failure(RestartManifestStatus::PAYLOAD_GENERATION_MISMATCH,
                std::string("payload generation mismatch: ")
                    +restartOwnerName(owner.owner));
        if(payload->consumedBytes!=payload->bytes.size())
            return failure(RestartManifestStatus::TRAILING_DATA,
                std::string("payload trailing data: ")
                    +restartOwnerName(owner.owner));
        if(owner.payloadBytes!=payload->bytes.size())
            return failure(RestartManifestStatus::PAYLOAD_SIZE_MISMATCH,
                std::string("payload size mismatch: ")
                    +restartOwnerName(owner.owner));
        if(owner.payloadSha256
            !=thermoProvenance::sha256Hex(payload->bytes))
            return failure(RestartManifestStatus::PAYLOAD_CHECKSUM_MISMATCH,
                std::string("payload checksum mismatch: ")
                    +restartOwnerName(owner.owner));
    }

    if(manifest.recoverySeeds.size()!=manifest.cellCount)
        return failure(RestartManifestStatus::INVALID_RECOVERY_SEED,
            "recovery seed count mismatch");
    for(std::size_t i=0;i<manifest.recoverySeeds.size();++i)
    {
        const RestartRecoverySeed& seed=manifest.recoverySeeds[i];
        if(seed.cell!=i||!std::isfinite(seed.pressure)||seed.pressure<=0.0
           ||!std::isfinite(seed.temperature)||seed.temperature<=0.0
           ||!validRestartSha256(seed.phaseIdentitySha256))
            return failure(RestartManifestStatus::INVALID_RECOVERY_SEED,
                "invalid or unordered rho-e recovery seed");
        if(seed.expectedPhase==RestartPhaseAssertion::THREE_PHASE_REQUIRED)
            return failure(
                RestartManifestStatus::UNSUPPORTED_PHASE_TOPOLOGY,
                "three-phase-required state is outside the restart envelope");
        if(!validPhaseAssertion(seed.expectedPhase))
            return failure(RestartManifestStatus::INVALID_RECOVERY_SEED,
                "unknown recovery phase assertion");
    }

    if(manifest.generationChecksumSha256
       !=restartManifestChecksum(manifest))
        return failure(
            RestartManifestStatus::GENERATION_CHECKSUM_MISMATCH,
            "restart generation checksum mismatch");

    RestartManifestValidationResult result;
    result.status=RestartManifestStatus::SUCCESS;
    result.detail="candidate-only manifest validated; no state published";
    result.publicationOrder.assign(roster.begin(),roster.end());
    result.generationState=manifest.generationState;
    return result;
}

} // End namespace Foam
