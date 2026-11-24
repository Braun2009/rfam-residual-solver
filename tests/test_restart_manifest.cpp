#include "RestartManifest.H"
#include "ThermoDataProvenance.H"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>

using namespace Foam;

namespace
{

int checks=0;
int failures=0;

void check(bool condition,const std::string& label)
{
    ++checks;
    std::cout << "  [" << (condition?"PASS":"FAIL") << "] "
              << label << '\n';
    if(!condition) ++failures;
}

std::string hash(const std::string& value)
{
    return thermoProvenance::sha256Hex(value);
}

RestartCandidatePackage validCandidate()
{
    RestartCandidatePackage candidate;
    RestartManifest& manifest=candidate.manifest;
    manifest.generationId="generation-20260723-0001";
    manifest.generationState.physicalTime=1.25e-5;
    manifest.generationState.acceptedIteration=42;
    manifest.generationState.previousAcceptedDtPresent=true;
    manifest.generationState.previousAcceptedDt=1.0e-8;
    manifest.byteOrder=nativeRestartByteOrder();
    manifest.cellCount=2;
    manifest.identities.configurationSha256=hash("configuration");
    manifest.identities.buildSha256=hash("build");
    manifest.identities.meshSha256=hash("mesh");
    manifest.identities.thermoModelId=
        "peng-robinson-rho-e:sha256="+hash("thermo");
    manifest.identities.kijManifestSha256=hash("kij");
    manifest.identities.mechanismSha256=hash("mechanism");
    manifest.identities.chemistryIntegratorSha256=hash("chemistry-integrator");
    manifest.identities.transportSha256=hash("transport");
    manifest.identities.topologyEnvelopeSha256=hash("topology-envelope");

    const auto& roster=restartOwnerRoster();
    manifest.owners.resize(roster.size());
    for(std::size_t i=0;i<roster.size();++i)
    {
        RestartOwnerRecord& owner=manifest.owners[i];
        owner.owner=roster[i];
        owner.publicationOrder=static_cast<std::uint32_t>(i);
        owner.generationId=manifest.generationId;
        owner.identitySha256=hash(
            std::string("owner-")+restartOwnerName(owner.owner));
        owner.payloadSha256=hash(std::string());
    }
    manifest.owners[0].disposition=RestartOwnerDisposition::STATE_PAYLOAD;
    manifest.owners[1].disposition=RestartOwnerDisposition::STATE_PAYLOAD;
    manifest.owners[2].disposition=RestartOwnerDisposition::DISABLED;
    manifest.owners[3].disposition=RestartOwnerDisposition::DISABLED;
    manifest.owners[4].disposition=RestartOwnerDisposition::IDENTITY_ONLY;
    manifest.owners[5].disposition=RestartOwnerDisposition::IDENTITY_ONLY;
    manifest.owners[6].disposition=RestartOwnerDisposition::IDENTITY_ONLY;

    candidate.payloads={
        {RestartOwner::FLUID,manifest.generationId,"fluid-state-v1",0},
        {RestartOwner::TURBULENCE,manifest.generationId,
            "turbulence-state-v1",0}
    };
    for(RestartPayloadCandidate& payload:candidate.payloads)
    {
        payload.consumedBytes=payload.bytes.size();
        RestartOwnerRecord& owner=
            manifest.owners[static_cast<std::size_t>(payload.owner)];
        owner.payloadBytes=payload.bytes.size();
        owner.payloadSha256=hash(payload.bytes);
    }

    manifest.recoverySeeds={
        {0,5.0e6,155.0,RestartPhaseAssertion::SINGLE_LIQUID,
            hash("phase-cell-0")},
        {1,5.0e6,180.0,RestartPhaseAssertion::TWO_PHASE,
            hash("phase-cell-1-quality-x-y")}
    };
    manifest.generationChecksumSha256=restartManifestChecksum(manifest);
    return candidate;
}

RestartManifestExpectation expectation(
    const RestartCandidatePackage& candidate)
{
    RestartManifestExpectation expected;
    expected.identities=candidate.manifest.identities;
    expected.byteOrder=nativeRestartByteOrder();
    expected.cellCount=candidate.manifest.cellCount;
    expected.requirePreviousAcceptedDt=true;
    const auto& roster=restartOwnerRoster();
    for(std::size_t i=0;i<roster.size();++i)
    {
        expected.owners[i].owner=roster[i];
        expected.owners[i].disposition=
            candidate.manifest.owners[i].disposition;
        expected.owners[i].identitySha256=
            candidate.manifest.owners[i].identitySha256;
    }
    return expected;
}

void refreshChecksum(RestartCandidatePackage& candidate)
{
    candidate.manifest.generationChecksumSha256=
        restartManifestChecksum(candidate.manifest);
}

void expectStatus(
    RestartCandidatePackage candidate,
    const RestartManifestExpectation& expected,
    RestartManifestStatus status,
    const std::string& label,
    bool refresh=true)
{
    if(refresh) refreshChecksum(candidate);
    const RestartManifestValidationResult result=
        RestartManifestValidator::validate(candidate,expected);
    check(result.status==status,
        label+" expected="+restartManifestStatusName(status)
            +" actual="+restartManifestStatusName(result.status));
    check(result.publicationOrder.empty(),
        label+" exposes no publication order");
    check(!result.generationState.has_value(),
        label+" exposes no generation state");
}

} // namespace

int main()
{
    const RestartCandidatePackage baseline=validCandidate();
    const RestartManifestExpectation expected=expectation(baseline);
    const RestartManifestValidationResult accepted=
        RestartManifestValidator::validate(baseline,expected);
    check(accepted.success(),"complete candidate-only manifest is accepted");
    check(accepted.publicationOrder==std::vector<RestartOwner>(
            restartOwnerRoster().begin(),restartOwnerRoster().end()),
        "success exposes only the canonical seven-owner publication order");
    check(accepted.detail.find("no state published")!=std::string::npos,
        "success result remains explicitly candidate-only");
    check(accepted.generationState.has_value()
            &&accepted.generationState->physicalTime
                ==baseline.manifest.generationState.physicalTime
            &&accepted.generationState->acceptedIteration==42
            &&accepted.generationState->previousAcceptedDtPresent
            &&accepted.generationState->previousAcceptedDt==1.0e-8,
        "success exposes validated clock state for last publication");

    RestartCandidatePackage candidate=baseline;
    candidate.manifest.schema="rfam.restart-manifest/v0";
    expectStatus(candidate,expected,RestartManifestStatus::INVALID_SCHEMA,
        "wrong schema");

    candidate=baseline;
    candidate.manifest.schema="rfam.restart-manifest/v1";
    candidate.manifest.version=1;
    expectStatus(candidate,expected,RestartManifestStatus::INVALID_SCHEMA,
        "superseded v1 schema");

    candidate=baseline; candidate.manifest.version=3;
    expectStatus(candidate,expected,RestartManifestStatus::UNSUPPORTED_VERSION,
        "future version");

    candidate=baseline; candidate.manifest.generationId="bad generation";
    for(RestartOwnerRecord& owner:candidate.manifest.owners)
        owner.generationId=candidate.manifest.generationId;
    for(RestartPayloadCandidate& payload:candidate.payloads)
        payload.generationId=candidate.manifest.generationId;
    expectStatus(candidate,expected,RestartManifestStatus::INVALID_GENERATION,
        "unsafe generation token");

    candidate=baseline;
    candidate.manifest.generationState.physicalTime=
        std::numeric_limits<double>::quiet_NaN();
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_GENERATION_STATE,
        "nonfinite physical time");
    candidate=baseline;
    candidate.manifest.generationState.physicalTime=-1.0;
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_GENERATION_STATE,
        "negative physical time");
    candidate=baseline;
    candidate.manifest.generationState.physicalTime=-0.0;
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_GENERATION_STATE,
        "negative-zero physical time");
    candidate=baseline;
    candidate.manifest.generationState.physicalTime=
        std::numeric_limits<double>::infinity();
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_GENERATION_STATE,
        "infinite physical time");
    candidate=baseline;
    candidate.manifest.generationState.acceptedIteration=
        static_cast<std::uint64_t>(std::numeric_limits<int>::max())+1U;
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_GENERATION_STATE,
        "iteration exceeds current solver range");
    candidate=baseline;
    candidate.manifest.generationState.previousAcceptedDtPresent=false;
    candidate.manifest.generationState.previousAcceptedDt=0.0;
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_GENERATION_STATE,
        "required previous accepted dt is missing");
    candidate=baseline;
    candidate.manifest.generationState.previousAcceptedDt=0.0;
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_GENERATION_STATE,
        "previous accepted dt is nonpositive");
    candidate=baseline;
    candidate.manifest.generationState.previousAcceptedDt=
        std::numeric_limits<double>::infinity();
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_GENERATION_STATE,
        "previous accepted dt is nonfinite");

    RestartCandidatePackage noPreviousDt=baseline;
    noPreviousDt.manifest.generationState.previousAcceptedDtPresent=false;
    noPreviousDt.manifest.generationState.previousAcceptedDt=0.0;
    refreshChecksum(noPreviousDt);
    RestartManifestExpectation noPreviousDtExpected=expected;
    noPreviousDtExpected.requirePreviousAcceptedDt=false;
    check(RestartManifestValidator::validate(
            noPreviousDt,noPreviousDtExpected).success(),
        "configuration without a dt-consuming boundary accepts canonical absence");
    candidate=baseline;
    expectStatus(candidate,noPreviousDtExpected,
        RestartManifestStatus::INVALID_GENERATION_STATE,
        "configuration rejects unexpected previous accepted dt");
    candidate=noPreviousDt;
    candidate.manifest.generationState.previousAcceptedDt=-0.0;
    expectStatus(candidate,noPreviousDtExpected,
        RestartManifestStatus::INVALID_GENERATION_STATE,
        "absent previous dt requires canonical positive zero");

    candidate=baseline;
    candidate.manifest.executionScope=RestartExecutionScope::MPI_PARTITIONED;
    expectStatus(candidate,expected,
        RestartManifestStatus::UNSUPPORTED_DISTRIBUTED_RESTART,
        "MPI execution scope");
    candidate=baseline; candidate.manifest.rankCount=2;
    expectStatus(candidate,expected,
        RestartManifestStatus::UNSUPPORTED_DISTRIBUTED_RESTART,
        "multiple ranks");
    candidate=baseline; candidate.manifest.rankIndex=1;
    expectStatus(candidate,expected,
        RestartManifestStatus::UNSUPPORTED_DISTRIBUTED_RESTART,
        "nonzero rank");
    candidate=baseline; candidate.manifest.collectiveCommit=true;
    expectStatus(candidate,expected,
        RestartManifestStatus::UNSUPPORTED_DISTRIBUTED_RESTART,
        "collective commit request");
    candidate=baseline; candidate.manifest.decompositionId="mpi-4-way";
    expectStatus(candidate,expected,
        RestartManifestStatus::UNSUPPORTED_DISTRIBUTED_RESTART,
        "partitioned decomposition");

    candidate=baseline;
    candidate.manifest.byteOrder=
        nativeRestartByteOrder()==RestartByteOrder::LITTLE_ENDIAN_IEEE754
        ?RestartByteOrder::BIG_ENDIAN_IEEE754
        :RestartByteOrder::LITTLE_ENDIAN_IEEE754;
    expectStatus(candidate,expected,RestartManifestStatus::UNSUPPORTED_ENCODING,
        "opposite byte order");
    candidate=baseline; candidate.manifest.floatingPointFormat="native-double";
    expectStatus(candidate,expected,RestartManifestStatus::UNSUPPORTED_ENCODING,
        "unbound floating-point format");

    candidate=baseline; candidate.manifestTrailingBytes=1;
    expectStatus(candidate,expected,RestartManifestStatus::TRAILING_DATA,
        "manifest trailing byte");
    candidate=baseline; candidate.manifest.trailingDataPermitted=true;
    expectStatus(candidate,expected,RestartManifestStatus::TRAILING_DATA,
        "manifest attempts to permit trailing data");

    candidate=baseline;
    candidate.manifest.identities.configurationSha256=hash("other-config");
    expectStatus(candidate,expected,RestartManifestStatus::IDENTITY_MISMATCH,
        "configuration identity mismatch");
    candidate=baseline;
    candidate.manifest.identities.buildSha256=hash("other-build");
    expectStatus(candidate,expected,RestartManifestStatus::IDENTITY_MISMATCH,
        "build identity mismatch");
    candidate=baseline;
    candidate.manifest.identities.meshSha256=hash("other-mesh");
    expectStatus(candidate,expected,RestartManifestStatus::IDENTITY_MISMATCH,
        "mesh identity mismatch");
    candidate=baseline;
    candidate.manifest.identities.thermoModelId="other-thermo";
    expectStatus(candidate,expected,RestartManifestStatus::IDENTITY_MISMATCH,
        "thermo identity mismatch");
    candidate=baseline;
    candidate.manifest.identities.kijManifestSha256=hash("other-kij");
    expectStatus(candidate,expected,RestartManifestStatus::IDENTITY_MISMATCH,
        "kij identity mismatch");
    candidate=baseline;
    candidate.manifest.identities.mechanismSha256=hash("other-mechanism");
    expectStatus(candidate,expected,RestartManifestStatus::IDENTITY_MISMATCH,
        "mechanism identity mismatch");
    candidate=baseline;
    candidate.manifest.identities.chemistryIntegratorSha256=
        hash("other-integrator");
    expectStatus(candidate,expected,RestartManifestStatus::IDENTITY_MISMATCH,
        "chemistry integrator identity mismatch");
    candidate=baseline;
    candidate.manifest.identities.transportSha256=hash("other-transport");
    expectStatus(candidate,expected,RestartManifestStatus::IDENTITY_MISMATCH,
        "transport identity mismatch");
    candidate=baseline;
    candidate.manifest.identities.topologyEnvelopeSha256=
        hash("other-topology");
    expectStatus(candidate,expected,RestartManifestStatus::IDENTITY_MISMATCH,
        "topology envelope identity mismatch");

    candidate=baseline; candidate.manifest.owners.pop_back();
    expectStatus(candidate,expected,RestartManifestStatus::INVALID_OWNER_ROSTER,
        "missing transport owner");
    candidate=baseline; candidate.manifest.owners[2].owner=RestartOwner::TURBULENCE;
    expectStatus(candidate,expected,RestartManifestStatus::INVALID_OWNER_ROSTER,
        "duplicate owner in roster");
    candidate=baseline; candidate.manifest.owners[2].publicationOrder=6;
    expectStatus(candidate,expected,RestartManifestStatus::INVALID_OWNER_ROSTER,
        "publication order drift");
    candidate=baseline;
    candidate.manifest.owners[3].disposition=
        RestartOwnerDisposition::IDENTITY_ONLY;
    expectStatus(candidate,expected,
        RestartManifestStatus::OWNER_DISPOSITION_MISMATCH,
        "owner disposition rebind");
    candidate=baseline;
    candidate.manifest.owners[4].generationId="other-generation";
    expectStatus(candidate,expected,
        RestartManifestStatus::OWNER_GENERATION_MISMATCH,
        "owner generation mismatch");
    candidate=baseline;
    candidate.manifest.owners[5].identitySha256=hash("other-owner-identity");
    expectStatus(candidate,expected,
        RestartManifestStatus::OWNER_IDENTITY_MISMATCH,
        "owner identity mismatch");
    candidate=baseline;
    candidate.manifest.owners[5].identitySha256="not-a-sha";
    expectStatus(candidate,expected,
        RestartManifestStatus::OWNER_IDENTITY_MISMATCH,
        "invalid owner identity");

    candidate=baseline; candidate.payloads.erase(candidate.payloads.begin());
    expectStatus(candidate,expected,RestartManifestStatus::PAYLOAD_MISSING,
        "missing fluid payload");
    candidate=baseline; candidate.payloads.push_back(candidate.payloads.front());
    expectStatus(candidate,expected,RestartManifestStatus::PAYLOAD_UNEXPECTED,
        "duplicate fluid payload");
    candidate=baseline;
    candidate.payloads.push_back({RestartOwner::CHT,
        candidate.manifest.generationId,"",0});
    expectStatus(candidate,expected,RestartManifestStatus::PAYLOAD_UNEXPECTED,
        "disabled CHT supplies payload");
    candidate=baseline;
    candidate.payloads.push_back({RestartOwner::CHEMISTRY,
        candidate.manifest.generationId,"identity-only-data",18});
    expectStatus(candidate,expected,RestartManifestStatus::PAYLOAD_UNEXPECTED,
        "identity-only chemistry supplies payload");
    candidate=baseline;
    candidate.payloads.push_back({static_cast<RestartOwner>(99),
        candidate.manifest.generationId,"unknown-owner",13});
    expectStatus(candidate,expected,RestartManifestStatus::PAYLOAD_UNEXPECTED,
        "unknown owner supplies payload");
    candidate=baseline; candidate.payloads[0].generationId="other-generation";
    expectStatus(candidate,expected,
        RestartManifestStatus::PAYLOAD_GENERATION_MISMATCH,
        "payload generation mismatch");
    candidate=baseline; candidate.payloads[0].consumedBytes-=1;
    expectStatus(candidate,expected,RestartManifestStatus::TRAILING_DATA,
        "fluid payload trailing byte");
    candidate=baseline; candidate.manifest.owners[0].payloadBytes+=1;
    expectStatus(candidate,expected,
        RestartManifestStatus::PAYLOAD_SIZE_MISMATCH,
        "fluid payload size mismatch");
    candidate=baseline;
    candidate.manifest.owners[0].payloadBytes=0;
    candidate.manifest.owners[0].payloadSha256=hash(std::string());
    expectStatus(candidate,expected,RestartManifestStatus::PAYLOAD_MISSING,
        "empty fluid state payload");
    candidate=baseline; candidate.payloads[0].bytes[0]='F';
    expectStatus(candidate,expected,
        RestartManifestStatus::PAYLOAD_CHECKSUM_MISMATCH,
        "fluid payload checksum mismatch");

    candidate=baseline; candidate.manifest.recoverySeeds.pop_back();
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_RECOVERY_SEED,
        "missing rho-e recovery seed");
    candidate=baseline; candidate.manifest.cellCount+=1;
    expectStatus(candidate,expected,RestartManifestStatus::IDENTITY_MISMATCH,
        "cell count identity mismatch");
    candidate=baseline; candidate.manifest.recoverySeeds[1].cell=0;
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_RECOVERY_SEED,
        "duplicate recovery seed cell");
    candidate=baseline;
    candidate.manifest.recoverySeeds[0].pressure=
        std::numeric_limits<double>::quiet_NaN();
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_RECOVERY_SEED,
        "nonfinite pressure seed");
    candidate=baseline; candidate.manifest.recoverySeeds[0].temperature=0.0;
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_RECOVERY_SEED,
        "nonpositive temperature seed");
    candidate=baseline;
    candidate.manifest.recoverySeeds[0].temperature=
        std::numeric_limits<double>::infinity();
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_RECOVERY_SEED,
        "nonfinite temperature seed");
    candidate=baseline;
    candidate.manifest.recoverySeeds[0].expectedPhase=
        RestartPhaseAssertion::UNKNOWN;
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_RECOVERY_SEED,
        "unknown phase assertion");
    candidate=baseline;
    candidate.manifest.recoverySeeds[1].expectedPhase=
        RestartPhaseAssertion::THREE_PHASE_REQUIRED;
    expectStatus(candidate,expected,
        RestartManifestStatus::UNSUPPORTED_PHASE_TOPOLOGY,
        "three-phase-required state");
    candidate=baseline;
    candidate.manifest.recoverySeeds[1].phaseIdentitySha256="not-a-sha";
    expectStatus(candidate,expected,
        RestartManifestStatus::INVALID_RECOVERY_SEED,
        "unbound phase identity");

    candidate=baseline;
    candidate.manifest.generationChecksumSha256=hash("wrong-generation");
    expectStatus(candidate,expected,
        RestartManifestStatus::GENERATION_CHECKSUM_MISMATCH,
        "generation checksum mismatch",false);

    RestartManifestExpectation invalidExpected=expected;
    invalidExpected.owners[0].owner=RestartOwner::TURBULENCE;
    const RestartManifestValidationResult invalidExpectation=
        RestartManifestValidator::validate(baseline,invalidExpected);
    check(invalidExpectation.status==RestartManifestStatus::INVALID_EXPECTATION,
        "invalid runtime expectation fails closed");
    check(invalidExpectation.publicationOrder.empty(),
        "invalid expectation exposes no publication order");
    check(!invalidExpectation.generationState.has_value(),
        "invalid expectation exposes no generation state");

    invalidExpected=expected;
    invalidExpected.owners[0].disposition=
        RestartOwnerDisposition::IDENTITY_ONLY;
    const RestartManifestValidationResult missingFluidExpectation=
        RestartManifestValidator::validate(baseline,invalidExpected);
    check(missingFluidExpectation.status
            ==RestartManifestStatus::INVALID_EXPECTATION,
        "runtime expectation cannot omit the fluid state payload");
    check(missingFluidExpectation.publicationOrder.empty(),
        "missing fluid expectation exposes no publication order");
    check(!missingFluidExpectation.generationState.has_value(),
        "missing fluid expectation exposes no generation state");

    invalidExpected=expected;
    invalidExpected.maximumAcceptedIteration=
        static_cast<std::uint64_t>(std::numeric_limits<int>::max())+1U;
    const RestartManifestValidationResult invalidIterationExpectation=
        RestartManifestValidator::validate(baseline,invalidExpected);
    check(invalidIterationExpectation.status
            ==RestartManifestStatus::INVALID_EXPECTATION,
        "runtime cannot advertise an iteration range it cannot publish");
    check(invalidIterationExpectation.publicationOrder.empty()
            &&!invalidIterationExpectation.generationState.has_value(),
        "invalid iteration expectation exposes no publication data");

    candidate=baseline;
    candidate.payloads[0].bytes="new-valid-fluid-state";
    candidate.payloads[0].consumedBytes=candidate.payloads[0].bytes.size();
    candidate.manifest.owners[0].payloadBytes=candidate.payloads[0].bytes.size();
    candidate.manifest.owners[0].payloadSha256=hash(candidate.payloads[0].bytes);
    refreshChecksum(candidate);
    check(RestartManifestValidator::validate(candidate,expected).success(),
        "self-consistent new state payload remains valid data, not identity drift");

    if(failures==0)
    {
        std::cout << "RESTART_MANIFEST_CANDIDATE_ONLY_PASS checks="
                  << checks << '\n';
        return 0;
    }
    std::cout << "RESTART_MANIFEST_CANDIDATE_ONLY_FAIL checks="
              << checks << " failures=" << failures << '\n';
    return 1;
}
