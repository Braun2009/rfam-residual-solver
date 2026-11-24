#include "RestartPublication.H"

#include "ThermoDataProvenance.H"

#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

namespace
{

using namespace Foam;

int checks = 0;

void check(bool condition, const std::string& label)
{
    if (!condition) throw std::runtime_error(label);
    ++checks;
}

std::string hash(const std::string& value)
{
    return thermoProvenance::sha256Hex(value);
}

RestartCandidatePackage candidate(
    const std::string& generationId,
    double physicalTime,
    std::uint64_t acceptedIteration)
{
    RestartCandidatePackage result;
    RestartManifest& manifest = result.manifest;
    manifest.generationId = generationId;
    manifest.generationState.physicalTime = physicalTime;
    manifest.generationState.acceptedIteration = acceptedIteration;
    manifest.generationState.previousAcceptedDtPresent = true;
    manifest.generationState.previousAcceptedDt = 1.0e-8;
    manifest.byteOrder = nativeRestartByteOrder();
    manifest.cellCount = 1U;
    manifest.identities.configurationSha256 = hash("configuration");
    manifest.identities.buildSha256 = hash("build");
    manifest.identities.meshSha256 = hash("mesh");
    manifest.identities.thermoModelId =
        "peng-robinson-rho-e:sha256=" + hash("thermo");
    manifest.identities.kijManifestSha256 = hash("kij");
    manifest.identities.mechanismSha256 = hash("mechanism");
    manifest.identities.chemistryIntegratorSha256 =
        hash("chemistry-integrator");
    manifest.identities.transportSha256 = hash("transport");
    manifest.identities.topologyEnvelopeSha256 = hash("topology");

    const auto& roster = restartOwnerRoster();
    manifest.owners.resize(roster.size());
    for (std::size_t index = 0; index < roster.size(); ++index)
    {
        RestartOwnerRecord& owner = manifest.owners[index];
        owner.owner = roster[index];
        owner.disposition = RestartOwnerDisposition::IDENTITY_ONLY;
        owner.publicationOrder = static_cast<std::uint32_t>(index);
        owner.generationId = generationId;
        owner.identitySha256 = hash(
            std::string("owner-") + restartOwnerName(owner.owner));
        owner.payloadSha256 = hash(std::string());
    }
    manifest.owners[static_cast<std::size_t>(RestartOwner::FLUID)]
        .disposition = RestartOwnerDisposition::STATE_PAYLOAD;
    manifest.owners[static_cast<std::size_t>(RestartOwner::TURBULENCE)]
        .disposition = RestartOwnerDisposition::STATE_PAYLOAD;
    manifest.owners[static_cast<std::size_t>(RestartOwner::CHT)]
        .disposition = RestartOwnerDisposition::DISABLED;
    manifest.owners[static_cast<std::size_t>(RestartOwner::CONDENSED_CARBON)]
        .disposition = RestartOwnerDisposition::DISABLED;

    result.payloads = {
        {RestartOwner::FLUID,generationId,
            "fluid-payload-" + generationId,0U},
        {RestartOwner::TURBULENCE,generationId,
            "turbulence-payload-" + generationId,0U}};
    for (RestartPayloadCandidate& payload : result.payloads)
    {
        payload.consumedBytes = payload.bytes.size();
        RestartOwnerRecord& owner = manifest.owners[
            static_cast<std::size_t>(payload.owner)];
        owner.payloadBytes = payload.bytes.size();
        owner.payloadSha256 = hash(payload.bytes);
    }
    manifest.recoverySeeds = {{
        0U,5.0e6,170.0,RestartPhaseAssertion::TWO_PHASE,
        hash("phase-cell-0-quality-x-y")}};
    manifest.generationChecksumSha256 = restartManifestChecksum(manifest);
    return result;
}

RestartManifestExpectation expectation(
    const RestartCandidatePackage& value)
{
    RestartManifestExpectation result;
    result.identities = value.manifest.identities;
    result.byteOrder = nativeRestartByteOrder();
    result.cellCount = value.manifest.cellCount;
    result.requirePreviousAcceptedDt = true;
    const auto& roster = restartOwnerRoster();
    for (std::size_t index = 0; index < roster.size(); ++index)
    {
        result.owners[index].owner = roster[index];
        result.owners[index].disposition =
            value.manifest.owners[index].disposition;
        result.owners[index].identitySha256 =
            value.manifest.owners[index].identitySha256;
    }
    return result;
}

RestartGenerationStoreLoadResult loadResult(
    const RestartCandidatePackage& value)
{
    RestartGenerationStoreLoadResult result;
    result.status = RestartGenerationStoreStatus::SUCCESS;
    result.candidate = value;
    return result;
}

RestartGenerationLedgerEntry acceptedLedgerEntry(
    const RestartCandidatePackage& value,
    std::uint64_t ordinal)
{
    RestartGenerationLedgerEntry result;
    result.disposition = RestartGenerationLedgerDisposition::ACCEPTED;
    result.generationId = value.manifest.generationId;
    result.ordinal = ordinal;
    result.generationChecksumSha256 =
        value.manifest.generationChecksumSha256;
    result.status = RestartPublicationStatus::SUCCESS;
    return result;
}

RestartCheckpointClockObservation checkpoint(
    const RestartCandidatePackage& value)
{
    return {
        value.manifest.generationState.physicalTime,
        static_cast<int>(value.manifest.generationState.acceptedIteration),
        1};
}

enum class StageMode
{
    ACCEPT,
    REJECT,
    THROW_EXCEPTION,
    BINDING_MISMATCH,
    MISSING_DOMAIN_LEDGER,
    WRONG_DOMAIN_LEDGER,
    NON_FLUID_DOMAIN_LEDGER
};

struct StageContext
{
    StageMode mode = StageMode::ACCEPT;
    std::uint64_t acceptedSteps = 0;
    std::uint64_t rejectedSteps = 0;
    int calls = 0;
};

class DummyDetachedState final : public RestartDetachedOwnerState
{
public:
    RestartOwner owner;
    std::string generationId;
    std::string payloadSha256;

    DummyDetachedState(
        RestartOwner ownerValue,
        const std::string& generation,
        const std::string& payload,
        std::optional<RestartDomainStepLedger> ledger)
        : RestartDetachedOwnerState(std::move(ledger)),
          owner(ownerValue),generationId(generation),payloadSha256(payload)
    {}
};

RestartOwnerStageResult stageOwner(
    const RestartOwnerStageInput& input,
    void* opaque)
{
    StageContext& context = *static_cast<StageContext*>(opaque);
    ++context.calls;
    if (context.mode == StageMode::THROW_EXCEPTION)
        throw std::runtime_error("injected owner staging exception");

    RestartOwnerStageResult result;
    result.owner = input.owner.owner;
    result.generationId = input.owner.generationId;
    result.identitySha256 = input.owner.identitySha256;
    result.payloadSha256 = input.owner.payloadSha256;
    if (context.mode == StageMode::REJECT)
    {
        result.detail = "injected detached owner rejection";
        return result;
    }
    if (context.mode == StageMode::BINDING_MISMATCH)
        result.identitySha256 = hash("wrong-staged-owner-identity");

    std::optional<RestartDomainStepLedger> ledger;
    if (input.owner.owner == RestartOwner::FLUID
        && context.mode != StageMode::MISSING_DOMAIN_LEDGER)
    {
        ledger = RestartDomainStepLedger{
            context.acceptedSteps,context.rejectedSteps};
        if (context.mode == StageMode::WRONG_DOMAIN_LEDGER)
            ++ledger->acceptedSteps;
    }
    else if (context.mode == StageMode::NON_FLUID_DOMAIN_LEDGER)
        ledger = RestartDomainStepLedger{
            context.acceptedSteps,context.rejectedSteps};

    result.state = std::make_shared<const DummyDetachedState>(
        input.owner.owner,input.owner.generationId,input.owner.payloadSha256,
        std::move(ledger));
    return result;
}

enum class CrossMode
{
    ACCEPT,
    REJECT,
    THROW_EXCEPTION
};

struct CrossContext
{
    CrossMode mode = CrossMode::ACCEPT;
    std::uint64_t expectedAcceptedSteps = 0;
    int calls = 0;
};

RestartCrossValidationResult crossValidate(
    const RestartPublicationRoot& root,
    void* opaque)
{
    CrossContext& context = *static_cast<CrossContext*>(opaque);
    ++context.calls;
    if (context.mode == CrossMode::THROW_EXCEPTION)
        throw std::runtime_error("injected cross-validation exception");
    if (context.mode == CrossMode::REJECT)
        return {false,"injected complete-root rejection"};

    const RestartPublishedOwner& fluid = root.owners()[
        static_cast<std::size_t>(RestartOwner::FLUID)];
    const auto* state = dynamic_cast<const DummyDetachedState*>(
        fluid.state.get());
    if (!root.valid() || state == nullptr
        || state->owner != RestartOwner::FLUID
        || state->generationId != root.manifest().generationId
        || !state->domainStepLedger().has_value()
        || state->domainStepLedger()->acceptedSteps
            != context.expectedAcceptedSteps)
        return {false,"cross-owner Domain/clock binding mismatch"};
    return {true,"complete detached root accepted"};
}

using StageContexts = std::array<StageContext,restartOwnerCount>;

void configureContexts(
    StageContexts& contexts,
    const RestartCandidatePackage& value,
    std::uint64_t rejectedSteps = 3U)
{
    for (StageContext& context : contexts)
    {
        context = StageContext{};
        context.acceptedSteps =
            value.manifest.generationState.acceptedIteration;
        context.rejectedSteps = rejectedSteps;
    }
}

RestartPublicationRequest request(
    const RestartCandidatePackage& value,
    std::uint64_t ordinal,
    StageContexts& contexts,
    CrossContext& cross)
{
    RestartPublicationRequest result;
    result.ledgerEntry = acceptedLedgerEntry(value,ordinal);
    const auto& roster = restartOwnerRoster();
    for (std::size_t index = 0; index < roster.size(); ++index)
    {
        result.ownerStagers[index].owner = roster[index];
        if (value.manifest.owners[index].disposition
            == RestartOwnerDisposition::STATE_PAYLOAD)
        {
            result.ownerStagers[index].stage = &stageOwner;
            result.ownerStagers[index].context = &contexts[index];
        }
    }
    result.crossValidate = &crossValidate;
    result.crossValidationContext = &cross;
    return result;
}

RestartPublicationStageResult stage(
    const std::shared_ptr<const RestartPublicationRoot>& base,
    const RestartCandidatePackage& value,
    RestartPublicationRequest requestValue)
{
    return RestartPublicationStager::stage(
        base,loadResult(value),expectation(value),checkpoint(value),
        requestValue);
}

void checkLiveUnchanged(
    const std::shared_ptr<const RestartPublicationRoot>& live,
    const std::shared_ptr<const RestartPublicationRoot>& before,
    const std::string& label)
{
    check(live.get() == before.get(),label + " keeps the live root pointer");
    if (before)
    {
        check(live->manifest().generationId
                == before->manifest().generationId
                && live->acceptedLedger().size()
                    == before->acceptedLedger().size()
                && live->owners()[0].state.get()
                    == before->owners()[0].state.get(),
            label + " keeps the immutable live-root contents");
    }
}

void checkStageFailure(
    RestartPublicationStageResult& result,
    RestartPublicationStatus expected,
    const std::shared_ptr<const RestartPublicationRoot>& live,
    const std::shared_ptr<const RestartPublicationRoot>& before,
    const std::string& label)
{
    check(result.status == expected && !result.success()
            && !result.publication,
        label + " returns its structured failure");
    check(result.rejectedLedgerEntry.has_value()
            && result.rejectedLedgerEntry->rejected()
            && result.rejectedLedgerEntry->status == expected,
        label + " emits only a diagnostic rejected-ledger entry");
    checkLiveUnchanged(live,before,label);
}

} // End anonymous namespace

int main()
{
    try
    {
        check(std::string(restartPublicationStatusName(
                RestartPublicationStatus::OWNER_LEDGER_MISMATCH))
                == "OWNER_LEDGER_MISMATCH",
            "Domain ledger mismatch has a stable status name");

        const RestartCandidatePackage first =
            candidate("generation-publication-1",1.0e-6,10U);
        StageContexts firstContexts;
        configureContexts(firstContexts,first,2U);
        CrossContext firstCross;
        firstCross.expectedAcceptedSteps = 10U;
        RestartPublicationRequest firstRequest =
            request(first,1U,firstContexts,firstCross);
        check(validRestartGenerationLedgerEntry(firstRequest.ledgerEntry),
            "first accepted ledger entry is valid");

        std::shared_ptr<const RestartPublicationRoot> live;
        RestartPublicationStageResult firstStage =
            stage(live,first,firstRequest);
        check(firstStage.success() && firstStage.publication->candidate().valid(),
            "first complete detached root stages from an empty base");
        check(firstCross.calls == 1
                && firstStage.publication->candidate().clock()
                    .acceptedIteration == 10,
            "cross-validation sees the coordinator clock before publication");
        const auto& firstDomainLedger =
            firstStage.publication->candidate().owners()[0]
                .state->domainStepLedger();
        check(firstDomainLedger.has_value()
                && firstDomainLedger->acceptedSteps == 10U
                && firstDomainLedger->rejectedSteps == 2U,
            "fluid Domain accepted/rejected ledger is staged immutably");

        const RestartPublicationRoot* firstCandidateAddress =
            &firstStage.publication->candidate();
        const RestartPublicationCommitResult firstCommit =
            firstStage.publication->publish(live);
        check(firstCommit.success() && live.get() == firstCandidateAddress,
            "initial publication installs the exact staged composite root");
        check(live->acceptedLedger().size() == 1U
                && live->acceptedLedger().back().ordinal == 1U,
            "accepted ledger becomes visible with the root swap");
        const std::shared_ptr<const RestartPublicationRoot> firstLive = live;
        const RestartPublicationCommitResult consumed =
            firstStage.publication->publish(live);
        check(consumed.status
                == RestartPublicationStatus::PUBLICATION_ALREADY_CONSUMED,
            "a successful staged publication is one-shot");
        checkLiveUnchanged(live,firstLive,"consumed publication retry");

        const RestartCandidatePackage slow =
            candidate("generation-publication-slow",2.0e-6,11U);
        const RestartCandidatePackage winner =
            candidate("generation-publication-winner",2.0e-6,11U);
        StageContexts slowContexts, winnerContexts;
        configureContexts(slowContexts,slow,4U);
        configureContexts(winnerContexts,winner,4U);
        CrossContext slowCross, winnerCross;
        slowCross.expectedAcceptedSteps = 11U;
        winnerCross.expectedAcceptedSteps = 11U;
        RestartPublicationStageResult slowStage = stage(
            live,slow,request(slow,2U,slowContexts,slowCross));
        RestartPublicationStageResult winnerStage = stage(
            live,winner,request(winner,2U,winnerContexts,winnerCross));
        check(slowStage.success() && winnerStage.success(),
            "two detached candidates may stage against one immutable base");
        const RestartPublicationRoot* winnerAddress =
            &winnerStage.publication->candidate();
        check(winnerStage.publication->publish(live).success()
                && live.get() == winnerAddress,
            "winner publishes by one composite-root pointer swap");
        const std::shared_ptr<const RestartPublicationRoot> winnerLive = live;
        const RestartPublicationCommitResult staleCommit =
            slowStage.publication->publish(live);
        check(staleCommit.status
                == RestartPublicationStatus::STALE_PUBLICATION_BASE
                && staleCommit.rejectedLedgerEntry.has_value(),
            "staged candidate rejects a stale publication base");
        checkLiveUnchanged(live,winnerLive,"stale-base rejection");
        check(firstLive->manifest().generationId
                == "generation-publication-1"
                && firstLive->acceptedLedger().size() == 1U,
            "older immutable root remains readable and unchanged");

        const RestartCandidatePackage conflict =
            candidate("generation-publication-conflict",3.0e-6,12U);
        StageContexts conflictContexts;
        configureContexts(conflictContexts,conflict);
        CrossContext conflictCross;
        conflictCross.expectedAcceptedSteps = 12U;
        RestartPublicationStageResult duplicateOrdinal = stage(
            live,conflict,request(
                conflict,2U,conflictContexts,conflictCross));
        checkStageFailure(duplicateOrdinal,
            RestartPublicationStatus::LEDGER_CONFLICT,
            live,winnerLive,"duplicate accepted-ledger ordinal");

        StageContexts duplicateIdContexts;
        configureContexts(duplicateIdContexts,winner);
        CrossContext duplicateIdCross;
        duplicateIdCross.expectedAcceptedSteps = 11U;
        RestartPublicationStageResult duplicateId = stage(
            live,winner,request(
                winner,3U,duplicateIdContexts,duplicateIdCross));
        checkStageFailure(duplicateId,
            RestartPublicationStatus::LEDGER_CONFLICT,
            live,winnerLive,"duplicate accepted generation ID");

        StageContexts rejectedContexts;
        configureContexts(rejectedContexts,conflict);
        CrossContext rejectedCross;
        RestartPublicationRequest rejectedRequest = request(
            conflict,3U,rejectedContexts,rejectedCross);
        rejectedRequest.ledgerEntry.disposition =
            RestartGenerationLedgerDisposition::REJECTED;
        rejectedRequest.ledgerEntry.status =
            RestartPublicationStatus::MANIFEST_REJECTED;
        rejectedRequest.ledgerEntry.detail = "prior offline rejection";
        RestartPublicationStageResult rejectedSelection = stage(
            live,conflict,rejectedRequest);
        checkStageFailure(rejectedSelection,
            RestartPublicationStatus::LEDGER_REJECTED,
            live,winnerLive,"rejected generation selection");

        RestartGenerationLedgerEntry invalidStatus =
            rejectedRequest.ledgerEntry;
        invalidStatus.status = static_cast<RestartPublicationStatus>(255);
        check(!validRestartGenerationLedgerEntry(invalidStatus),
            "unknown rejected-ledger status fails closed");
        std::vector<RestartGenerationLedgerEntry> duplicatedAccepted =
            live->acceptedLedger();
        duplicatedAccepted.push_back(duplicatedAccepted.back());
        check(!validRestartAcceptedGenerationLedger(duplicatedAccepted),
            "accepted ledger rejects duplicate generation and ordinal");

        StageContexts exceptionContexts;
        configureContexts(exceptionContexts,conflict);
        exceptionContexts[static_cast<std::size_t>(RestartOwner::TURBULENCE)]
            .mode = StageMode::THROW_EXCEPTION;
        CrossContext exceptionCross;
        exceptionCross.expectedAcceptedSteps = 12U;
        RestartPublicationStageResult ownerException = stage(
            live,conflict,request(
                conflict,3U,exceptionContexts,exceptionCross));
        checkStageFailure(ownerException,
            RestartPublicationStatus::OWNER_STAGE_EXCEPTION,
            live,winnerLive,"owner staging exception");
        check(ownerException.rejectedLedgerEntry->owner
                == std::optional<RestartOwner>(RestartOwner::TURBULENCE),
            "owner staging exception records the failing owner");

        StageContexts bindingContexts;
        configureContexts(bindingContexts,conflict);
        bindingContexts[static_cast<std::size_t>(RestartOwner::TURBULENCE)]
            .mode = StageMode::BINDING_MISMATCH;
        CrossContext bindingCross;
        bindingCross.expectedAcceptedSteps = 12U;
        RestartPublicationStageResult bindingMismatch = stage(
            live,conflict,request(
                conflict,3U,bindingContexts,bindingCross));
        checkStageFailure(bindingMismatch,
            RestartPublicationStatus::OWNER_BINDING_MISMATCH,
            live,winnerLive,"owner identity binding mismatch");

        StageContexts clockLedgerContexts;
        configureContexts(clockLedgerContexts,conflict);
        clockLedgerContexts[static_cast<std::size_t>(RestartOwner::FLUID)]
            .mode = StageMode::WRONG_DOMAIN_LEDGER;
        CrossContext clockLedgerCross;
        clockLedgerCross.expectedAcceptedSteps = 12U;
        RestartPublicationStageResult clockLedgerMismatch = stage(
            live,conflict,request(
                conflict,3U,clockLedgerContexts,clockLedgerCross));
        checkStageFailure(clockLedgerMismatch,
            RestartPublicationStatus::OWNER_LEDGER_MISMATCH,
            live,winnerLive,"coordinator clock and Domain ledger mismatch");
        check(clockLedgerCross.calls == 0,
            "Domain ledger mismatch rejects before cross-validation");

        StageContexts missingLedgerContexts;
        configureContexts(missingLedgerContexts,conflict);
        missingLedgerContexts[static_cast<std::size_t>(RestartOwner::FLUID)]
            .mode = StageMode::MISSING_DOMAIN_LEDGER;
        CrossContext missingLedgerCross;
        missingLedgerCross.expectedAcceptedSteps = 12U;
        RestartPublicationStageResult missingLedger = stage(
            live,conflict,request(
                conflict,3U,missingLedgerContexts,missingLedgerCross));
        checkStageFailure(missingLedger,
            RestartPublicationStatus::OWNER_LEDGER_MISMATCH,
            live,winnerLive,"missing fluid Domain ledger");

        StageContexts foreignLedgerContexts;
        configureContexts(foreignLedgerContexts,conflict);
        foreignLedgerContexts[
            static_cast<std::size_t>(RestartOwner::TURBULENCE)]
                .mode = StageMode::NON_FLUID_DOMAIN_LEDGER;
        CrossContext foreignLedgerCross;
        foreignLedgerCross.expectedAcceptedSteps = 12U;
        RestartPublicationStageResult foreignLedger = stage(
            live,conflict,request(
                conflict,3U,foreignLedgerContexts,foreignLedgerCross));
        checkStageFailure(foreignLedger,
            RestartPublicationStatus::OWNER_LEDGER_MISMATCH,
            live,winnerLive,"non-fluid Domain ledger injection");

        StageContexts crossRejectContexts;
        configureContexts(crossRejectContexts,conflict);
        CrossContext crossReject;
        crossReject.mode = CrossMode::REJECT;
        crossReject.expectedAcceptedSteps = 12U;
        RestartPublicationStageResult crossRejected = stage(
            live,conflict,request(
                conflict,3U,crossRejectContexts,crossReject));
        checkStageFailure(crossRejected,
            RestartPublicationStatus::CROSS_VALIDATION_REJECTED,
            live,winnerLive,"complete-root cross-validation rejection");
        check(crossReject.calls == 1,
            "cross-validation rejection occurs after owner staging");

        StageContexts crossThrowContexts;
        configureContexts(crossThrowContexts,conflict);
        CrossContext crossThrow;
        crossThrow.mode = CrossMode::THROW_EXCEPTION;
        RestartPublicationStageResult crossException = stage(
            live,conflict,request(
                conflict,3U,crossThrowContexts,crossThrow));
        checkStageFailure(crossException,
            RestartPublicationStatus::CROSS_VALIDATION_EXCEPTION,
            live,winnerLive,"complete-root cross-validation exception");

        StageContexts finalContexts;
        configureContexts(finalContexts,conflict,9U);
        CrossContext finalCross;
        finalCross.expectedAcceptedSteps = 12U;
        RestartPublicationStageResult finalStage = stage(
            live,conflict,request(
                conflict,3U,finalContexts,finalCross));
        check(finalStage.success(),
            "fully consistent successor stages after rejected attempts");
        const RestartPublicationRoot* finalAddress =
            &finalStage.publication->candidate();
        const RestartPublicationCommitResult finalCommit =
            finalStage.publication->publish(live);
        check(finalCommit.success() && live.get() == finalAddress
                && live.get() != winnerLive.get(),
            "success performs one exact successor-root pointer swap");
        check(live->manifest().generationId
                == "generation-publication-conflict"
                && live->clock().acceptedIteration == 12
                && live->acceptedLedger().size() == 3U
                && live->acceptedLedger().back().ordinal == 3U,
            "generation clock and accepted ledger publish together");
        const auto& finalDomainLedger =
            live->owners()[static_cast<std::size_t>(RestartOwner::FLUID)]
                .state->domainStepLedger();
        check(finalDomainLedger.has_value()
                && finalDomainLedger->acceptedSteps
                    == static_cast<std::uint64_t>(
                        live->clock().acceptedIteration)
                && finalDomainLedger->rejectedSteps == 9U,
            "published Domain ledger remains consistent with the clock");
        check(winnerLive->manifest().generationId
                == "generation-publication-winner"
                && winnerLive->acceptedLedger().size() == 2U,
            "successful successor publication does not mutate the prior root");

        RestartRejectedGenerationLedger diagnostics;
        check(diagnostics.record(*ownerException.rejectedLedgerEntry)
                && diagnostics.record(*crossRejected.rejectedLedgerEntry)
                && diagnostics.entries().size() == 2U
                && live->acceptedLedger().size() == 3U,
            "diagnostic rejected ledger cannot alter selection authority");

        std::cout
            << "RESTART_PUBLICATION_CANDIDATE_PASS checks=" << checks
            << " scope=offline-immutable-composite-root"
            << " solver_integration=0 production_ready=0\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "RESTART_PUBLICATION_CANDIDATE_FAIL "
                  << error.what() << '\n';
        return 1;
    }
}
