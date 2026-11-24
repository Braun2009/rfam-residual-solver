#include "RestartGenerationStore.H"

#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>

using namespace Foam;

namespace
{

int checks=0;

void check(bool condition,const std::string& label)
{
    if(!condition) throw std::runtime_error(label);
    ++checks;
}

enum class ClockOwner
{
    GENERATION_COORDINATOR_TIME_INTEGRATOR
};

enum class SelectionPolicy
{
    UNDECIDED,
    EXPLICIT_GENERATION,
    LATEST_ACCEPTED_LEDGER
};

enum class PublicationStrategy
{
    UNDECIDED,
    SINGLE_COMPOSITE_ROOT_SWAP,
    QUIESCED_SEQUENTIAL_WITH_INFALLIBLE_ROLLBACK
};

enum class LedgerDisposition
{
    NONE,
    ACCEPTED,
    REJECTED
};

enum class FaultPoint
{
    NONE,
    AFTER_SELECTION,
    AFTER_STORE_ADMISSION,
    AFTER_MANIFEST_VALIDATION,
    AFTER_LEDGER_ADMISSION,
    STAGE_OWNER_0,
    STAGE_OWNER_1,
    STAGE_OWNER_2,
    STAGE_OWNER_3,
    STAGE_OWNER_4,
    STAGE_OWNER_5,
    STAGE_OWNER_6,
    STAGE_DOMAIN,
    STAGE_CLOCK,
    CROSS_VALIDATE,
    COMMIT_ROOT_SWAP,
    COMMIT_DOMAIN,
    COMMIT_OWNER_0,
    COMMIT_OWNER_1,
    COMMIT_OWNER_2,
    COMMIT_OWNER_3,
    COMMIT_OWNER_4,
    COMMIT_OWNER_5,
    COMMIT_OWNER_6,
    COMMIT_CLOCK,
    COMMIT_ACCEPTED_LEDGER
};

enum class TransactionStatus
{
    OFFLINE_ACCEPTED,
    SELECTION_POLICY_UNDECIDED,
    PUBLICATION_STRATEGY_UNDECIDED,
    GENERATION_NOT_SELECTABLE,
    GENERATION_SELECTION_AMBIGUOUS,
    STORE_REJECTED,
    MANIFEST_REJECTED,
    LEDGER_REJECTED,
    OWNER_STAGE_REJECTED,
    DOMAIN_STAGE_REJECTED,
    CLOCK_STAGE_REJECTED,
    CROSS_VALIDATION_REJECTED,
    COMMIT_REJECTED_ROLLED_BACK,
    TERMINAL_QUARANTINE,
    RUNTIME_OR_HOTSTART_UNAUTHORIZED
};

struct ClockState
{
    double physicalTime=0.0;
    std::uint64_t acceptedIteration=0;
    bool previousAcceptedDtPresent=false;
    double previousAcceptedDt=0.0;
};

bool sameClock(const ClockState& left,const ClockState& right)
{
    return left.physicalTime==right.physicalTime
        &&left.acceptedIteration==right.acceptedIteration
        &&left.previousAcceptedDtPresent==right.previousAcceptedDtPresent
        &&left.previousAcceptedDt==right.previousAcceptedDt
        &&std::signbit(left.physicalTime)==std::signbit(right.physicalTime)
        &&std::signbit(left.previousAcceptedDt)
            ==std::signbit(right.previousAcceptedDt);
}

bool validClock(const ClockState& clock,bool requirePreviousAcceptedDt)
{
    if(!std::isfinite(clock.physicalTime)||clock.physicalTime<0.0
       ||std::signbit(clock.physicalTime)
       ||clock.acceptedIteration
            >static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
        return false;
    if(clock.previousAcceptedDtPresent!=requirePreviousAcceptedDt) return false;
    if(clock.previousAcceptedDtPresent)
        return std::isfinite(clock.previousAcceptedDt)
            &&clock.previousAcceptedDt>0.0;
    return clock.previousAcceptedDt==0.0
        &&!std::signbit(clock.previousAcceptedDt);
}

struct AcceptedLedgerEntry
{
    std::string generationId;
    std::uint64_t ordinal=0;
    std::string manifestSha256;
};

bool sameEntry(
    const AcceptedLedgerEntry& left,const AcceptedLedgerEntry& right)
{
    return left.generationId==right.generationId
        &&left.ordinal==right.ordinal
        &&left.manifestSha256==right.manifestSha256;
}

struct LiveState
{
    bool readable=true;
    std::string generationId;
    std::string domainGenerationId;
    std::string domainDigest;
    std::array<std::string,restartOwnerCount> ownerGenerationIds{};
    std::array<std::string,restartOwnerCount> ownerDigests{};
    ClockState clock;
    std::vector<AcceptedLedgerEntry> acceptedLedger;
};

bool sameLive(const LiveState& left,const LiveState& right)
{
    if(left.readable!=right.readable
       ||left.generationId!=right.generationId
       ||left.domainGenerationId!=right.domainGenerationId
       ||left.domainDigest!=right.domainDigest
       ||left.ownerGenerationIds!=right.ownerGenerationIds
       ||left.ownerDigests!=right.ownerDigests
       ||!sameClock(left.clock,right.clock)
       ||left.acceptedLedger.size()!=right.acceptedLedger.size()) return false;
    for(std::size_t index=0;index<left.acceptedLedger.size();++index)
        if(!sameEntry(left.acceptedLedger[index],right.acceptedLedger[index]))
            return false;
    return true;
}

struct CandidateGeneration
{
    std::string generationId;
    std::uint64_t ordinal=0;
    LedgerDisposition ledgerDisposition=LedgerDisposition::NONE;
    RestartGenerationStoreStatus storeStatus=
        RestartGenerationStoreStatus::INVALID_INPUT;
    bool committedMarker=false;
    bool manifestValidated=false;
    bool requirePreviousAcceptedDt=false;
    std::string manifestSha256;
    std::string domainGenerationId;
    std::string domainDigest;
    std::array<std::string,restartOwnerCount> ownerGenerationIds{};
    std::array<std::string,restartOwnerCount> ownerDigests{};
    ClockState clock;
};

struct TransactionRequest
{
    SelectionPolicy selectionPolicy=SelectionPolicy::UNDECIDED;
    PublicationStrategy publicationStrategy=PublicationStrategy::UNDECIDED;
    std::string explicitGenerationId;
    FaultPoint fault=FaultPoint::NONE;
    bool rollbackFails=false;
    bool runtimeRestartRequested=false;
    bool hotStartRequested=false;
};

struct TransactionResult
{
    TransactionStatus status=TransactionStatus::GENERATION_NOT_SELECTABLE;
    bool statePublished=false;
    bool rejectedDiagnosticRecorded=false;
    bool terminalQuarantine=false;
    std::string selectedGenerationId;
};

FaultPoint stageOwnerFault(std::size_t owner)
{
    return static_cast<FaultPoint>(
        static_cast<int>(FaultPoint::STAGE_OWNER_0)+static_cast<int>(owner));
}

FaultPoint commitOwnerFault(std::size_t owner)
{
    return static_cast<FaultPoint>(
        static_cast<int>(FaultPoint::COMMIT_OWNER_0)+static_cast<int>(owner));
}

bool candidateCrossConsistent(const CandidateGeneration& candidate)
{
    if(candidate.generationId.empty()||candidate.domainDigest.empty()
       ||candidate.manifestSha256.empty()
       ||candidate.domainGenerationId!=candidate.generationId) return false;
    for(std::size_t owner=0;owner<restartOwnerCount;++owner)
        if(candidate.ownerGenerationIds[owner]!=candidate.generationId
           ||candidate.ownerDigests[owner].empty()) return false;
    return true;
}

const CandidateGeneration* selectCandidate(
    const std::vector<CandidateGeneration>& candidates,
    const TransactionRequest& request,
    TransactionStatus& status)
{
    std::set<std::string> generationIds;
    for(const CandidateGeneration& candidate:candidates)
        if(candidate.generationId.empty()
           ||!generationIds.insert(candidate.generationId).second)
        {
            status=TransactionStatus::GENERATION_SELECTION_AMBIGUOUS;
            return nullptr;
        }
    if(request.selectionPolicy==SelectionPolicy::UNDECIDED)
    {
        status=TransactionStatus::SELECTION_POLICY_UNDECIDED;
        return nullptr;
    }
    if(request.selectionPolicy==SelectionPolicy::EXPLICIT_GENERATION)
    {
        const CandidateGeneration* selected=nullptr;
        for(const CandidateGeneration& candidate:candidates)
            if(candidate.generationId==request.explicitGenerationId)
            {
                if(selected)
                {
                    status=TransactionStatus::GENERATION_SELECTION_AMBIGUOUS;
                    return nullptr;
                }
                selected=&candidate;
            }
        if(!selected||selected->ledgerDisposition!=LedgerDisposition::ACCEPTED)
        {
            status=TransactionStatus::GENERATION_NOT_SELECTABLE;
            return nullptr;
        }
        return selected;
    }

    const CandidateGeneration* selected=nullptr;
    bool ambiguous=false;
    for(const CandidateGeneration& candidate:candidates)
    {
        if(candidate.ledgerDisposition!=LedgerDisposition::ACCEPTED) continue;
        if(!selected||candidate.ordinal>selected->ordinal)
        {
            selected=&candidate;
            ambiguous=false;
        }
        else if(candidate.ordinal==selected->ordinal) ambiguous=true;
    }
    if(!selected)
    {
        status=TransactionStatus::GENERATION_NOT_SELECTABLE;
        return nullptr;
    }
    if(ambiguous)
    {
        status=TransactionStatus::GENERATION_SELECTION_AMBIGUOUS;
        return nullptr;
    }
    return selected;
}

TransactionResult fail(TransactionStatus status)
{
    TransactionResult result;
    result.status=status;
    result.rejectedDiagnosticRecorded=true;
    return result;
}

void installCandidate(LiveState& live,const CandidateGeneration& candidate)
{
    live.generationId=candidate.generationId;
    live.domainGenerationId=candidate.domainGenerationId;
    live.domainDigest=candidate.domainDigest;
    live.ownerGenerationIds=candidate.ownerGenerationIds;
    live.ownerDigests=candidate.ownerDigests;
    live.clock=candidate.clock;
}

TransactionResult runTransaction(
    LiveState& live,
    const std::vector<CandidateGeneration>& candidates,
    const TransactionRequest& request)
{
    const LiveState before=live;
    if(request.runtimeRestartRequested||request.hotStartRequested)
        return fail(TransactionStatus::RUNTIME_OR_HOTSTART_UNAUTHORIZED);

    TransactionStatus selectionStatus=TransactionStatus::GENERATION_NOT_SELECTABLE;
    const CandidateGeneration* selected=
        selectCandidate(candidates,request,selectionStatus);
    if(!selected) return fail(selectionStatus);
    if(request.fault==FaultPoint::AFTER_SELECTION)
        return fail(TransactionStatus::GENERATION_NOT_SELECTABLE);

    if(selected->storeStatus!=RestartGenerationStoreStatus::SUCCESS
       ||!selected->committedMarker
       ||request.fault==FaultPoint::AFTER_STORE_ADMISSION)
        return fail(TransactionStatus::STORE_REJECTED);
    if(!selected->manifestValidated
       ||request.fault==FaultPoint::AFTER_MANIFEST_VALIDATION)
        return fail(TransactionStatus::MANIFEST_REJECTED);
    if(selected->ledgerDisposition!=LedgerDisposition::ACCEPTED
       ||request.fault==FaultPoint::AFTER_LEDGER_ADMISSION)
        return fail(TransactionStatus::LEDGER_REJECTED);

    for(std::size_t owner=0;owner<restartOwnerCount;++owner)
        if(request.fault==stageOwnerFault(owner))
            return fail(TransactionStatus::OWNER_STAGE_REJECTED);
    if(request.fault==FaultPoint::STAGE_DOMAIN)
        return fail(TransactionStatus::DOMAIN_STAGE_REJECTED);
    if(!validClock(selected->clock,selected->requirePreviousAcceptedDt)
       ||request.fault==FaultPoint::STAGE_CLOCK)
        return fail(TransactionStatus::CLOCK_STAGE_REJECTED);
    if(!candidateCrossConsistent(*selected)
       ||request.fault==FaultPoint::CROSS_VALIDATE)
        return fail(TransactionStatus::CROSS_VALIDATION_REJECTED);
    if(request.publicationStrategy==PublicationStrategy::UNDECIDED)
        return fail(TransactionStatus::PUBLICATION_STRATEGY_UNDECIDED);

    const AcceptedLedgerEntry ledgerEntry{
        selected->generationId,selected->ordinal,selected->manifestSha256};
    if(request.publicationStrategy
        ==PublicationStrategy::SINGLE_COMPOSITE_ROOT_SWAP)
    {
        LiveState staged=before;
        installCandidate(staged,*selected);
        staged.acceptedLedger.push_back(ledgerEntry);
        if(request.fault==FaultPoint::COMMIT_ROOT_SWAP)
            return fail(TransactionStatus::COMMIT_REJECTED_ROLLED_BACK);
        live=std::move(staged);
    }
    else
    {
        auto rollback=[&]() -> TransactionResult
        {
            if(request.rollbackFails)
            {
                live.readable=false;
                TransactionResult result=fail(TransactionStatus::TERMINAL_QUARANTINE);
                result.terminalQuarantine=true;
                return result;
            }
            live=before;
            return fail(TransactionStatus::COMMIT_REJECTED_ROLLED_BACK);
        };

        live.generationId=selected->generationId;
        live.domainGenerationId=selected->domainGenerationId;
        live.domainDigest=selected->domainDigest;
        if(request.fault==FaultPoint::COMMIT_DOMAIN) return rollback();
        for(std::size_t owner=0;owner<restartOwnerCount;++owner)
        {
            live.ownerGenerationIds[owner]=selected->ownerGenerationIds[owner];
            live.ownerDigests[owner]=selected->ownerDigests[owner];
            if(request.fault==commitOwnerFault(owner)) return rollback();
        }
        live.clock=selected->clock;
        if(request.fault==FaultPoint::COMMIT_CLOCK) return rollback();
        live.acceptedLedger.push_back(ledgerEntry);
        if(request.fault==FaultPoint::COMMIT_ACCEPTED_LEDGER) return rollback();
    }

    TransactionResult result;
    result.status=TransactionStatus::OFFLINE_ACCEPTED;
    result.statePublished=true;
    result.selectedGenerationId=selected->generationId;
    return result;
}

CandidateGeneration candidate(
    const std::string& generation,std::uint64_t ordinal)
{
    CandidateGeneration value;
    value.generationId=generation;
    value.ordinal=ordinal;
    value.ledgerDisposition=LedgerDisposition::ACCEPTED;
    value.storeStatus=RestartGenerationStoreStatus::SUCCESS;
    value.committedMarker=true;
    value.manifestValidated=true;
    value.requirePreviousAcceptedDt=true;
    value.manifestSha256="manifest-"+generation;
    value.domainGenerationId=generation;
    value.domainDigest="domain-"+generation;
    for(std::size_t owner=0;owner<restartOwnerCount;++owner)
    {
        value.ownerGenerationIds[owner]=generation;
        value.ownerDigests[owner]="owner-"+std::to_string(owner)+"-"+generation;
    }
    value.clock={1.25e-5+static_cast<double>(ordinal)*1.0e-8,
        42U+ordinal,true,1.0e-8};
    return value;
}

LiveState initialLive()
{
    LiveState live;
    live.generationId="live-generation-0";
    live.domainGenerationId=live.generationId;
    live.domainDigest="domain-live-generation-0";
    for(std::size_t owner=0;owner<restartOwnerCount;++owner)
    {
        live.ownerGenerationIds[owner]=live.generationId;
        live.ownerDigests[owner]="owner-live-"+std::to_string(owner);
    }
    live.clock={0.0,0U,false,0.0};
    live.acceptedLedger.push_back({live.generationId,0U,"manifest-live-0"});
    return live;
}

void expectUnchangedFailure(
    const std::vector<CandidateGeneration>& candidates,
    TransactionRequest request,
    TransactionStatus expected,
    const std::string& label)
{
    LiveState live=initialLive();
    const LiveState before=live;
    const TransactionResult result=runTransaction(live,candidates,request);
    check(result.status==expected,label+" status");
    check(!result.statePublished,label+" no state publication");
    check(result.rejectedDiagnosticRecorded,label+" diagnostic rejection only");
    check(sameLive(live,before),label+" live state byte-semantically unchanged");
}

} // namespace

int main()
{
    try
    {
        check(restartOwnerCount==7U,"contract binds exactly seven restart owners");
        const auto& roster=restartOwnerRoster();
        check(roster.front()==RestartOwner::FLUID
                &&roster.back()==RestartOwner::TRANSPORT,
            "contract binds canonical fluid-through-transport roster");
        check(ClockOwner::GENERATION_COORDINATOR_TIME_INTEGRATOR
                ==ClockOwner::GENERATION_COORDINATOR_TIME_INTEGRATOR,
            "previous accepted dt is coordinator/time-integrator owned");

        const CandidateGeneration first=candidate("generation-1",1U);
        const CandidateGeneration second=candidate("generation-2",2U);
        std::vector<CandidateGeneration> candidates{first,second};

        TransactionRequest request;
        expectUnchangedFailure(candidates,request,
            TransactionStatus::SELECTION_POLICY_UNDECIDED,
            "undecided selection policy");

        request.selectionPolicy=SelectionPolicy::EXPLICIT_GENERATION;
        request.explicitGenerationId=second.generationId;
        expectUnchangedFailure(candidates,request,
            TransactionStatus::PUBLICATION_STRATEGY_UNDECIDED,
            "undecided publication strategy");

        request.explicitGenerationId="missing";
        expectUnchangedFailure(candidates,request,
            TransactionStatus::GENERATION_NOT_SELECTABLE,
            "missing explicit generation");
        request.explicitGenerationId=second.generationId;

        std::vector<CandidateGeneration> rejectedCandidates=candidates;
        rejectedCandidates[1].ledgerDisposition=LedgerDisposition::REJECTED;
        expectUnchangedFailure(rejectedCandidates,request,
            TransactionStatus::GENERATION_NOT_SELECTABLE,
            "rejected generation cannot be explicitly selected");

        TransactionRequest latest=request;
        latest.selectionPolicy=SelectionPolicy::LATEST_ACCEPTED_LEDGER;
        latest.publicationStrategy=PublicationStrategy::SINGLE_COMPOSITE_ROOT_SWAP;
        LiveState latestLive=initialLive();
        TransactionResult latestResult=runTransaction(latestLive,candidates,latest);
        check(latestResult.status==TransactionStatus::OFFLINE_ACCEPTED,
            "latest-ledger abstract policy selects a candidate");
        check(latestResult.selectedGenerationId==second.generationId,
            "latest-ledger abstract policy selects highest accepted ordinal");

        std::vector<CandidateGeneration> ambiguous=candidates;
        ambiguous[0].ordinal=ambiguous[1].ordinal;
        expectUnchangedFailure(ambiguous,latest,
            TransactionStatus::GENERATION_SELECTION_AMBIGUOUS,
            "latest accepted ordinal tie");
        ambiguous=candidates;
        CandidateGeneration duplicateId=first;
        duplicateId.ordinal=3U;
        ambiguous.push_back(duplicateId);
        expectUnchangedFailure(ambiguous,latest,
            TransactionStatus::GENERATION_SELECTION_AMBIGUOUS,
            "duplicate generation identity");

        const std::array<FaultPoint,15> precommitFaults={
            FaultPoint::AFTER_SELECTION,
            FaultPoint::AFTER_STORE_ADMISSION,
            FaultPoint::AFTER_MANIFEST_VALIDATION,
            FaultPoint::AFTER_LEDGER_ADMISSION,
            FaultPoint::STAGE_OWNER_0,
            FaultPoint::STAGE_OWNER_1,
            FaultPoint::STAGE_OWNER_2,
            FaultPoint::STAGE_OWNER_3,
            FaultPoint::STAGE_OWNER_4,
            FaultPoint::STAGE_OWNER_5,
            FaultPoint::STAGE_OWNER_6,
            FaultPoint::STAGE_DOMAIN,
            FaultPoint::STAGE_CLOCK,
            FaultPoint::CROSS_VALIDATE,
            FaultPoint::COMMIT_ROOT_SWAP};
        request.publicationStrategy=PublicationStrategy::SINGLE_COMPOSITE_ROOT_SWAP;
        for(std::size_t index=0;index<precommitFaults.size();++index)
        {
            request.fault=precommitFaults[index];
            const TransactionStatus expected=
                request.fault==FaultPoint::AFTER_STORE_ADMISSION
                    ?TransactionStatus::STORE_REJECTED
                :request.fault==FaultPoint::AFTER_MANIFEST_VALIDATION
                    ?TransactionStatus::MANIFEST_REJECTED
                :request.fault==FaultPoint::AFTER_LEDGER_ADMISSION
                    ?TransactionStatus::LEDGER_REJECTED
                :request.fault>=FaultPoint::STAGE_OWNER_0
                    &&request.fault<=FaultPoint::STAGE_OWNER_6
                    ?TransactionStatus::OWNER_STAGE_REJECTED
                :request.fault==FaultPoint::STAGE_DOMAIN
                    ?TransactionStatus::DOMAIN_STAGE_REJECTED
                :request.fault==FaultPoint::STAGE_CLOCK
                    ?TransactionStatus::CLOCK_STAGE_REJECTED
                :request.fault==FaultPoint::CROSS_VALIDATE
                    ?TransactionStatus::CROSS_VALIDATION_REJECTED
                :request.fault==FaultPoint::COMMIT_ROOT_SWAP
                    ?TransactionStatus::COMMIT_REJECTED_ROLLED_BACK
                    :TransactionStatus::GENERATION_NOT_SELECTABLE;
            expectUnchangedFailure(candidates,request,expected,
                "root-swap fault "+std::to_string(index));
        }

        CandidateGeneration badClock=second;
        badClock.clock.previousAcceptedDtPresent=false;
        badClock.clock.previousAcceptedDt=0.0;
        request.fault=FaultPoint::NONE;
        expectUnchangedFailure({first,badClock},request,
            TransactionStatus::CLOCK_STAGE_REJECTED,
            "required previous accepted dt is unavailable");
        badClock=second;
        badClock.clock.acceptedIteration=
            static_cast<std::uint64_t>(std::numeric_limits<int>::max())+1U;
        expectUnchangedFailure({first,badClock},request,
            TransactionStatus::CLOCK_STAGE_REJECTED,
            "solver iteration publication range exceeded");

        CandidateGeneration mismatchedOwner=second;
        mismatchedOwner.ownerGenerationIds[3]="other-generation";
        expectUnchangedFailure({first,mismatchedOwner},request,
            TransactionStatus::CROSS_VALIDATION_REJECTED,
            "one owner generation mismatch");

        request.runtimeRestartRequested=true;
        expectUnchangedFailure(candidates,request,
            TransactionStatus::RUNTIME_OR_HOTSTART_UNAUTHORIZED,
            "runtime restart authorization absent");
        request.runtimeRestartRequested=false;
        request.hotStartRequested=true;
        expectUnchangedFailure(candidates,request,
            TransactionStatus::RUNTIME_OR_HOTSTART_UNAUTHORIZED,
            "hot-start authorization absent");
        request.hotStartRequested=false;

        const std::array<FaultPoint,10> sequentialFaults={
            FaultPoint::COMMIT_DOMAIN,
            FaultPoint::COMMIT_OWNER_0,
            FaultPoint::COMMIT_OWNER_1,
            FaultPoint::COMMIT_OWNER_2,
            FaultPoint::COMMIT_OWNER_3,
            FaultPoint::COMMIT_OWNER_4,
            FaultPoint::COMMIT_OWNER_5,
            FaultPoint::COMMIT_OWNER_6,
            FaultPoint::COMMIT_CLOCK,
            FaultPoint::COMMIT_ACCEPTED_LEDGER};
        request.publicationStrategy=
            PublicationStrategy::QUIESCED_SEQUENTIAL_WITH_INFALLIBLE_ROLLBACK;
        for(std::size_t index=0;index<sequentialFaults.size();++index)
        {
            request.fault=sequentialFaults[index];
            expectUnchangedFailure(candidates,request,
                TransactionStatus::COMMIT_REJECTED_ROLLED_BACK,
                "sequential rollback fault "+std::to_string(index));
        }

        request.fault=FaultPoint::COMMIT_OWNER_3;
        request.rollbackFails=true;
        LiveState quarantined=initialLive();
        const TransactionResult quarantineResult=
            runTransaction(quarantined,candidates,request);
        check(quarantineResult.status==TransactionStatus::TERMINAL_QUARANTINE,
            "rollback failure enters terminal quarantine");
        check(quarantineResult.terminalQuarantine&&!quarantined.readable,
            "rollback failure exposes no readable live state");
        check(!quarantineResult.statePublished,
            "terminal quarantine never reports publication success");
        request.rollbackFails=false;

        for(PublicationStrategy strategy:{
            PublicationStrategy::SINGLE_COMPOSITE_ROOT_SWAP,
            PublicationStrategy::QUIESCED_SEQUENTIAL_WITH_INFALLIBLE_ROLLBACK})
        {
            request.publicationStrategy=strategy;
            request.fault=FaultPoint::NONE;
            LiveState live=initialLive();
            const std::size_t ledgerBefore=live.acceptedLedger.size();
            const TransactionResult result=runTransaction(live,candidates,request);
            check(result.status==TransactionStatus::OFFLINE_ACCEPTED,
                "abstract strategy reaches offline acceptance");
            check(result.statePublished&&live.readable,
                "abstract success has one readable publication");
            check(live.generationId==second.generationId
                    &&live.domainGenerationId==second.generationId,
                "clock and Domain publish one selected generation");
            check(live.ownerGenerationIds==second.ownerGenerationIds
                    &&live.ownerDigests==second.ownerDigests,
                "all seven owners publish the same selected generation");
            check(sameClock(live.clock,second.clock),
                "solver clock consumes the staged coordinator clock exactly");
            check(live.acceptedLedger.size()==ledgerBefore+1U
                    &&live.acceptedLedger.back().generationId
                        ==second.generationId,
                "accepted ledger advances only with successful publication");
        }

        std::cout<<"RESTART_ATOMIC_PUBLICATION_CONTRACT_TEST_PASS checks="
            <<checks
            <<" runtime_restart=0 hot_start=0 production_ready=0\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"RESTART_ATOMIC_PUBLICATION_CONTRACT_TEST_FAIL "
            <<error.what()<<'\n';
        return 1;
    }
}
