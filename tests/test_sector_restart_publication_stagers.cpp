#include "framework/SectorRestartPublicationStagers.H"
#include "framework/RestartLiveRoot.H"
#include "mesh/SectorMeshMetrics.H"
#include "thermo/MixingRules/BinaryInteractionDB.H"
#include "thermo/MixingRules/SpeciesCriticalDB.H"
#include "thermo/PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "thermo/SpeciesRegistry.H"
#include "ThermoDataProvenance.H"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using namespace Foam;
namespace fs=std::filesystem;

int checks=0;

struct QuiescenceContext
{
    bool accepted=true;
    int calls=0;
};

RestartRuntimeQuiescenceResult observeQuiescence(void* opaque)
{
    auto* context=static_cast<QuiescenceContext*>(opaque);
    if(!context) return {false,"missing quiescence context"};
    ++context->calls;
    return {context->accepted,
        context->accepted ? "runtime is quiescent" : "active callback"};
}

RestartCrossValidationResult rejectLiveRoot(
    const RestartPublicationRoot&,
    void*)
{
    return {false,"injected prepublication rejection"};
}

void check(bool condition,const std::string& label)
{
    if(!condition) throw std::runtime_error(label);
    ++checks;
}

std::string hash(const std::string& value)
{
    return thermoProvenance::sha256Hex(value);
}

BinaryInteractionDB completeInteractions(const SpeciesRegistry& registry)
{
    BinaryInteractionDB result=BinaryInteractionDB::createDefault();
    const std::vector<std::string> names=registry.names();
    for(std::size_t left=0;left<names.size();++left)
        for(std::size_t right=left+1;right<names.size();++right)
            if(!result.has(names[left],names[right]))
                result.add(names[left],names[right],0.0);
    return result;
}

phaseEquilibrium::RhoEFlashOptions flashOptions()
{
    phaseEquilibrium::RhoEFlashOptions options;
    options.minimumPressure=1.0e6;
    options.maximumPressure=4.0e6;
    options.minimumTemperature=250.0;
    options.maximumTemperature=1500.0;
    options.initialPressure=2.0e6;
    options.initialTemperature=300.0;
    options.pressureScanIntervals=64;
    options.maximumIterations=120;
    options.finiteDifferenceStep=1.0e-4;
    options.logStepLimit=0.75;
    options.densityAbsoluteTolerance=1.0e-7;
    options.densityRelativeTolerance=2.0e-8;
    options.energyAbsoluteTolerance=1.0e-4;
    options.energyRelativeTolerance=2.0e-10;
    return options;
}

RestartPhaseAssertion restartAssertion(thermo::EquilibriumRegion region)
{
    switch(region)
    {
        case thermo::EquilibriumRegion::SINGLE_LIQUID:
            return RestartPhaseAssertion::SINGLE_LIQUID;
        case thermo::EquilibriumRegion::SINGLE_VAPOR:
            return RestartPhaseAssertion::SINGLE_VAPOR;
        case thermo::EquilibriumRegion::SUPERCRITICAL:
            return RestartPhaseAssertion::SUPERCRITICAL;
        case thermo::EquilibriumRegion::TWO_PHASE:
            return RestartPhaseAssertion::TWO_PHASE;
        case thermo::EquilibriumRegion::UNKNOWN:
        case thermo::EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE:
            break;
    }
    throw std::runtime_error("unknown equilibrium region");
}

RestartCandidatePackage makeCandidate(
    const Domain& source,
    const SectorRestartEncodeResult& fluid,
    const SectorRestartEncodeResult& turbulence,
    const std::string& generationId="sector-generation-0001")
{
    RestartCandidatePackage candidate;
    RestartManifest& manifest=candidate.manifest;
    manifest.generationId=generationId;
    manifest.generationState.physicalTime=2.0e-6;
    manifest.generationState.acceptedIteration=source.acceptedSteps();
    manifest.generationState.previousAcceptedDtPresent=true;
    manifest.generationState.previousAcceptedDt=1.0e-8;
    manifest.byteOrder=nativeRestartByteOrder();
    manifest.cellCount=static_cast<std::uint64_t>(source.N());
    manifest.identities.configurationSha256=hash("configuration");
    manifest.identities.buildSha256=hash("build");
    manifest.identities.meshSha256=hash("mesh");
    manifest.identities.thermoModelId=source.thermoModelId();
    manifest.identities.kijManifestSha256=hash("kij");
    manifest.identities.mechanismSha256=hash("mechanism");
    manifest.identities.chemistryIntegratorSha256=hash("integrator");
    manifest.identities.transportSha256=hash("transport");
    manifest.identities.topologyEnvelopeSha256=hash("topology-envelope");

    const auto& roster=restartOwnerRoster();
    manifest.owners.resize(roster.size());
    for(std::size_t index=0;index<roster.size();++index)
    {
        RestartOwnerRecord& owner=manifest.owners[index];
        owner.owner=roster[index];
        owner.disposition=RestartOwnerDisposition::IDENTITY_ONLY;
        owner.publicationOrder=static_cast<std::uint32_t>(index);
        owner.generationId=manifest.generationId;
        owner.identitySha256=hash(
            std::string("owner-")+restartOwnerName(owner.owner));
        owner.payloadSha256=hash(std::string());
    }
    manifest.owners[static_cast<std::size_t>(RestartOwner::FLUID)]
        .disposition=RestartOwnerDisposition::STATE_PAYLOAD;
    manifest.owners[static_cast<std::size_t>(RestartOwner::TURBULENCE)]
        .disposition=RestartOwnerDisposition::STATE_PAYLOAD;
    manifest.owners[static_cast<std::size_t>(RestartOwner::CHT)]
        .disposition=RestartOwnerDisposition::DISABLED;
    manifest.owners[static_cast<std::size_t>(RestartOwner::CONDENSED_CARBON)]
        .disposition=RestartOwnerDisposition::DISABLED;

    candidate.payloads={{RestartOwner::FLUID,manifest.generationId,
            fluid.bytes,static_cast<std::uint64_t>(fluid.bytes.size())},
        {RestartOwner::TURBULENCE,manifest.generationId,
            turbulence.bytes,
            static_cast<std::uint64_t>(turbulence.bytes.size())}};
    for(const RestartPayloadCandidate& payload:candidate.payloads)
    {
        RestartOwnerRecord& owner=manifest.owners[
            static_cast<std::size_t>(payload.owner)];
        owner.payloadBytes=payload.bytes.size();
        owner.payloadSha256=hash(payload.bytes);
    }

    manifest.recoverySeeds.reserve(static_cast<std::size_t>(source.N()));
    for(int cell=0;cell<source.N();++cell)
    {
        const CompleteCellConversionResult canonical=
            CompleteStateConversion::fromConserved(
                source.U(cell),source.W(cell),source.thermoModel());
        check(canonical.success(),
            "source rho-e-Y has a canonical manifest recovery");
        const thermo::EquilibriumState& equilibrium=canonical.state->Q;
        const std::optional<std::string> phaseIdentity=
            SectorRestartOwnerCodec::phaseIdentityForState(
                equilibrium,source.speciesNames().size());
        check(phaseIdentity.has_value(),
            "source state has a canonical restart phase identity");
        manifest.recoverySeeds.push_back({
            static_cast<std::uint64_t>(cell),equilibrium.p,
            equilibrium.T,restartAssertion(equilibrium.region),
            *phaseIdentity});
    }
    manifest.generationChecksumSha256=restartManifestChecksum(manifest);
    return candidate;
}

RestartManifestExpectation expectation(
    const RestartCandidatePackage& candidate)
{
    RestartManifestExpectation result;
    result.identities=candidate.manifest.identities;
    result.byteOrder=nativeRestartByteOrder();
    result.cellCount=candidate.manifest.cellCount;
    result.requirePreviousAcceptedDt=true;
    const auto& roster=restartOwnerRoster();
    for(std::size_t index=0;index<roster.size();++index)
    {
        result.owners[index].owner=roster[index];
        result.owners[index].disposition=
            candidate.manifest.owners[index].disposition;
        result.owners[index].identitySha256=
            candidate.manifest.owners[index].identitySha256;
    }
    return result;
}

RestartPublicationStageResult stage(
    const RestartCandidatePackage& candidate,
    const RestartManifestExpectation& expected,
    SectorRestartFluidStageContext& fluidContext,
    SectorRestartTurbulenceStageContext& turbulenceContext)
{
    RestartGenerationStoreLoadResult load;
    load.status=RestartGenerationStoreStatus::SUCCESS;
    load.candidate=candidate;

    RestartPublicationRequest request;
    request.ledgerEntry.disposition=
        RestartGenerationLedgerDisposition::ACCEPTED;
    request.ledgerEntry.generationId=candidate.manifest.generationId;
    request.ledgerEntry.ordinal=1U;
    request.ledgerEntry.generationChecksumSha256=
        candidate.manifest.generationChecksumSha256;
    request.ledgerEntry.status=RestartPublicationStatus::SUCCESS;
    const auto& roster=restartOwnerRoster();
    for(std::size_t index=0;index<roster.size();++index)
        request.ownerStagers[index].owner=roster[index];
    request.ownerStagers[static_cast<std::size_t>(RestartOwner::FLUID)]={
        RestartOwner::FLUID,&stageSectorRestartFluidOwner,&fluidContext};
    request.ownerStagers[static_cast<std::size_t>(RestartOwner::TURBULENCE)]={
        RestartOwner::TURBULENCE,&stageSectorRestartTurbulenceOwner,
        &turbulenceContext};
    request.crossValidate=&validateSectorRestartPublicationRoot;

    const RestartCheckpointClockObservation checkpoint{
        candidate.manifest.generationState.physicalTime,
        static_cast<int>(
            candidate.manifest.generationState.acceptedIteration),1};
    return RestartPublicationStager::stage(
        {},load,expected,checkpoint,request);
}

void refreshChecksum(RestartCandidatePackage& candidate)
{
    candidate.manifest.generationChecksumSha256=
        restartManifestChecksum(candidate.manifest);
}

class ScopedDirectory
{
    fs::path path_;

public:
    explicit ScopedDirectory(const std::string& label)
    {
        const auto tick=std::chrono::high_resolution_clock::now()
            .time_since_epoch().count();
        path_=fs::temp_directory_path()/(
            "rfam-sector-restart-publication-"+label+"-"
            +std::to_string(tick));
        std::error_code error;
        fs::remove_all(path_,error);
        error.clear();
        if(!fs::create_directory(path_,error)||error)
            throw std::runtime_error(
                "cannot create sector restart publication store");
    }

    ~ScopedDirectory()
    {
        std::error_code error;
        fs::remove_all(path_,error);
    }

    const fs::path& path() const noexcept { return path_; }
};

std::string readBytes(const fs::path& path)
{
    std::ifstream stream(path,std::ios::binary);
    if(!stream) throw std::runtime_error("cannot read store payload");
    return std::string(
        std::istreambuf_iterator<char>(stream),
        std::istreambuf_iterator<char>());
}

void writeBytes(const fs::path& path,const std::string& bytes)
{
    std::ofstream stream(path,std::ios::binary|std::ios::trunc);
    if(!stream) throw std::runtime_error("cannot write store payload");
    stream.write(bytes.data(),static_cast<std::streamsize>(bytes.size()));
    if(!stream) throw std::runtime_error("store payload write failed");
}

SectorRestartExplicitGenerationStageRequest explicitRequest(
    const fs::path& storeRoot,
    const RestartCandidatePackage& candidate,
    const RestartManifestExpectation& expected,
    const Domain& target,
    const SectorTurbulenceOperator& targetTurbulence,
    std::uint64_t ordinal=1U)
{
    SectorRestartExplicitGenerationStageRequest request;
    request.storeRoot=storeRoot;
    request.explicitGenerationId=candidate.manifest.generationId;
    request.selectedLedgerEntry.disposition=
        RestartGenerationLedgerDisposition::ACCEPTED;
    request.selectedLedgerEntry.generationId=
        candidate.manifest.generationId;
    request.selectedLedgerEntry.ordinal=ordinal;
    request.selectedLedgerEntry.generationChecksumSha256=
        candidate.manifest.generationChecksumSha256;
    request.selectedLedgerEntry.status=RestartPublicationStatus::SUCCESS;
    request.expected=expected;
    request.checkpoint={
        candidate.manifest.generationState.physicalTime,
        static_cast<int>(
            candidate.manifest.generationState.acceptedIteration),1};
    request.targetTemplate=&target;
    request.targetTurbulence=&targetTurbulence;
    return request;
}

} // namespace

int main()
{
    try
    {
        SpeciesRegistry registry;
        const SpeciesCriticalDB critical=SpeciesCriticalDB::createDefault();
        const BinaryInteractionDB interactions=completeInteractions(registry);
        const auto model=std::make_shared<const
            thermo::PREquilibriumThermoModel>(
                critical,interactions,registry.names(),flashOptions());

        std::vector<double> composition(
            static_cast<std::size_t>(registry.size()),0.0);
        phaseEquilibrium::PRTPFlash feedFlash(
            critical,interactions,{"O2","C3H8"});
        const phaseEquilibrium::TPFlashResult feed=
            feedFlash.flashMoleFractions(2.0e6,300.0,{0.2,0.8});
        check(feed.valid()
                &&feed.region==phaseEquilibrium::EquilibriumRegion::TWO_PHASE,
            "fixture feed is two phase");
        constexpr double trace=1.0e-4;
        const double reactantScale=1.0-2.0*trace;
        composition[registry.index("O2")]=
            reactantScale*feed.feedMassFractions[0];
        composition[registry.index("C3H8")]=
            reactantScale*feed.feedMassFractions[1];
        composition[registry.index("CO")]=trace;
        composition[registry.index("H2O")]=trace;
        const thermo::EquilibriumResult initial=
            model->initializePTY(2.0e6,300.0,composition);
        check(initial.success()
                &&initial.state->region==thermo::EquilibriumRegion::TWO_PHASE,
            "fixture model supplies a complete HEM state");
        const std::optional<std::string> phaseCaloricIdentity=
            SectorRestartOwnerCodec::phaseIdentityForState(
                *initial.state,registry.names().size());
        thermo::EquilibriumState changedCaloric=*initial.state;
        changedCaloric.h=std::nextafter(
            changedCaloric.h,std::numeric_limits<double>::infinity());
        const std::optional<std::string> changedCaloricIdentity=
            SectorRestartOwnerCodec::phaseIdentityForState(
                changedCaloric,registry.names().size());
        thermo::EquilibriumState changedPhaseCaloric=*initial.state;
        check(changedPhaseCaloric.liquid.has_value(),
            "two-phase fixture has a liquid caloric payload");
        changedPhaseCaloric.liquid->h=std::nextafter(
            changedPhaseCaloric.liquid->h,
            std::numeric_limits<double>::infinity());
        const std::optional<std::string> changedPhaseCaloricIdentity=
            SectorRestartOwnerCodec::phaseIdentityForState(
                changedPhaseCaloric,registry.names().size());
        check(phaseCaloricIdentity.has_value()
                &&changedCaloricIdentity.has_value()
                &&changedPhaseCaloricIdentity.has_value()
                &&*phaseCaloricIdentity!=*changedCaloricIdentity
                &&*phaseCaloricIdentity!=*changedPhaseCaloricIdentity,
            "restart fingerprint binds mixture and phase caloric enthalpy");

        constexpr int extent=4;
        SectorMeshMetrics mesh(
            vortexProfilePath(),extent,extent,extent,120.0);
        Domain source(
            extent,extent,extent,mesh.length(),0.065,mesh.sectorAngle(),
            registry.size(),model);
        source.setSpecies(registry.names());
        PrimitiveVars primitive(registry.size());
        primitive.p=initial.state->p;
        primitive.T=initial.state->T;
        primitive.U={0.5,-0.25,0.125};
        primitive.Y=composition;
        check(source.replaceCompleteField(
                std::vector<PrimitiveVars>(source.N(),primitive),
                std::vector<thermo::EquilibriumState>(
                    source.N(),*initial.state)),
            "source Domain publishes authoritative rho-e-Y/Q");
        check(source.beginStepTransaction(),"accepted step starts");
        source.acceptStepTransaction();
        check(source.beginStepTransaction(),"rejected step starts");
        source.rollbackStepTransaction();

        SectorTurbulenceOperator sourceTurbulence(mesh,registry);
        std::vector<double> sourceK=sourceTurbulence.kField();
        std::vector<double> sourceOmega=sourceTurbulence.omegaField();
        for(double& value:sourceK) value*=1.5;
        for(double& value:sourceOmega) value*=1.25;
        sourceTurbulence.restoreFields(sourceK,sourceOmega);
        const SectorRestartEncodeResult fluid=
            SectorRestartOwnerCodec::encodeFluid(source);
        const SectorRestartEncodeResult turbulence=
            SectorRestartOwnerCodec::encodeTurbulence(
                source,sourceTurbulence);
        check(fluid.success()&&turbulence.success(),
            "fluid and turbulence payloads encode");

        Domain target(
            extent,extent,extent,mesh.length(),0.065,mesh.sectorAngle(),
            registry.size(),model);
        target.setSpecies(registry.names());
        SectorTurbulenceOperator targetTurbulence(mesh,registry);
        const std::vector<double> targetK=targetTurbulence.kField();
        const std::vector<double> targetOmega=targetTurbulence.omegaField();
        SectorRestartFluidStageContext fluidContext{&target};
        SectorRestartTurbulenceStageContext turbulenceContext{
            &target,&targetTurbulence};

        const RestartCandidatePackage canonical=
            makeCandidate(source,fluid,turbulence);
        const RestartManifestExpectation expected=expectation(canonical);
        ScopedDirectory store("actual-load");
        const RestartGenerationStorePublishResult stored=
            RestartGenerationStore::publish(
                store.path(),canonical,expected);
        check(stored.success(),
            "canonical generation commits to the actual store");
        const SectorRestartExplicitGenerationStageRequest canonicalRequest=
            explicitRequest(
                store.path(),canonical,expected,target,targetTurbulence);
        SectorRestartExplicitGenerationStageResult loadedStage=
            stageSectorRestartExplicitGeneration({},canonicalRequest);
        RestartPublicationStageResult& staged=loadedStage.staging;
        if(!loadedStage.success())
            std::cerr<<"canonical stage status="
                     <<restartPublicationStatusName(staged.status)
                     <<" detail="<<staged.detail
                     <<" load="
                     <<restartGenerationStoreStatusName(
                        loadedStage.loadStatus)
                     <<" load-detail="<<loadedStage.loadDetail<<'\n';
        check(loadedStage.success()
                &&staged.publication->candidate().valid(),
            "actual generation load stages real owners in one immutable root");
        const auto* stagedFluid=dynamic_cast<
            const SectorRestartDetachedFluidState*>(
                staged.publication->candidate().owners()[
                    static_cast<std::size_t>(RestartOwner::FLUID)]
                    .state.get());
        const auto* stagedTurbulence=dynamic_cast<
            const SectorRestartDetachedTurbulenceState*>(
                staged.publication->candidate().owners()[
                    static_cast<std::size_t>(RestartOwner::TURBULENCE)]
                    .state.get());
        check(stagedFluid&&stagedTurbulence
                &&stagedFluid->domain().N()==source.N()
                &&stagedFluid->domain().acceptedSteps()
                    ==source.acceptedSteps()
                &&stagedFluid->domain().rejectedSteps()
                    ==source.rejectedSteps()
                &&stagedTurbulence->fields().k==sourceK
                &&stagedTurbulence->fields().omega==sourceOmega
                &&stagedFluid->record().identitySha256
                    ==canonical.manifest.owners[
                        static_cast<std::size_t>(RestartOwner::FLUID)]
                        .identitySha256
                &&stagedTurbulence->record().payloadSha256
                    ==canonical.manifest.owners[
                        static_cast<std::size_t>(RestartOwner::TURBULENCE)]
                        .payloadSha256,
            "loaded candidate contains exact detached owners and bindings");
        check(!target.stateReady()
                &&targetTurbulence.kField()==targetK
                &&targetTurbulence.omegaField()==targetOmega,
            "staging does not mutate live templates");

        RestartLiveRoot routedLive;
        QuiescenceContext quiescence;
        RestartLiveRootPublicationRequest routedRequest;
        routedRequest.observeQuiescence=&observeQuiescence;
        routedRequest.quiescenceContext=&quiescence;
        routedRequest.prepublishValidate=
            &validateSectorRestartPublicationRoot;

        SectorRestartExplicitGenerationStageResult routedStage=
            stageSectorRestartExplicitGeneration(
                routedLive.stagingBase(),canonicalRequest);
        check(routedStage.success(),
            "live-root fixture stages from the routed empty base");
        {
            RestartLiveRootReadLease emptyRead=routedLive.read();
            check(!emptyRead.available(),
                "empty live-root read lease exposes no partial candidate");
            const RestartLiveRootPublicationResult readersActive=
                routedLive.tryPublish(
                    *routedStage.staging.publication,routedRequest);
            check(readersActive.status
                    ==RestartLiveRootStatus::READERS_OR_WRITER_ACTIVE,
                "active routed reader rejects publication without consuming it");
        }
        check(!routedLive.stagingBase(),
            "reader-conflicted publication leaves the routed root absent");

        quiescence.accepted=false;
        const RestartLiveRootPublicationResult notQuiescent=
            routedLive.tryPublish(
                *routedStage.staging.publication,routedRequest);
        check(notQuiescent.status
                ==RestartLiveRootStatus::QUIESCENCE_REJECTED
                &&!routedLive.stagingBase(),
            "external non-quiescence leaves the routed root absent");
        quiescence.accepted=true;

        RestartLiveRootPublicationRequest rejectedRequest=routedRequest;
        rejectedRequest.prepublishValidate=&rejectLiveRoot;
        const RestartLiveRootPublicationResult prepublishRejected=
            routedLive.tryPublish(
                *routedStage.staging.publication,rejectedRequest);
        check(prepublishRejected.status
                ==RestartLiveRootStatus::PREPUBLICATION_REJECTED
                &&!routedLive.stagingBase(),
            "prepublication phase/caloric rejection leaves no routed root");

        RestartLiveRootPublicationRequest rollbackRequest=routedRequest;
        rollbackRequest.fault=
            RestartLiveRootFaultPoint::AFTER_ROOT_SWAP;
        const RestartLiveRootPublicationResult rolledBack=
            routedLive.tryPublish(
                *routedStage.staging.publication,rollbackRequest);
        check(rolledBack.status
                ==RestartLiveRootStatus::TEST_FAULT_ROLLED_BACK
                &&rolledBack.rollbackCompleted
                &&!routedLive.stagingBase(),
            "post-swap fault rolls the root back before readers can enter");

        SectorRestartExplicitGenerationStageResult routedSuccessStage=
            stageSectorRestartExplicitGeneration(
                routedLive.stagingBase(),canonicalRequest);
        SectorRestartExplicitGenerationStageResult routedStaleStage=
            stageSectorRestartExplicitGeneration(
                routedLive.stagingBase(),canonicalRequest);
        check(routedSuccessStage.success()&&routedStaleStage.success(),
            "success and stale probes stage against the same immutable base");
        const RestartLiveRootPublicationResult routedPublished=
            routedLive.tryPublish(
                *routedSuccessStage.staging.publication,routedRequest);
        check(routedPublished.success(),
            "quiescent publication installs one complete routed root");
        {
            RestartLiveRootReadLease liveRead=routedLive.read();
            const RestartPublicationRoot* visible=liveRead.get();
            const auto* visibleFluid=visible ? dynamic_cast<
                const SectorRestartDetachedFluidState*>(
                    visible->owners()[static_cast<std::size_t>(
                        RestartOwner::FLUID)].state.get()) : nullptr;
            check(liveRead.available()&&visible&&visibleFluid
                    &&visible->manifest().generationId
                        ==canonical.manifest.generationId
                    &&visible->clock().physicalTime
                        ==canonical.manifest.generationState.physicalTime
                    &&visible->clock().acceptedIteration
                        ==static_cast<int>(canonical.manifest.generationState
                            .acceptedIteration)
                    &&visible->clock().previousAcceptedDt
                        ==std::optional<double>(canonical.manifest
                            .generationState.previousAcceptedDt)
                    &&visibleFluid->domain().acceptedSteps()
                        ==canonical.manifest.generationState.acceptedIteration,
                "one read lease observes matching generation, clock, dt, and Domain");
        }
        const std::shared_ptr<const RestartPublicationRoot> routedWinner=
            routedLive.stagingBase();
        const RestartLiveRootPublicationResult staleRouted=
            routedLive.tryPublish(
                *routedStaleStage.staging.publication,routedRequest);
        check(staleRouted.status
                ==RestartLiveRootStatus::STAGED_PUBLICATION_REJECTED
                &&routedLive.stagingBase().get()==routedWinner.get(),
            "stale staged base cannot replace the routed winner");

        check(source.beginStepTransaction(),
            "successor fixture starts one accepted step");
        source.acceptStepTransaction();
        const SectorRestartEncodeResult successorFluid=
            SectorRestartOwnerCodec::encodeFluid(source);
        const SectorRestartEncodeResult successorTurbulence=
            SectorRestartOwnerCodec::encodeTurbulence(
                source,sourceTurbulence);
        check(successorFluid.success()&&successorTurbulence.success(),
            "successor fluid and turbulence payloads encode");
        const RestartCandidatePackage successor=makeCandidate(
            source,successorFluid,successorTurbulence,
            "sector-generation-0002");
        const RestartManifestExpectation successorExpected=
            expectation(successor);
        check(RestartGenerationStore::publish(
                store.path(),successor,successorExpected).success(),
            "successor generation commits to the actual store");
        SectorRestartExplicitGenerationStageRequest successorRequest=
            explicitRequest(
                store.path(),successor,successorExpected,
                target,targetTurbulence);
        successorRequest.selectedLedgerEntry.ordinal=2U;
        SectorRestartExplicitGenerationStageResult successorStage=
            stageSectorRestartExplicitGeneration(
                routedLive.stagingBase(),successorRequest);
        check(successorStage.success(),
            "successor stages against the current routed winner");
        const RestartLiveRootPublicationResult successorRollback=
            routedLive.tryPublish(
                *successorStage.staging.publication,rollbackRequest);
        const std::shared_ptr<const RestartPublicationRoot> afterRollback=
            routedLive.stagingBase();
        check(successorRollback.status
                ==RestartLiveRootStatus::TEST_FAULT_ROLLED_BACK
                &&successorRollback.rollbackCompleted
                &&afterRollback.get()==routedWinner.get()
                &&afterRollback->manifest().generationId
                    ==canonical.manifest.generationId
                &&afterRollback->clock().acceptedIteration
                    ==static_cast<int>(canonical.manifest.generationState
                        .acceptedIteration),
            "successor fault restores the exact nonempty root and clock");

        std::shared_ptr<const RestartPublicationRoot> live;
        SectorRestartExplicitGenerationStageRequest missingRequest=
            canonicalRequest;
        missingRequest.explicitGenerationId="missing-generation";
        missingRequest.selectedLedgerEntry.generationId=
            missingRequest.explicitGenerationId;
        missingRequest.selectedLedgerEntry.generationChecksumSha256=
            hash("missing-generation-manifest");
        const SectorRestartExplicitGenerationStageResult missing=
            stageSectorRestartExplicitGeneration(live,missingRequest);
        check(missing.loadAttempted
                &&missing.loadStatus
                    ==RestartGenerationStoreStatus::NOT_COMMITTED
                &&missing.staging.status
                    ==RestartPublicationStatus::INVALID_STORE_RESULT
                &&!missing.staging.publication&&!live,
            "missing explicit generation exposes no live publication");

        SectorRestartExplicitGenerationStageRequest identityRequest=
            canonicalRequest;
        identityRequest.expected.identities.buildSha256=
            hash("wrong-live-build-identity");
        const SectorRestartExplicitGenerationStageResult identityRejected=
            stageSectorRestartExplicitGeneration(live,identityRequest);
        check(identityRejected.loadAttempted
                &&identityRejected.loadStatus
                ==RestartGenerationStoreStatus::MANIFEST_REJECTED
                &&identityRejected.staging.status
                    ==RestartPublicationStatus::INVALID_STORE_RESULT
                &&!identityRejected.staging.publication&&!live,
            "independent owner/manifest expectation mismatch fails at load");

        SectorRestartExplicitGenerationStageRequest ledgerRequest=
            canonicalRequest;
        ledgerRequest.selectedLedgerEntry.generationChecksumSha256=
            hash("wrong-selected-ledger-checksum");
        const SectorRestartExplicitGenerationStageResult ledgerBindingRejected=
            stageSectorRestartExplicitGeneration(live,ledgerRequest);
        check(ledgerBindingRejected.loadAttempted
                &&ledgerBindingRejected.loadStatus
                ==RestartGenerationStoreStatus::SUCCESS
                &&ledgerBindingRejected.staging.status
                    ==RestartPublicationStatus::LEDGER_REJECTED
                &&!ledgerBindingRejected.staging.publication&&!live,
            "loaded generation must match the selected ledger binding");

        SectorRestartExplicitGenerationStageRequest clockRequest=
            canonicalRequest;
        clockRequest.checkpoint.physicalTime+=1.0e-12;
        const SectorRestartExplicitGenerationStageResult clockRejected=
            stageSectorRestartExplicitGeneration(live,clockRequest);
        check(clockRejected.loadAttempted
                &&clockRejected.loadStatus
                    ==RestartGenerationStoreStatus::SUCCESS
                &&clockRejected.staging.status
                    ==RestartPublicationStatus::CLOCK_REJECTED
                &&!clockRejected.staging.publication&&!live,
            "store/manifest/checkpoint clock mismatch exposes no root");

        const fs::path fluidPath=
            RestartGenerationStore::generationDirectory(
                store.path(),canonical.manifest.generationId)
            /RestartGenerationStore::ownerPayloadFileName(
                RestartOwner::FLUID);
        const std::string storedFluid=readBytes(fluidPath);
        std::string corruptedFluid=storedFluid;
        corruptedFluid[corruptedFluid.size()/2U]^=1;
        writeBytes(fluidPath,corruptedFluid);
        const SectorRestartExplicitGenerationStageResult payloadRejected=
            stageSectorRestartExplicitGeneration(live,canonicalRequest);
        check(payloadRejected.loadAttempted
                &&payloadRejected.loadStatus
                ==RestartGenerationStoreStatus::OWNER_CHECKSUM_MISMATCH
                &&payloadRejected.staging.status
                    ==RestartPublicationStatus::INVALID_STORE_RESULT
                &&!payloadRejected.staging.publication&&!live,
            "tampered on-disk fluid payload fails before owner staging");
        writeBytes(fluidPath,storedFluid);

        Domain wrongTarget(
            extent,extent,extent,mesh.length(),0.065,mesh.sectorAngle(),
            registry.size());
        wrongTarget.setSpecies(registry.names());
        SectorTurbulenceOperator wrongTargetTurbulence(mesh,registry);
        SectorRestartExplicitGenerationStageRequest targetRequest=
            canonicalRequest;
        targetRequest.targetTemplate=&wrongTarget;
        targetRequest.targetTurbulence=&wrongTargetTurbulence;
        const SectorRestartExplicitGenerationStageResult targetRejected=
            stageSectorRestartExplicitGeneration(live,targetRequest);
        check(targetRejected.loadAttempted
                &&targetRejected.loadStatus
                    ==RestartGenerationStoreStatus::SUCCESS
                &&targetRejected.staging.status
                    ==RestartPublicationStatus::OWNER_STAGE_REJECTED
                &&!targetRejected.staging.publication&&!live,
            "target thermo identity mismatch fails detached fluid staging");

        std::string wrongGeneration="wrong-cross-generation";
        SectorRestartPublicationCrossValidationContext wrongCrossContext{
            &wrongGeneration,&canonicalRequest.expected,
            &target,&targetTurbulence};
        check(!validateSectorRestartPublicationRoot(
                staged.publication->candidate(),&wrongCrossContext).accepted,
            "cross validator independently binds the explicit generation");

        check(!live,
            "all actual-load failures leave the live root uncreated");
        check(staged.publication->publish(live).success()
                &&live&&live->valid(),
            "one root swap publishes the complete detached candidate");
        const std::shared_ptr<const RestartPublicationRoot> publishedLive=live;

        SectorRestartExplicitGenerationStageRequest mismatchedSelection=
            canonicalRequest;
        mismatchedSelection.storeRoot=fs::path(
            "this-store-path-must-not-be-read");
        mismatchedSelection.explicitGenerationId="different-generation";
        mismatchedSelection.selectedLedgerEntry.ordinal=2U;
        const SectorRestartExplicitGenerationStageResult selectionRejected=
            stageSectorRestartExplicitGeneration(
                live,mismatchedSelection);
        check(!selectionRejected.loadAttempted
                &&selectionRejected.staging.status
                    ==RestartPublicationStatus::LEDGER_REJECTED
                &&!selectionRejected.staging.publication
                &&live.get()==publishedLive.get(),
            "generation/ledger selection mismatch rejects before store I/O");

        SectorRestartExplicitGenerationStageRequest duplicateSelection=
            canonicalRequest;
        duplicateSelection.selectedLedgerEntry.ordinal=2U;
        const SectorRestartExplicitGenerationStageResult duplicateRejected=
            stageSectorRestartExplicitGeneration(
                live,duplicateSelection);
        check(!duplicateRejected.loadAttempted
                &&duplicateRejected.staging.status
                    ==RestartPublicationStatus::LEDGER_CONFLICT
                &&!duplicateRejected.staging.publication
                &&live.get()==publishedLive.get(),
            "live-ledger generation conflict rejects before store I/O");

        SectorRestartExplicitGenerationStageRequest postPublishMissing=
            missingRequest;
        postPublishMissing.selectedLedgerEntry.ordinal=2U;
        const SectorRestartExplicitGenerationStageResult postFailure=
            stageSectorRestartExplicitGeneration(
                live,postPublishMissing);
        check(postFailure.loadAttempted
                &&postFailure.loadStatus
                ==RestartGenerationStoreStatus::NOT_COMMITTED
                &&!postFailure.staging.publication
                &&live.get()==publishedLive.get(),
            "failed successor load leaves the published root pointer exact");

        RestartCandidatePackage wrongPressure=canonical;
        wrongPressure.manifest.recoverySeeds[0].pressure*=1.0001;
        refreshChecksum(wrongPressure);
        RestartPublicationStageResult pressureRejected=stage(
            wrongPressure,expectation(wrongPressure),
            fluidContext,turbulenceContext);
        check(pressureRejected.status
                ==RestartPublicationStatus::OWNER_STAGE_REJECTED
                &&!pressureRejected.publication,
            "manifest pressure seed mismatch fails in the fluid stager");

        RestartCandidatePackage wrongTemperature=canonical;
        wrongTemperature.manifest.recoverySeeds[0].temperature+=1.0;
        refreshChecksum(wrongTemperature);
        RestartPublicationStageResult temperatureRejected=stage(
            wrongTemperature,expectation(wrongTemperature),
            fluidContext,turbulenceContext);
        check(temperatureRejected.status
                ==RestartPublicationStatus::OWNER_STAGE_REJECTED,
            "manifest temperature seed mismatch fails in the fluid stager");

        RestartCandidatePackage wrongRegion=canonical;
        wrongRegion.manifest.recoverySeeds[0].expectedPhase=
            RestartPhaseAssertion::SINGLE_VAPOR;
        refreshChecksum(wrongRegion);
        RestartPublicationStageResult regionRejected=stage(
            wrongRegion,expectation(wrongRegion),
            fluidContext,turbulenceContext);
        check(regionRejected.status
                ==RestartPublicationStatus::OWNER_STAGE_REJECTED,
            "manifest phase-region mismatch fails in the fluid stager");

        RestartCandidatePackage wrongFingerprint=canonical;
        wrongFingerprint.manifest.recoverySeeds[0].phaseIdentitySha256=
            hash("wrong-phase-identity");
        refreshChecksum(wrongFingerprint);
        RestartPublicationStageResult fingerprintRejected=stage(
            wrongFingerprint,expectation(wrongFingerprint),
            fluidContext,turbulenceContext);
        check(fingerprintRejected.status
                ==RestartPublicationStatus::OWNER_STAGE_REJECTED,
            "manifest phase fingerprint mismatch fails in the fluid stager");

        RestartCandidatePackage wrongLedger=canonical;
        ++wrongLedger.manifest.generationState.acceptedIteration;
        refreshChecksum(wrongLedger);
        RestartPublicationStageResult ledgerRejected=stage(
            wrongLedger,expectation(wrongLedger),
            fluidContext,turbulenceContext);
        if(ledgerRejected.status
            !=RestartPublicationStatus::OWNER_STAGE_REJECTED)
            std::cerr<<"ledger stage status="
                     <<restartPublicationStatusName(ledgerRejected.status)
                     <<" detail="<<ledgerRejected.detail<<'\n';
        check(ledgerRejected.status
                ==RestartPublicationStatus::OWNER_STAGE_REJECTED,
            "decoded Domain ledger must match the coordinator clock");
        check(!target.stateReady()
                &&targetTurbulence.kField()==targetK
                &&targetTurbulence.omegaField()==targetOmega,
            "all rejected candidates leave live templates unchanged");

        std::cout
            <<"SECTOR_RESTART_PUBLICATION_STAGERS_PASS checks="<<checks
            <<" scope=offline-two-owner-hem-live-root production_ready=0\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"SECTOR_RESTART_PUBLICATION_STAGERS_FAIL "
                 <<error.what()<<'\n';
        return 1;
    }
}
