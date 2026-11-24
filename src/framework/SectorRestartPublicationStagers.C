#include "SectorRestartPublicationStagers.H"

#include "SectorTurbulenceOperator.H"

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Foam
{
namespace
{

std::optional<thermo::EquilibriumRegion> equilibriumRegion(
    RestartPhaseAssertion assertion) noexcept
{
    switch(assertion)
    {
        case RestartPhaseAssertion::SINGLE_LIQUID:
            return thermo::EquilibriumRegion::SINGLE_LIQUID;
        case RestartPhaseAssertion::SINGLE_VAPOR:
            return thermo::EquilibriumRegion::SINGLE_VAPOR;
        case RestartPhaseAssertion::SUPERCRITICAL:
            return thermo::EquilibriumRegion::SUPERCRITICAL;
        case RestartPhaseAssertion::TWO_PHASE:
            return thermo::EquilibriumRegion::TWO_PHASE;
        case RestartPhaseAssertion::UNKNOWN:
        case RestartPhaseAssertion::THREE_PHASE_REQUIRED:
            return std::nullopt;
    }
    return std::nullopt;
}

bool sameScalar(double left,double right) noexcept
{
    constexpr double absoluteTolerance=1.0e-12;
    constexpr double relativeTolerance=1.0e-12;
    return std::isfinite(left)&&std::isfinite(right)
        &&std::abs(left-right)<=absoluteTolerance
            +relativeTolerance*std::max(std::abs(left),std::abs(right));
}

bool sameIdentities(
    const RestartIdentitySet& left,
    const RestartIdentitySet& right) noexcept
{
    return left.configurationSha256==right.configurationSha256
        &&left.buildSha256==right.buildSha256
        &&left.meshSha256==right.meshSha256
        &&left.thermoModelId==right.thermoModelId
        &&left.kijManifestSha256==right.kijManifestSha256
        &&left.mechanismSha256==right.mechanismSha256
        &&left.chemistryIntegratorSha256
            ==right.chemistryIntegratorSha256
        &&left.transportSha256==right.transportSha256
        &&left.topologyEnvelopeSha256==right.topologyEnvelopeSha256;
}

bool sameOwnerRecord(
    const RestartOwnerRecord& left,
    const RestartOwnerRecord& right) noexcept
{
    return left.owner==right.owner
        &&left.disposition==right.disposition
        &&left.publicationOrder==right.publicationOrder
        &&left.generationId==right.generationId
        &&left.identitySha256==right.identitySha256
        &&left.payloadBytes==right.payloadBytes
        &&left.payloadSha256==right.payloadSha256;
}

bool matchesExpectation(
    const RestartManifest& manifest,
    const RestartManifestExpectation& expected) noexcept
{
    if(!sameIdentities(manifest.identities,expected.identities)
       ||manifest.byteOrder!=expected.byteOrder
       ||manifest.floatingPointFormat!=expected.floatingPointFormat
       ||manifest.decompositionId!=expected.decompositionId
       ||manifest.cellCount!=expected.cellCount
       ||manifest.generationState.previousAcceptedDtPresent
            !=expected.requirePreviousAcceptedDt
       ||manifest.owners.size()!=restartOwnerCount)
        return false;
    for(std::size_t index=0;index<restartOwnerCount;++index)
    {
        const RestartOwnerRecord& owner=manifest.owners[index];
        const RestartOwnerExpectation& ownerExpected=expected.owners[index];
        if(owner.owner!=ownerExpected.owner
           ||owner.disposition!=ownerExpected.disposition
           ||owner.identitySha256!=ownerExpected.identitySha256)
            return false;
    }
    return true;
}

RestartDomainStepLedger requiredDomainLedger(
    const std::unique_ptr<Domain>& domain)
{
    if(!domain)
        throw std::invalid_argument(
            "detached fluid publication state requires a Domain");
    return {domain->acceptedSteps(),domain->rejectedSteps()};
}

void validateStateRecord(
    const RestartOwnerRecord& record,
    RestartOwner expectedOwner)
{
    if(record.owner!=expectedOwner
       ||record.disposition!=RestartOwnerDisposition::STATE_PAYLOAD
       ||!RestartGenerationStore::validGenerationId(record.generationId)
       ||!validRestartSha256(record.identitySha256)
       ||record.payloadBytes==0U
       ||!validRestartSha256(record.payloadSha256))
        throw std::invalid_argument(
            "detached owner publication record is invalid");
}

bool buildAssertions(
    const RestartManifest& manifest,
    std::vector<SectorRestartFluidRecoveryAssertion>& assertions)
{
    if(manifest.recoverySeeds.size()!=manifest.cellCount) return false;
    assertions.clear();
    assertions.reserve(manifest.recoverySeeds.size());
    for(std::size_t index=0;index<manifest.recoverySeeds.size();++index)
    {
        const RestartRecoverySeed& seed=manifest.recoverySeeds[index];
        const std::optional<thermo::EquilibriumRegion> region=
            equilibriumRegion(seed.expectedPhase);
        if(seed.cell!=index||!region.has_value()) return false;
        SectorRestartFluidRecoveryAssertion assertion;
        assertion.cell=seed.cell;
        assertion.seed.pressure=seed.pressure;
        assertion.seed.temperature=seed.temperature;
        assertion.expectedRegion=*region;
        assertion.phaseIdentitySha256=seed.phaseIdentitySha256;
        assertions.push_back(std::move(assertion));
    }
    return true;
}

RestartOwnerStageResult rejected(
    RestartOwner owner,
    const std::string& detail)
{
    RestartOwnerStageResult result;
    result.owner=owner;
    result.detail=detail;
    return result;
}

RestartPublicationStageResult explicitSelectionMismatch(
    const RestartGenerationLedgerEntry& selected,
    const std::string& explicitGenerationId)
{
    RestartPublicationStageResult result;
    result.status=RestartPublicationStatus::LEDGER_REJECTED;
    result.detail=
        "explicit generation ID does not match the selected ledger entry";
    if(validRestartGenerationLedgerEntry(selected)&&selected.accepted())
    {
        RestartGenerationLedgerEntry rejectedEntry=selected;
        rejectedEntry.disposition=
            RestartGenerationLedgerDisposition::REJECTED;
        rejectedEntry.status=RestartPublicationStatus::LEDGER_REJECTED;
        rejectedEntry.owner.reset();
        rejectedEntry.detail=result.detail+": "+explicitGenerationId;
        if(validRestartGenerationLedgerEntry(rejectedEntry))
            result.rejectedLedgerEntry=std::move(rejectedEntry);
    }
    return result;
}

bool validateFluidAgainstManifest(
    const Domain& domain,
    const RestartManifest& manifest)
{
    if(domain.tripped()||!domain.stateReady()
       ||domain.stepTransactionActive()
       ||static_cast<std::uint64_t>(domain.N())!=manifest.cellCount
       ||domain.thermoModelId()!=manifest.identities.thermoModelId
       ||domain.acceptedSteps()!=manifest.generationState.acceptedIteration
       ||manifest.recoverySeeds.size()!=manifest.cellCount)
        return false;

    for(int cell=0;cell<domain.N();++cell)
    {
        const RestartRecoverySeed& seed=
            manifest.recoverySeeds[static_cast<std::size_t>(cell)];
        const std::optional<thermo::EquilibriumRegion> region=
            equilibriumRegion(seed.expectedPhase);
        const std::optional<std::string> fingerprint=
            SectorRestartOwnerCodec::phaseIdentityForState(
                domain.Q(cell),domain.speciesNames().size());
        if(seed.cell!=static_cast<std::uint64_t>(cell)
           ||!region.has_value()||domain.Q(cell).region!=*region
           ||!sameScalar(domain.Q(cell).p,seed.pressure)
           ||!sameScalar(domain.Q(cell).T,seed.temperature)
           ||!fingerprint.has_value()
           ||*fingerprint!=seed.phaseIdentitySha256)
            return false;
    }
    return true;
}

} // namespace

SectorRestartDetachedFluidState::SectorRestartDetachedFluidState(
    std::unique_ptr<Domain> domain,
    RestartOwnerRecord record)
    : RestartDetachedOwnerState(requiredDomainLedger(domain)),
      domain_(std::move(domain)),record_(std::move(record))
{
    validateStateRecord(record_,RestartOwner::FLUID);
}

SectorRestartDetachedTurbulenceState::
SectorRestartDetachedTurbulenceState(
    std::unique_ptr<SectorRestartTurbulenceState> fields,
    RestartOwnerRecord record)
    : RestartDetachedOwnerState(),fields_(std::move(fields)),
      record_(std::move(record))
{
    if(!fields_)
        throw std::invalid_argument(
            "detached turbulence publication state requires fields");
    validateStateRecord(record_,RestartOwner::TURBULENCE);
}

RestartOwnerStageResult stageSectorRestartFluidOwner(
    const RestartOwnerStageInput& input,
    void* opaque)
{
    if(input.owner.owner!=RestartOwner::FLUID
       ||input.payload.owner!=RestartOwner::FLUID)
        return rejected(input.owner.owner,
            "fluid stager received a non-fluid owner");
    const auto* context=static_cast<const SectorRestartFluidStageContext*>(
        opaque);
    if(!context||!context->targetTemplate)
        return rejected(RestartOwner::FLUID,
            "fluid stager has no detached target template");

    std::vector<SectorRestartFluidRecoveryAssertion> assertions;
    if(!buildAssertions(input.manifest,assertions))
        return rejected(RestartOwner::FLUID,
            "manifest recovery assertions are not canonical");
    SectorRestartFluidDecodeResult decoded=
        SectorRestartOwnerCodec::decodeFluid(
            input.payload.bytes,*context->targetTemplate,assertions);
    if(!decoded.success())
        return rejected(RestartOwner::FLUID,
            std::string("fluid detached decode rejected: ")
                +sectorRestartOwnerCodecStatusName(decoded.status)
                +": "+decoded.detail);
    if(!validateFluidAgainstManifest(*decoded.stagedDomain,input.manifest))
        return rejected(RestartOwner::FLUID,
            "fluid detached state does not match manifest clock/recovery "
            "assertions");

    RestartOwnerStageResult result;
    result.owner=RestartOwner::FLUID;
    result.generationId=input.owner.generationId;
    result.identitySha256=input.owner.identitySha256;
    result.payloadSha256=input.owner.payloadSha256;
    result.state=std::make_shared<const SectorRestartDetachedFluidState>(
        std::move(decoded.stagedDomain),input.owner);
    return result;
}

RestartOwnerStageResult stageSectorRestartTurbulenceOwner(
    const RestartOwnerStageInput& input,
    void* opaque)
{
    if(input.owner.owner!=RestartOwner::TURBULENCE
       ||input.payload.owner!=RestartOwner::TURBULENCE)
        return rejected(input.owner.owner,
            "turbulence stager received a non-turbulence owner");
    const auto* context=
        static_cast<const SectorRestartTurbulenceStageContext*>(opaque);
    if(!context||!context->targetDomain||!context->targetTurbulence)
        return rejected(RestartOwner::TURBULENCE,
            "turbulence stager has no detached target template");

    SectorRestartTurbulenceDecodeResult decoded=
        SectorRestartOwnerCodec::decodeTurbulence(
            input.payload.bytes,*context->targetDomain,
            *context->targetTurbulence);
    if(!decoded.success())
        return rejected(RestartOwner::TURBULENCE,
            std::string("turbulence detached decode rejected: ")
                +sectorRestartOwnerCodecStatusName(decoded.status)
                +": "+decoded.detail);

    RestartOwnerStageResult result;
    result.owner=RestartOwner::TURBULENCE;
    result.generationId=input.owner.generationId;
    result.identitySha256=input.owner.identitySha256;
    result.payloadSha256=input.owner.payloadSha256;
    result.state=
        std::make_shared<const SectorRestartDetachedTurbulenceState>(
            std::move(decoded.stagedState),input.owner);
    return result;
}

RestartCrossValidationResult validateSectorRestartPublicationRoot(
    const RestartPublicationRoot& root,
    void* opaque)
{
    if(!root.valid()) return {false,"publication root is invalid"};
    const auto* context=static_cast<const
        SectorRestartPublicationCrossValidationContext*>(opaque);
    if(context)
    {
        if(!context->explicitGenerationId||!context->expected
           ||!context->targetTemplate||!context->targetTurbulence)
            return {false,"sector cross-validation context is incomplete"};
        if(root.manifest().generationId!=*context->explicitGenerationId
           ||!matchesExpectation(root.manifest(),*context->expected)
           ||context->targetTemplate->thermoModelId()
                !=root.manifest().identities.thermoModelId
           ||static_cast<std::uint64_t>(context->targetTemplate->N())
                !=root.manifest().cellCount
           ||!context->targetTurbulence->meshMatches(
                *context->targetTemplate))
            return {false,
                "explicit generation, expectation, or target identity mismatch"};
    }
    const auto& owners=root.owners();
    const RestartPublishedOwner& fluid=
        owners[static_cast<std::size_t>(RestartOwner::FLUID)];
    const auto* fluidState=dynamic_cast<
        const SectorRestartDetachedFluidState*>(fluid.state.get());
    if(!fluidState
       ||!sameOwnerRecord(fluidState->record(),fluid.record)
       ||!validateFluidAgainstManifest(fluidState->domain(),root.manifest()))
        return {false,"fluid Domain does not match manifest recovery assertions"};

    for(const RestartPublishedOwner& owner:owners)
    {
        if(owner.record.disposition!=RestartOwnerDisposition::STATE_PAYLOAD)
            continue;
        if(owner.record.owner==RestartOwner::FLUID) continue;
        if(owner.record.owner==RestartOwner::TURBULENCE)
        {
            const auto* turbulence=dynamic_cast<
                const SectorRestartDetachedTurbulenceState*>(
                    owner.state.get());
            if(!turbulence
               ||!sameOwnerRecord(turbulence->record(),owner.record)
               ||turbulence->fields().k.size()
                    !=static_cast<std::size_t>(root.manifest().cellCount)
               ||turbulence->fields().omega.size()
                    !=static_cast<std::size_t>(root.manifest().cellCount))
                return {false,"turbulence payload does not match fluid cells"};
            continue;
        }
        return {false,
            "sector publication candidate has no decoder for an enabled owner"};
    }
    return {true,"sector detached owners match one manifest"};
}

SectorRestartExplicitGenerationStageResult
stageSectorRestartExplicitGeneration(
    const std::shared_ptr<const RestartPublicationRoot>& currentRoot,
    const SectorRestartExplicitGenerationStageRequest& request) noexcept
{
    SectorRestartExplicitGenerationStageResult result;
    try
    {
        SectorRestartFluidStageContext fluidContext{
            request.targetTemplate};
        SectorRestartTurbulenceStageContext turbulenceContext{
            request.targetTemplate,request.targetTurbulence};
        SectorRestartPublicationCrossValidationContext crossContext{
            &request.explicitGenerationId,&request.expected,
            request.targetTemplate,request.targetTurbulence};

        RestartPublicationRequest publicationRequest;
        publicationRequest.ledgerEntry=request.selectedLedgerEntry;
        const auto& roster=restartOwnerRoster();
        for(std::size_t index=0;index<roster.size();++index)
            publicationRequest.ownerStagers[index].owner=roster[index];
        const std::size_t fluidIndex=
            static_cast<std::size_t>(RestartOwner::FLUID);
        const std::size_t turbulenceIndex=
            static_cast<std::size_t>(RestartOwner::TURBULENCE);
        if(request.expected.owners[fluidIndex].disposition
            ==RestartOwnerDisposition::STATE_PAYLOAD)
            publicationRequest.ownerStagers[fluidIndex]={
                RestartOwner::FLUID,&stageSectorRestartFluidOwner,
                &fluidContext};
        if(request.expected.owners[turbulenceIndex].disposition
            ==RestartOwnerDisposition::STATE_PAYLOAD)
            publicationRequest.ownerStagers[turbulenceIndex]={
                RestartOwner::TURBULENCE,
                &stageSectorRestartTurbulenceOwner,&turbulenceContext};
        publicationRequest.crossValidate=
            &validateSectorRestartPublicationRoot;
        publicationRequest.crossValidationContext=&crossContext;

        if(request.explicitGenerationId
            !=request.selectedLedgerEntry.generationId)
        {
            result.loadDetail="selection rejected before generation-store load";
            result.staging=explicitSelectionMismatch(
                request.selectedLedgerEntry,request.explicitGenerationId);
            return result;
        }

        RestartGenerationStoreLoadResult selectionProbe;
        selectionProbe.status=RestartGenerationStoreStatus::INVALID_INPUT;
        selectionProbe.detail=
            "selection-only preflight intentionally has no loaded candidate";
        RestartPublicationStageResult selectionAdmission=
            RestartPublicationStager::stage(
                currentRoot,selectionProbe,request.expected,
                request.checkpoint,publicationRequest);
        if(selectionAdmission.status
            !=RestartPublicationStatus::INVALID_STORE_RESULT)
        {
            result.loadDetail="selection rejected before generation-store load";
            result.staging=std::move(selectionAdmission);
            return result;
        }

        result.loadAttempted=true;
        const RestartGenerationStoreLoadResult load=
            RestartGenerationStore::load(
                request.storeRoot,request.explicitGenerationId,
                request.expected);
        result.loadStatus=load.status;
        result.loadDetail=load.detail;
        result.staging=RestartPublicationStager::stage(
            currentRoot,load,request.expected,request.checkpoint,
            publicationRequest);
        return result;
    }
    catch(const std::exception& error)
    {
        result.loadStatus=RestartGenerationStoreStatus::INTERNAL_ERROR;
        try
        {
            result.loadDetail=std::string(
                "unexpected explicit-generation staging exception: ")
                +error.what();
        }
        catch(...)
        {}
        return result;
    }
    catch(...)
    {
        result.loadStatus=RestartGenerationStoreStatus::INTERNAL_ERROR;
        return result;
    }
}

} // End namespace Foam
