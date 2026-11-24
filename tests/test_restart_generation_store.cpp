#include "RestartGenerationStore.H"
#include "ThermoDataProvenance.H"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using namespace Foam;

namespace
{

namespace fs = std::filesystem;

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

RestartCandidatePackage validCandidate(
    const std::string& generationId,
    const std::string& payloadTag)
{
    RestartCandidatePackage candidate;
    RestartManifest& manifest=candidate.manifest;
    manifest.generationId=generationId;
    manifest.generationState.physicalTime=1.25e-5;
    manifest.generationState.acceptedIteration=42U;
    manifest.generationState.previousAcceptedDtPresent=true;
    manifest.generationState.previousAcceptedDt=1.0e-8;
    manifest.byteOrder=nativeRestartByteOrder();
    manifest.cellCount=2U;
    manifest.identities.configurationSha256=hash("configuration");
    manifest.identities.buildSha256=hash("build");
    manifest.identities.meshSha256=hash("mesh");
    manifest.identities.thermoModelId=
        "peng-robinson-rho-e:sha256="+hash("thermo");
    manifest.identities.kijManifestSha256=hash("kij");
    manifest.identities.mechanismSha256=hash("mechanism");
    manifest.identities.chemistryIntegratorSha256=
        hash("chemistry-integrator");
    manifest.identities.transportSha256=hash("transport");
    manifest.identities.topologyEnvelopeSha256=hash("topology-envelope");

    const auto& roster=restartOwnerRoster();
    manifest.owners.resize(roster.size());
    for(std::size_t i=0;i<roster.size();++i)
    {
        RestartOwnerRecord& owner=manifest.owners[i];
        owner.owner=roster[i];
        owner.publicationOrder=static_cast<std::uint32_t>(i);
        owner.generationId=generationId;
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
        {RestartOwner::FLUID,generationId,
            std::string("fluid\0state\xff",12)+payloadTag,0U},
        {RestartOwner::TURBULENCE,generationId,
            "turbulence-state-v2-"+payloadTag,0U}
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
        {0U,5.0e6,155.0,RestartPhaseAssertion::SINGLE_LIQUID,
            hash("phase-cell-0")},
        {1U,5.0e6,180.0,RestartPhaseAssertion::TWO_PHASE,
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

class ScopedDirectory
{
    fs::path path_;

public:
    explicit ScopedDirectory(const std::string& label)
    {
        static std::uint64_t sequence=0U;
        const auto tick=std::chrono::high_resolution_clock::now()
            .time_since_epoch().count();
        path_=fs::temp_directory_path()/(
            "rfam-restart-store-"+label+"-"+std::to_string(tick)
            +"-"+std::to_string(++sequence));
        std::error_code error;
        fs::remove_all(path_,error);
        error.clear();
        if(!fs::create_directory(path_,error)||error)
            throw std::runtime_error("cannot create test directory");
    }

    ~ScopedDirectory()
    {
        std::error_code error;
        fs::remove_all(path_,error);
    }

    const fs::path& path() const
    {
        return path_;
    }
};

void writeBytes(const fs::path& path,const std::string& bytes)
{
    std::ofstream output(path,std::ios::binary|std::ios::trunc);
    if(!output.is_open()) throw std::runtime_error("test write open failed");
    output.write(bytes.data(),static_cast<std::streamsize>(bytes.size()));
    output.flush();
    output.close();
    if(output.fail()) throw std::runtime_error("test write failed");
}

std::string readBytes(const fs::path& path)
{
    std::ifstream input(path,std::ios::binary);
    if(!input.is_open()) throw std::runtime_error("test read open failed");
    return std::string(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>());
}

std::map<std::string,std::string> directorySnapshot(const fs::path& path)
{
    std::map<std::string,std::string> snapshot;
    for(const fs::directory_entry& entry:fs::directory_iterator(path))
    {
        if(entry.is_regular_file())
            snapshot[entry.path().filename().string()]=readBytes(entry.path());
    }
    return snapshot;
}

fs::path generationPath(
    const fs::path& root,
    const RestartCandidatePackage& candidate)
{
    return RestartGenerationStore::generationDirectory(
        root,candidate.manifest.generationId);
}

fs::path manifestPath(
    const fs::path& root,
    const RestartCandidatePackage& candidate)
{
    return generationPath(root,candidate)
        /RestartGenerationStore::manifestFileName;
}

fs::path payloadPath(
    const fs::path& root,
    const RestartCandidatePackage& candidate,
    RestartOwner owner)
{
    return generationPath(root,candidate)
        /RestartGenerationStore::ownerPayloadFileName(owner);
}

void checkNoCandidate(
    const RestartGenerationStoreLoadResult& result,
    const std::string& label)
{
    check(!result.success()&&!result.candidate.has_value(),label);
}

void checkFault(
    RestartGenerationStoreFaultPoint point,
    RestartOwner owner,
    const std::string& label,
    RestartGenerationStoreStatus expectedStatus)
{
    ScopedDirectory root("fault-"+label);
    const RestartCandidatePackage candidate=validCandidate(
        "generation-fault-"+label,"fault-"+label);
    const RestartManifestExpectation expected=expectation(candidate);
    RestartGenerationStoreFault fault;
    fault.point=point;
    fault.owner=owner;
    const RestartGenerationStorePublishResult publish=
        RestartGenerationStore::publish(
            root.path(),candidate,expected,fault);
    check(!publish.success()&&publish.status==expectedStatus,
        label+" returns the exact structured failure");
    check(!fs::exists(manifestPath(root.path(),candidate)),
        label+" leaves no final manifest marker");
    checkNoCandidate(RestartGenerationStore::load(
        root.path(),candidate.manifest.generationId,expected),
        label+" exposes no candidate");

    const fs::path generation=generationPath(root.path(),candidate);
    bool temporaryReopenable=true;
    if(fs::exists(generation))
    {
        for(const fs::directory_entry& entry:fs::directory_iterator(generation))
        {
            if(entry.path().extension()==".tmp")
            {
                std::ofstream append(entry.path(),std::ios::binary|std::ios::app);
                temporaryReopenable=temporaryReopenable&&append.is_open();
                append.close();
            }
        }
    }
    check(temporaryReopenable,label+" leaves no open staged-file handle");
}

struct RootSwapContext
{
    fs::path admittedPath;
    fs::path movedPath;
    std::error_code error;
    bool swapped = false;
};

void replaceRootIdentity(void* opaque) noexcept
{
    RootSwapContext* context=static_cast<RootSwapContext*>(opaque);
    if(context==nullptr) return;
    fs::rename(context->admittedPath,context->movedPath,context->error);
    if(context->error) return;
    if(!fs::create_directory(context->admittedPath,context->error)
       ||context->error)
        return;
    context->swapped=true;
}

} // End anonymous namespace

int main()
{
    try
    {
        const RestartCandidatePackage baseline=validCandidate(
            "generation-store-20260727-0001","baseline");
        const RestartManifestExpectation expected=expectation(baseline);

        check(std::string(restartGenerationStoreStatusName(
            RestartGenerationStoreStatus::SUCCESS))=="SUCCESS",
            "store status names are stable");
        check(std::string(restartGenerationStoreStatusName(
            RestartGenerationStoreStatus::STORE_ROOT_CHANGED))
                =="STORE_ROOT_CHANGED",
            "root-identity failure status name is stable");
        check(std::string(restartGenerationStoreStatusName(
            RestartGenerationStoreStatus::COMMIT_STATE_UNCERTAIN))
                =="COMMIT_STATE_UNCERTAIN",
            "commit-uncertain status name is stable");
        check(RestartGenerationStore::ownerPayloadFileName(
            RestartOwner::FLUID)=="00-fluid.payload",
            "fluid owner filename binds roster order and identity");
        check(RestartGenerationStore::ownerPayloadFileName(
            RestartOwner::TRANSPORT)=="06-transport.payload",
            "transport owner filename binds the seventh owner");
        check(RestartGenerationStore::ownerPayloadFileName(
            RestartOwner::COUNT).empty(),
            "invalid owner has no store filename");

        ScopedDirectory successRoot("success");
        const RestartGenerationStorePublishResult published=
            RestartGenerationStore::publish(
                successRoot.path(),baseline,expected);
        if(!published.success())
            std::cerr << "INITIAL_PUBLISH_FAILURE "
                      << restartGenerationStoreStatusName(published.status)
                      << ' ' << published.detail << '\n';
        check(published.success(),"valid generation publishes");
        check(published.manifestPath==manifestPath(successRoot.path(),baseline)
              &&fs::is_regular_file(published.manifestPath),
            "manifest is the final visible commit marker");
        check(fs::is_regular_file(payloadPath(
                successRoot.path(),baseline,RestartOwner::FLUID))
              &&fs::is_regular_file(payloadPath(
                successRoot.path(),baseline,RestartOwner::TURBULENCE)),
            "all state owner files are published before completion");
        check(!fs::exists(payloadPath(
                successRoot.path(),baseline,RestartOwner::CHT))
              &&!fs::exists(payloadPath(
                successRoot.path(),baseline,RestartOwner::CHEMISTRY)),
            "disabled and identity-only owners have no payload files");

        const RestartGenerationStoreLoadResult loaded=
            RestartGenerationStore::load(
                successRoot.path(),baseline.manifest.generationId,expected);
        check(loaded.success(),"committed generation loads");
        const RestartCodecEncodeResult baselineBytes=
            RestartManifestCodec::encode(baseline,expected);
        const RestartCodecEncodeResult loadedBytes=loaded.success()
            ?RestartManifestCodec::encode(*loaded.candidate,expected)
            :RestartCodecEncodeResult{};
        check(baselineBytes.success()&&loadedBytes.success()
              &&baselineBytes.bytes==loadedBytes.bytes,
            "loaded external owner package is byte-identical to codec input");
        check(readBytes(published.manifestPath)==baselineBytes.bytes,
            "commit marker contains the exact accepted v2 codec bytes");

#ifdef __linux__
        const fs::path successTrailingSeparator(
            successRoot.path().string()+"/");
        check(RestartGenerationStore::load(
                successTrailingSeparator,baseline.manifest.generationId,
                expected).success(),
            "a trailing separator does not weaken a direct-root load");
        const RestartCandidatePackage trailingSeparatorCandidate=
            validCandidate(
                "generation-direct-root-trailing-separator",
                "direct-root-trailing-separator");
        check(RestartGenerationStore::publish(
                successTrailingSeparator,trailingSeparatorCandidate,
                expectation(trailingSeparatorCandidate)).success(),
            "a trailing separator preserves direct-root publish semantics");
#endif

        const auto beforeConflict=directorySnapshot(
            generationPath(successRoot.path(),baseline));
        const RestartGenerationStorePublishResult conflict=
            RestartGenerationStore::publish(
                successRoot.path(),baseline,expected);
        check(!conflict.success()
              &&conflict.status
                ==RestartGenerationStoreStatus::GENERATION_CONFLICT,
            "an existing committed generation is immutable");
        check(beforeConflict==directorySnapshot(
                generationPath(successRoot.path(),baseline)),
            "conflicting publish changes no committed byte");

        fs::path manifestTemp=manifestPath(successRoot.path(),baseline);
        manifestTemp += ".tmp";
        writeBytes(manifestTemp,"stale-manifest-temp");
        fs::path payloadTemp=payloadPath(
            successRoot.path(),baseline,RestartOwner::FLUID);
        payloadTemp += ".tmp";
        writeBytes(payloadTemp,"stale-owner-temp");
        check(RestartGenerationStore::load(
                successRoot.path(),baseline.manifest.generationId,expected)
                .success(),
            "canonical stale temporary files are ignored after commit");
        fs::remove(manifestTemp);
        fs::remove(payloadTemp);

        const fs::path unknown=generationPath(successRoot.path(),baseline)
            /"unknown.payload";
        writeBytes(unknown,"unexpected");
        const RestartGenerationStoreLoadResult unknownLayout=
            RestartGenerationStore::load(
                successRoot.path(),baseline.manifest.generationId,expected);
        checkNoCandidate(unknownLayout,
            "unexpected generation file fails closed");
        check(unknownLayout.status
                ==RestartGenerationStoreStatus::FILE_LAYOUT_REJECTED,
            "unexpected generation file has exact layout status");
        fs::remove(unknown);

        ScopedDirectory deterministicRoot("deterministic");
        const RestartGenerationStorePublishResult deterministicPublish=
            RestartGenerationStore::publish(
                deterministicRoot.path(),baseline,expected);
        check(deterministicPublish.success(),
            "identical candidate publishes in an independent root");
        check(directorySnapshot(generationPath(successRoot.path(),baseline))
                ==directorySnapshot(
                    generationPath(deterministicRoot.path(),baseline)),
            "final generation bytes are deterministic across store roots");

        ScopedDirectory staleRoot("stale-before-publish");
        const fs::path staleGeneration=
            generationPath(staleRoot.path(),baseline);
        fs::create_directory(staleGeneration);
        fs::path stalePayloadTemp=payloadPath(
            staleRoot.path(),baseline,RestartOwner::FLUID);
        stalePayloadTemp += ".tmp";
        writeBytes(stalePayloadTemp,"orphaned-temp");
        checkNoCandidate(RestartGenerationStore::load(
                staleRoot.path(),baseline.manifest.generationId,expected),
            "temporary-only generation is not committed");
        const RestartGenerationStorePublishResult staleReuse=
            RestartGenerationStore::publish(
                staleRoot.path(),baseline,expected);
        check(!staleReuse.success()
              &&staleReuse.status
                ==RestartGenerationStoreStatus::INCOMPLETE_GENERATION,
            "publisher never reuses a stale generation directory");
        check(readBytes(stalePayloadTemp)=="orphaned-temp",
            "failed generation reuse does not modify stale bytes");

        ScopedDirectory missingRoot("missing-root");
        const RestartGenerationStorePublishResult missingStore=
            RestartGenerationStore::publish(
                missingRoot.path()/"not-created",baseline,expected);
        check(!missingStore.success()
              &&missingStore.status
                ==RestartGenerationStoreStatus::STORE_ROOT_UNAVAILABLE,
            "store root must already exist");

        ScopedDirectory symlinkTarget("root-symlink-target");
        check(RestartGenerationStore::publish(
                symlinkTarget.path(),baseline,expected).success(),
            "root-symlink target fixture publishes directly");
        fs::path symlinkRoot=symlinkTarget.path();
        symlinkRoot += "-link";
        std::error_code symlinkError;
        fs::create_directory_symlink(
            symlinkTarget.path(),symlinkRoot,symlinkError);
        check(!symlinkError,
            "root-symlink negative fixture is created");
        const RestartGenerationStoreLoadResult symlinkLoad=
            RestartGenerationStore::load(
                symlinkRoot,baseline.manifest.generationId,expected);
        checkNoCandidate(symlinkLoad,
            "load rejects a symlinked store root before manifest admission");
        check(symlinkLoad.status
                ==RestartGenerationStoreStatus::STORE_ROOT_UNAVAILABLE,
            "symlinked load root has the exact admission status");
        const RestartCandidatePackage symlinkPublishCandidate=validCandidate(
            "generation-store-root-symlink-publish","root-symlink");
        const RestartGenerationStorePublishResult symlinkPublish=
            RestartGenerationStore::publish(
                symlinkRoot,symlinkPublishCandidate,
                expectation(symlinkPublishCandidate));
        check(!symlinkPublish.success()
              &&symlinkPublish.status
                ==RestartGenerationStoreStatus::STORE_ROOT_UNAVAILABLE,
            "publish rejects a symlinked store root before generation lookup");
        check(!fs::exists(generationPath(
                symlinkTarget.path(),symlinkPublishCandidate)),
            "symlinked root rejection creates no generation");

#ifdef __linux__
        const std::vector<fs::path> symlinkAliases={
            fs::path(symlinkRoot.string()+"/"),
            symlinkRoot/".",
            symlinkRoot/"child"/".."
        };
        fs::create_directory(symlinkTarget.path()/"child");
        for(std::size_t aliasIndex=0U;
            aliasIndex<symlinkAliases.size();++aliasIndex)
        {
            const std::string label="symlink-root-alias-"
                +std::to_string(aliasIndex);
            const RestartGenerationStoreLoadResult aliasLoad=
                RestartGenerationStore::load(
                    symlinkAliases[aliasIndex],
                    baseline.manifest.generationId,expected);
            checkNoCandidate(aliasLoad,label+" load exposes no candidate");
            check(aliasLoad.status
                    ==RestartGenerationStoreStatus::STORE_ROOT_UNAVAILABLE,
                label+" load fails root admission");

            const RestartCandidatePackage aliasCandidate=validCandidate(
                "generation-"+label,label);
            const RestartGenerationStorePublishResult aliasPublish=
                RestartGenerationStore::publish(
                    symlinkAliases[aliasIndex],aliasCandidate,
                    expectation(aliasCandidate));
            check(!aliasPublish.success()
                  &&aliasPublish.status
                    ==RestartGenerationStoreStatus::STORE_ROOT_UNAVAILABLE,
                label+" publish fails root admission");
            check(!fs::exists(generationPath(
                    symlinkTarget.path(),aliasCandidate)),
                label+" publish creates no target generation");
        }

        ScopedDirectory ancestorTarget("root-symlink-ancestor-target");
        const fs::path nestedStore=ancestorTarget.path()/"nested-store";
        fs::create_directory(nestedStore);
        check(RestartGenerationStore::publish(
                nestedStore,baseline,expected).success(),
            "symlinked-ancestor target fixture publishes directly");
        fs::path ancestorLink=ancestorTarget.path();
        ancestorLink += "-link";
        fs::create_directory_symlink(
            ancestorTarget.path(),ancestorLink,symlinkError);
        check(!symlinkError,
            "symlinked-ancestor negative fixture is created");
        const fs::path nestedAlias=ancestorLink/"nested-store";
        const RestartGenerationStoreLoadResult ancestorLoad=
            RestartGenerationStore::load(
                nestedAlias,baseline.manifest.generationId,expected);
        checkNoCandidate(ancestorLoad,
            "load rejects a symlinked store-root ancestor");
        check(ancestorLoad.status
                ==RestartGenerationStoreStatus::STORE_ROOT_UNAVAILABLE,
            "symlinked-ancestor load has exact admission status");
        const RestartCandidatePackage ancestorCandidate=validCandidate(
            "generation-store-symlinked-ancestor","symlinked-ancestor");
        const RestartGenerationStorePublishResult ancestorPublish=
            RestartGenerationStore::publish(
                nestedAlias,ancestorCandidate,expectation(ancestorCandidate));
        check(!ancestorPublish.success()
              &&ancestorPublish.status
                ==RestartGenerationStoreStatus::STORE_ROOT_UNAVAILABLE,
            "publish rejects a symlinked store-root ancestor");
        check(!fs::exists(generationPath(nestedStore,ancestorCandidate)),
            "symlinked-ancestor publish creates no target generation");
        fs::remove(ancestorLink,symlinkError);
#endif
        fs::remove(symlinkRoot,symlinkError);

        RestartCandidatePackage unsafe=baseline;
        unsafe.manifest.generationId="../escape";
        for(RestartOwnerRecord& owner:unsafe.manifest.owners)
            owner.generationId=unsafe.manifest.generationId;
        for(RestartPayloadCandidate& payload:unsafe.payloads)
            payload.generationId=unsafe.manifest.generationId;
        unsafe.manifest.generationChecksumSha256=
            restartManifestChecksum(unsafe.manifest);
        check(RestartGenerationStore::publish(
                missingRoot.path(),unsafe,expectation(unsafe)).status
                ==RestartGenerationStoreStatus::UNSAFE_GENERATION_ID,
            "path-separator generation ID is rejected before file creation");
        check(RestartGenerationStore::load(
                missingRoot.path(),"../escape",expected).status
                ==RestartGenerationStoreStatus::UNSAFE_GENERATION_ID,
            "unsafe requested generation is rejected before lookup");

        ScopedDirectory incompleteRoot("incomplete");
        fs::create_directory(generationPath(incompleteRoot.path(),baseline));
        writeBytes(payloadPath(
            incompleteRoot.path(),baseline,RestartOwner::FLUID),"orphan");
        const RestartGenerationStorePublishResult incomplete=
            RestartGenerationStore::publish(
                incompleteRoot.path(),baseline,expected);
        check(!incomplete.success()
              &&incomplete.status
                ==RestartGenerationStoreStatus::INCOMPLETE_GENERATION,
            "orphan final owner file blocks generation reuse");
        checkNoCandidate(RestartGenerationStore::load(
                incompleteRoot.path(),baseline.manifest.generationId,expected),
            "orphan final owner file is not a candidate");

        const std::vector<std::pair<RestartGenerationStoreFaultPoint,
            std::pair<std::string,RestartGenerationStoreStatus>>>
            payloadFaults={
                {RestartGenerationStoreFaultPoint::PAYLOAD_OPEN,
                    {"payload-open",RestartGenerationStoreStatus::STAGE_OPEN_FAILED}},
                {RestartGenerationStoreFaultPoint::PAYLOAD_WRITE,
                    {"payload-write",RestartGenerationStoreStatus::STAGE_WRITE_FAILED}},
                {RestartGenerationStoreFaultPoint::PAYLOAD_FLUSH,
                    {"payload-flush",RestartGenerationStoreStatus::STAGE_FLUSH_FAILED}},
                {RestartGenerationStoreFaultPoint::PAYLOAD_CLOSE,
                    {"payload-close",RestartGenerationStoreStatus::STAGE_CLOSE_FAILED}},
                {RestartGenerationStoreFaultPoint::PAYLOAD_FILE_SYNC,
                    {"payload-file-sync",RestartGenerationStoreStatus::FILE_SYNC_FAILED}},
                {RestartGenerationStoreFaultPoint::PAYLOAD_READBACK,
                    {"payload-readback",RestartGenerationStoreStatus::STAGE_READBACK_FAILED}},
                {RestartGenerationStoreFaultPoint::PAYLOAD_RENAME,
                    {"payload-rename",RestartGenerationStoreStatus::ATOMIC_RENAME_FAILED}},
                {RestartGenerationStoreFaultPoint::PAYLOAD_DIRECTORY_SYNC,
                    {"payload-directory-sync",RestartGenerationStoreStatus::DIRECTORY_SYNC_FAILED}}
            };
        for(const auto& fault:payloadFaults)
            checkFault(fault.first,RestartOwner::TURBULENCE,
                fault.second.first,fault.second.second);

        const std::vector<std::pair<RestartGenerationStoreFaultPoint,
            std::pair<std::string,RestartGenerationStoreStatus>>>
            manifestFaults={
                {RestartGenerationStoreFaultPoint::MANIFEST_OPEN,
                    {"manifest-open",RestartGenerationStoreStatus::STAGE_OPEN_FAILED}},
                {RestartGenerationStoreFaultPoint::MANIFEST_WRITE,
                    {"manifest-write",RestartGenerationStoreStatus::STAGE_WRITE_FAILED}},
                {RestartGenerationStoreFaultPoint::MANIFEST_FLUSH,
                    {"manifest-flush",RestartGenerationStoreStatus::STAGE_FLUSH_FAILED}},
                {RestartGenerationStoreFaultPoint::MANIFEST_CLOSE,
                    {"manifest-close",RestartGenerationStoreStatus::STAGE_CLOSE_FAILED}},
                {RestartGenerationStoreFaultPoint::MANIFEST_FILE_SYNC,
                    {"manifest-file-sync",RestartGenerationStoreStatus::FILE_SYNC_FAILED}},
                {RestartGenerationStoreFaultPoint::MANIFEST_READBACK,
                    {"manifest-readback",RestartGenerationStoreStatus::STAGE_READBACK_FAILED}},
                {RestartGenerationStoreFaultPoint::MANIFEST_RENAME,
                    {"manifest-rename",RestartGenerationStoreStatus::ATOMIC_RENAME_FAILED}},
                {RestartGenerationStoreFaultPoint::MANIFEST_DIRECTORY_SYNC,
                    {"manifest-directory-sync",RestartGenerationStoreStatus::DIRECTORY_SYNC_FAILED}}
            };
        for(const auto& fault:manifestFaults)
            checkFault(fault.first,RestartOwner::FLUID,
                fault.second.first,fault.second.second);

        ScopedDirectory rollbackRemoveRoot("rollback-remove-failure");
        const RestartCandidatePackage rollbackRemoveCandidate=validCandidate(
            "generation-rollback-remove-failure","rollback-remove");
        const RestartManifestExpectation rollbackRemoveExpected=
            expectation(rollbackRemoveCandidate);
        RestartGenerationStoreFault rollbackRemoveFault;
        rollbackRemoveFault.point=
            RestartGenerationStoreFaultPoint::MANIFEST_ROLLBACK_REMOVE;
        const RestartGenerationStorePublishResult rollbackRemove=
            RestartGenerationStore::publish(
                rollbackRemoveRoot.path(),rollbackRemoveCandidate,
                rollbackRemoveExpected,rollbackRemoveFault);
        check(!rollbackRemove.success()
              &&rollbackRemove.status
                ==RestartGenerationStoreStatus::COMMIT_STATE_UNCERTAIN,
            "marker-cleanup failure returns commit-state-uncertain");
        check(fs::is_regular_file(manifestPath(
                rollbackRemoveRoot.path(),rollbackRemoveCandidate)),
            "marker-cleanup failure may leave the final marker visible");
        check(RestartGenerationStore::load(
                rollbackRemoveRoot.path(),
                rollbackRemoveCandidate.manifest.generationId,
                rollbackRemoveExpected).success(),
            "visible marker after cleanup failure is never reported uncommitted");

        ScopedDirectory rollbackSyncRoot("rollback-sync-failure");
        const RestartCandidatePackage rollbackSyncCandidate=validCandidate(
            "generation-rollback-sync-failure","rollback-sync");
        const RestartManifestExpectation rollbackSyncExpected=
            expectation(rollbackSyncCandidate);
        RestartGenerationStoreFault rollbackSyncFault;
        rollbackSyncFault.point=RestartGenerationStoreFaultPoint::
            MANIFEST_ROLLBACK_DIRECTORY_SYNC;
        const RestartGenerationStorePublishResult rollbackSync=
            RestartGenerationStore::publish(
                rollbackSyncRoot.path(),rollbackSyncCandidate,
                rollbackSyncExpected,rollbackSyncFault);
        check(!rollbackSync.success()
              &&rollbackSync.status
                ==RestartGenerationStoreStatus::COMMIT_STATE_UNCERTAIN,
            "cleanup-directory-sync failure returns commit-state-uncertain");
        check(!fs::exists(manifestPath(
                rollbackSyncRoot.path(),rollbackSyncCandidate)),
            "cleanup-directory-sync failure removes the live marker");
        const RestartGenerationStoreLoadResult rollbackSyncLoad=
            RestartGenerationStore::load(
                rollbackSyncRoot.path(),
                rollbackSyncCandidate.manifest.generationId,
                rollbackSyncExpected);
        checkNoCandidate(rollbackSyncLoad,
            "cleanup-directory-sync failure exposes no live candidate");
        check(rollbackSyncLoad.status
                ==RestartGenerationStoreStatus::NOT_COMMITTED,
            "cleanup-directory-sync failure remains crash-persistence uncertain");

        ScopedDirectory rollbackSwapRoot("rollback-root-identity-swap");
        const RestartCandidatePackage rollbackSwapCandidate=validCandidate(
            "generation-rollback-root-identity-swap","rollback-root-swap");
        const RestartManifestExpectation rollbackSwapExpected=
            expectation(rollbackSwapCandidate);
        fs::path rollbackMovedRoot=rollbackSwapRoot.path();
        rollbackMovedRoot += "-retained";
        std::error_code rollbackSwapCleanupError;
        fs::remove_all(rollbackMovedRoot,rollbackSwapCleanupError);
        RootSwapContext rollbackSwapContext{
            rollbackSwapRoot.path(),rollbackMovedRoot,{},false};
        RestartGenerationStoreFault rollbackSwapFault;
        rollbackSwapFault.point=
            RestartGenerationStoreFaultPoint::MANIFEST_DIRECTORY_SYNC;
        rollbackSwapFault.beforeFinalIdentityCheck=replaceRootIdentity;
        rollbackSwapFault.hookContext=&rollbackSwapContext;
        const RestartGenerationStorePublishResult rollbackSwap=
            RestartGenerationStore::publish(
                rollbackSwapRoot.path(),rollbackSwapCandidate,
                rollbackSwapExpected,rollbackSwapFault);
        check(rollbackSwapContext.swapped&&!rollbackSwapContext.error,
            "post-marker rollback test hook replaces the admitted root");
        check(!rollbackSwap.success()
              &&rollbackSwap.status
                ==RestartGenerationStoreStatus::COMMIT_STATE_UNCERTAIN,
            "directory-sync failure plus namespace swap is commit-uncertain");
        check(!fs::exists(manifestPath(
                rollbackMovedRoot,rollbackSwapCandidate)),
            "identity-changing rollback removed the retained-root marker");
        const RestartGenerationStoreLoadResult rollbackSwapLoad=
            RestartGenerationStore::load(
                rollbackMovedRoot,
                rollbackSwapCandidate.manifest.generationId,
                rollbackSwapExpected);
        checkNoCandidate(rollbackSwapLoad,
            "identity-changing rollback exposes no retained-root candidate");
        check(rollbackSwapLoad.status
                ==RestartGenerationStoreStatus::NOT_COMMITTED,
            "identity-changing rollback leaves no live commit marker");
        fs::remove_all(rollbackMovedRoot,rollbackSwapCleanupError);

        ScopedDirectory orderRoot("commit-order");
        RestartGenerationStoreFault manifestRename;
        manifestRename.point=
            RestartGenerationStoreFaultPoint::MANIFEST_RENAME;
        const RestartGenerationStorePublishResult orderFailure=
            RestartGenerationStore::publish(
                orderRoot.path(),baseline,expected,manifestRename);
        check(!orderFailure.success()
              &&fs::is_regular_file(payloadPath(
                orderRoot.path(),baseline,RestartOwner::FLUID))
              &&fs::is_regular_file(payloadPath(
                orderRoot.path(),baseline,RestartOwner::TURBULENCE)),
            "all owner payloads are final before manifest rename is attempted");
        fs::path stagedManifest=manifestPath(orderRoot.path(),baseline);
        stagedManifest += ".tmp";
        check(fs::is_regular_file(stagedManifest)
              &&stagedManifest.parent_path()
                ==manifestPath(orderRoot.path(),baseline).parent_path(),
            "manifest staging and final paths share one directory");

        ScopedDirectory corruptionRoot("corruption");
        check(RestartGenerationStore::publish(
                corruptionRoot.path(),baseline,expected).success(),
            "corruption fixture publishes");
        const fs::path fluidPath=payloadPath(
            corruptionRoot.path(),baseline,RestartOwner::FLUID);
        const std::string fluidOriginal=readBytes(fluidPath);
        writeBytes(fluidPath,fluidOriginal.substr(0,fluidOriginal.size()-1U));
        const RestartGenerationStoreLoadResult truncatedPayload=
            RestartGenerationStore::load(
                corruptionRoot.path(),baseline.manifest.generationId,expected);
        checkNoCandidate(truncatedPayload,"truncated owner payload fails closed");
        check(truncatedPayload.status
                ==RestartGenerationStoreStatus::OWNER_SIZE_MISMATCH,
            "truncated owner payload reports size mismatch");
        std::string corruptFluid=fluidOriginal;
        corruptFluid.back()=static_cast<char>(corruptFluid.back()^0x01);
        writeBytes(fluidPath,corruptFluid);
        const RestartGenerationStoreLoadResult corruptPayload=
            RestartGenerationStore::load(
                corruptionRoot.path(),baseline.manifest.generationId,expected);
        checkNoCandidate(corruptPayload,"same-size owner corruption fails closed");
        check(corruptPayload.status
                ==RestartGenerationStoreStatus::OWNER_CHECKSUM_MISMATCH,
            "same-size owner corruption reports SHA256 mismatch");
        writeBytes(fluidPath,fluidOriginal);

        const fs::path committedManifest=manifestPath(corruptionRoot.path(),baseline);
        const std::string manifestOriginal=readBytes(committedManifest);
        writeBytes(committedManifest,
            manifestOriginal.substr(0,manifestOriginal.size()-1U));
        checkNoCandidate(RestartGenerationStore::load(
                corruptionRoot.path(),baseline.manifest.generationId,expected),
            "truncated manifest fails exact codec decode");
        std::string corruptManifest=manifestOriginal;
        corruptManifest[corruptManifest.size()/2U]=static_cast<char>(
            corruptManifest[corruptManifest.size()/2U]^0x01);
        writeBytes(committedManifest,corruptManifest);
        checkNoCandidate(RestartGenerationStore::load(
                corruptionRoot.path(),baseline.manifest.generationId,expected),
            "single-bit manifest corruption fails closed");
        writeBytes(committedManifest,manifestOriginal);

        const fs::path missingPayload=payloadPath(
            corruptionRoot.path(),baseline,RestartOwner::TURBULENCE);
        const std::string turbulenceOriginal=readBytes(missingPayload);
        fs::remove(missingPayload);
        const RestartGenerationStoreLoadResult missingOwner=
            RestartGenerationStore::load(
                corruptionRoot.path(),baseline.manifest.generationId,expected);
        checkNoCandidate(missingOwner,"missing state owner fails closed");
        check(missingOwner.status
                ==RestartGenerationStoreStatus::OWNER_FILE_MISSING,
            "missing state owner reports exact status");
        writeBytes(missingPayload,turbulenceOriginal);

        const fs::path unexpectedOwner=payloadPath(
            corruptionRoot.path(),baseline,RestartOwner::CHT);
        writeBytes(unexpectedOwner,"unexpected-disabled-owner");
        checkNoCandidate(RestartGenerationStore::load(
                corruptionRoot.path(),baseline.manifest.generationId,expected),
            "payload for disabled owner fails closed");
        fs::remove(unexpectedOwner);

        ScopedDirectory mixedRoot("mixed-source");
        const RestartCandidatePackage other=validCandidate(
            "generation-store-20260727-0002","baseline");
        const RestartManifestExpectation otherExpected=expectation(other);
        check(RestartGenerationStore::publish(
                mixedRoot.path(),other,otherExpected).success(),
            "second generation fixture publishes");
        writeBytes(fluidPath,readBytes(payloadPath(
            mixedRoot.path(),other,RestartOwner::FLUID)));
        const RestartGenerationStoreLoadResult mixedPayload=
            RestartGenerationStore::load(
                corruptionRoot.path(),baseline.manifest.generationId,expected);
        checkNoCandidate(mixedPayload,"mixed-generation owner file fails closed");
        check(mixedPayload.status
                ==RestartGenerationStoreStatus::OWNER_GENERATION_MISMATCH,
            "identical raw payload from another generation is rejected");
        writeBytes(fluidPath,fluidOriginal);

        writeBytes(committedManifest,readBytes(manifestPath(mixedRoot.path(),other)));
        const RestartGenerationStoreLoadResult mixedManifest=
            RestartGenerationStore::load(
                corruptionRoot.path(),baseline.manifest.generationId,expected);
        checkNoCandidate(mixedManifest,"mixed-generation manifest fails closed");
        check(mixedManifest.status
                ==RestartGenerationStoreStatus::GENERATION_MISMATCH,
            "mixed-generation manifest is rejected by requested generation");
        writeBytes(committedManifest,manifestOriginal);

        RestartManifestExpectation wrongExpected=expected;
        wrongExpected.identities.buildSha256=hash("different-build");
        checkNoCandidate(RestartGenerationStore::load(
                corruptionRoot.path(),baseline.manifest.generationId,
                wrongExpected),
            "independently supplied expectation mismatch fails closed");
        check(RestartGenerationStore::load(
                corruptionRoot.path(),baseline.manifest.generationId,expected)
                .success(),
            "restored generation remains readable after negative tests");

        const fs::path hardLinkSource=
            corruptionRoot.path()/"fluid-hard-link-source";
        fs::rename(fluidPath,hardLinkSource);
        fs::create_hard_link(hardLinkSource,fluidPath);
        checkNoCandidate(RestartGenerationStore::load(
                corruptionRoot.path(),baseline.manifest.generationId,expected),
            "hard-linked owner substitution fails stable-handle admission");
        fs::remove(fluidPath);
        fs::rename(hardLinkSource,fluidPath);

        const fs::path symlinkSource=
            corruptionRoot.path()/"fluid-symlink-source";
        fs::rename(fluidPath,symlinkSource);
        fs::create_symlink(symlinkSource,fluidPath);
        checkNoCandidate(RestartGenerationStore::load(
                corruptionRoot.path(),baseline.manifest.generationId,expected),
            "symlinked owner substitution fails no-follow admission");
        fs::remove(fluidPath);
        fs::rename(symlinkSource,fluidPath);

        const fs::path manifestHardLinkSource=
            corruptionRoot.path()/"manifest-hard-link-source";
        fs::rename(committedManifest,manifestHardLinkSource);
        fs::create_hard_link(manifestHardLinkSource,committedManifest);
        const RestartGenerationStoreLoadResult hardLinkedManifest=
            RestartGenerationStore::load(
                corruptionRoot.path(),baseline.manifest.generationId,expected);
        checkNoCandidate(hardLinkedManifest,
            "hard-linked manifest substitution fails stable-handle admission");
        check(hardLinkedManifest.status
                ==RestartGenerationStoreStatus::MANIFEST_READ_FAILED,
            "hard-linked manifest has exact read failure status");
        fs::remove(committedManifest);
        fs::rename(manifestHardLinkSource,committedManifest);

        ScopedDirectory publishSwapRoot("publish-root-identity-swap");
        const RestartCandidatePackage publishSwapCandidate=validCandidate(
            "generation-publish-root-identity-swap","publish-root-swap");
        const RestartManifestExpectation publishSwapExpected=
            expectation(publishSwapCandidate);
        fs::path publishMovedRoot=publishSwapRoot.path();
        publishMovedRoot += "-retained";
        std::error_code swapCleanupError;
        fs::remove_all(publishMovedRoot,swapCleanupError);
        RootSwapContext publishSwapContext{
            publishSwapRoot.path(),publishMovedRoot,{},false};
        RestartGenerationStoreFault publishSwapFault;
        publishSwapFault.beforeFinalIdentityCheck=replaceRootIdentity;
        publishSwapFault.hookContext=&publishSwapContext;
        const RestartGenerationStorePublishResult publishSwap=
            RestartGenerationStore::publish(
                publishSwapRoot.path(),publishSwapCandidate,
                publishSwapExpected,publishSwapFault);
        check(publishSwapContext.swapped&&!publishSwapContext.error,
            "publish test hook replaces the admitted root identity");
        check(!publishSwap.success()
              &&publishSwap.status
                ==RestartGenerationStoreStatus::COMMIT_STATE_UNCERTAIN,
            "publish never reports success after a root identity change");
        check(RestartGenerationStore::load(
                publishMovedRoot,publishSwapCandidate.manifest.generationId,
                publishSwapExpected).success(),
            "root-anchored publish remains complete under the retained root");
        checkNoCandidate(RestartGenerationStore::load(
                publishSwapRoot.path(),
                publishSwapCandidate.manifest.generationId,
                publishSwapExpected),
            "replacement root does not inherit the anchored generation");
        fs::remove_all(publishMovedRoot,swapCleanupError);

        ScopedDirectory loadSwapRoot("load-root-identity-swap");
        const RestartCandidatePackage loadSwapCandidate=validCandidate(
            "generation-load-root-identity-swap","load-root-swap");
        const RestartManifestExpectation loadSwapExpected=
            expectation(loadSwapCandidate);
        check(RestartGenerationStore::publish(
                loadSwapRoot.path(),loadSwapCandidate,loadSwapExpected)
                .success(),
            "load root-swap fixture publishes");
        fs::path loadMovedRoot=loadSwapRoot.path();
        loadMovedRoot += "-retained";
        fs::remove_all(loadMovedRoot,swapCleanupError);
        RootSwapContext loadSwapContext{
            loadSwapRoot.path(),loadMovedRoot,{},false};
        RestartGenerationStoreFault loadSwapFault;
        loadSwapFault.beforeFinalIdentityCheck=replaceRootIdentity;
        loadSwapFault.hookContext=&loadSwapContext;
        const RestartGenerationStoreLoadResult loadSwap=
            RestartGenerationStore::load(
                loadSwapRoot.path(),loadSwapCandidate.manifest.generationId,
                loadSwapExpected,loadSwapFault);
        check(loadSwapContext.swapped&&!loadSwapContext.error,
            "load test hook replaces the admitted root identity");
        checkNoCandidate(loadSwap,
            "load exposes no candidate after a root identity change");
        check(loadSwap.status
                ==RestartGenerationStoreStatus::STORE_ROOT_CHANGED,
            "load root identity change has the exact structured status");
        check(RestartGenerationStore::load(
                loadMovedRoot,loadSwapCandidate.manifest.generationId,
                loadSwapExpected).success(),
            "retained root remains independently readable after load rejection");
        fs::remove_all(loadMovedRoot,swapCleanupError);

        const auto acceptedBeforeFailure=directorySnapshot(
            generationPath(successRoot.path(),baseline));
        const RestartCandidatePackage failedNewGeneration=validCandidate(
            "generation-store-20260727-failed-new","failed-new");
        const RestartManifestExpectation failedNewExpected=
            expectation(failedNewGeneration);
        RestartGenerationStoreFault failNewMarker;
        failNewMarker.point=
            RestartGenerationStoreFaultPoint::MANIFEST_RENAME;
        check(!RestartGenerationStore::publish(
                successRoot.path(),failedNewGeneration,failedNewExpected,
                failNewMarker).success(),
            "new generation manifest fault is reported");
        check(acceptedBeforeFailure==directorySnapshot(
                generationPath(successRoot.path(),baseline)),
            "failed new generation never modifies prior accepted generation");
        check(RestartGenerationStore::load(
                successRoot.path(),baseline.manifest.generationId,expected)
                .success(),
            "prior accepted generation remains readable after new write failure");

        const std::vector<std::pair<RestartGenerationStoreFaultPoint,
            std::string>> immutableBoundaryFaults={
                {RestartGenerationStoreFaultPoint::PAYLOAD_FILE_SYNC,
                    "old-safe-payload-file-sync"},
                {RestartGenerationStoreFaultPoint::PAYLOAD_RENAME,
                    "old-safe-payload-rename"},
                {RestartGenerationStoreFaultPoint::PAYLOAD_DIRECTORY_SYNC,
                    "old-safe-payload-dir-sync"},
                {RestartGenerationStoreFaultPoint::MANIFEST_FILE_SYNC,
                    "old-safe-manifest-file-sync"},
                {RestartGenerationStoreFaultPoint::MANIFEST_DIRECTORY_SYNC,
                    "old-safe-manifest-dir-sync"}
            };
        for(const auto& boundary:immutableBoundaryFaults)
        {
            const RestartCandidatePackage attempted=validCandidate(
                "generation-store-20260727-"+boundary.second,boundary.second);
            RestartGenerationStoreFault fault;
            fault.point=boundary.first;
            fault.owner=RestartOwner::TURBULENCE;
            check(!RestartGenerationStore::publish(
                    successRoot.path(),attempted,expectation(attempted),fault)
                    .success(),
                boundary.second+" rejects the new generation");
            check(acceptedBeforeFailure==directorySnapshot(
                    generationPath(successRoot.path(),baseline)),
                boundary.second+" preserves every prior accepted byte");
            check(RestartGenerationStore::load(
                    successRoot.path(),baseline.manifest.generationId,expected)
                    .success(),
                boundary.second+" leaves the prior generation readable");
        }

        ScopedDirectory concurrentRoot("concurrent");
        const RestartCandidatePackage concurrentCandidate=validCandidate(
            "generation-store-20260727-concurrent","concurrent");
        const RestartManifestExpectation concurrentExpected=
            expectation(concurrentCandidate);
        std::atomic<bool> start{false};
        RestartGenerationStorePublishResult concurrentPublish;
        RestartGenerationStoreLoadResult concurrentLoad;
        std::thread writer([&]()
        {
            while(!start.load(std::memory_order_acquire)) {}
            concurrentPublish=RestartGenerationStore::publish(
                concurrentRoot.path(),concurrentCandidate,concurrentExpected);
        });
        std::thread reader([&]()
        {
            while(!start.load(std::memory_order_acquire)) {}
            concurrentLoad=RestartGenerationStore::load(
                concurrentRoot.path(),concurrentCandidate.manifest.generationId,
                concurrentExpected);
        });
        start.store(true,std::memory_order_release);
        writer.join();
        reader.join();
        check(concurrentPublish.success(),
            "serialized concurrent writer commits the generation");
        check(concurrentLoad.success()
              ||(concurrentLoad.status
                    ==RestartGenerationStoreStatus::NOT_COMMITTED
                 &&!concurrentLoad.candidate.has_value()),
            "concurrent reader observes either no marker or a complete candidate");
        check(RestartGenerationStore::load(
                concurrentRoot.path(),
                concurrentCandidate.manifest.generationId,
                concurrentExpected).success(),
            "candidate is complete after concurrent publication joins");

        std::cout << "RESTART_GENERATION_STORE_MANIFEST_SHA256 "
                  << hash(baselineBytes.bytes) << '\n';
    }
    catch(const std::exception& error)
    {
        std::cerr << "RESTART_GENERATION_STORE_TEST_EXCEPTION "
                  << error.what() << '\n';
        return 2;
    }

    if(failures!=0)
    {
        std::cerr << "RESTART_GENERATION_STORE_TEST_FAIL checks=" << checks
                  << " failures=" << failures << '\n';
        return 1;
    }
    std::cout << "RESTART_GENERATION_STORE_TEST_PASS checks=" << checks << '\n';
    return 0;
}
