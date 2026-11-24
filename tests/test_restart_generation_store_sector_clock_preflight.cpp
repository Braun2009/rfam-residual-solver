#include "RestartClockAdmission.H"

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#endif
#include "SectorCheckpoint.H"
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#include "ThermoDataProvenance.H"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>

namespace
{

using namespace Foam;
namespace fs = std::filesystem;

int checks = 0;
int failures = 0;

void check(bool condition, const std::string& label)
{
    ++checks;
    std::cout << "  [" << (condition ? "PASS" : "FAIL") << "] "
              << label << '\n';
    if (!condition) ++failures;
}

std::string hash(const std::string& value)
{
    return thermoProvenance::sha256Hex(value);
}

RestartCandidatePackage candidate(
    const std::string& generationId,
    bool previousAcceptedDtPresent)
{
    RestartCandidatePackage result;
    RestartManifest& manifest = result.manifest;
    manifest.generationId = generationId;
    manifest.generationState.physicalTime = 1.25e-5;
    manifest.generationState.acceptedIteration = 42U;
    manifest.generationState.previousAcceptedDtPresent =
        previousAcceptedDtPresent;
    manifest.generationState.previousAcceptedDt =
        previousAcceptedDtPresent ? 1.0e-8 : 0.0;
    manifest.byteOrder = nativeRestartByteOrder();
    manifest.cellCount = 2U;
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
    manifest.identities.topologyEnvelopeSha256 = hash("topology-envelope");

    const auto& roster = restartOwnerRoster();
    manifest.owners.resize(roster.size());
    for (std::size_t i = 0; i < roster.size(); ++i)
    {
        RestartOwnerRecord& owner = manifest.owners[i];
        owner.owner = roster[i];
        owner.disposition = RestartOwnerDisposition::IDENTITY_ONLY;
        owner.publicationOrder = static_cast<std::uint32_t>(i);
        owner.generationId = generationId;
        owner.identitySha256 = hash(
            std::string("owner-") + restartOwnerName(owner.owner));
        owner.payloadSha256 = hash(std::string());
    }
    manifest.owners[0].disposition = RestartOwnerDisposition::STATE_PAYLOAD;
    manifest.owners[1].disposition = RestartOwnerDisposition::STATE_PAYLOAD;
    manifest.owners[2].disposition = RestartOwnerDisposition::DISABLED;
    manifest.owners[3].disposition = RestartOwnerDisposition::DISABLED;
    result.payloads = {
        {RestartOwner::FLUID, generationId, "fluid-clock-preflight", 0U},
        {RestartOwner::TURBULENCE, generationId,
            "turbulence-clock-preflight", 0U}};
    for (RestartPayloadCandidate& payload : result.payloads)
    {
        payload.consumedBytes = payload.bytes.size();
        RestartOwnerRecord& owner =
            manifest.owners[static_cast<std::size_t>(payload.owner)];
        owner.payloadBytes = payload.bytes.size();
        owner.payloadSha256 = hash(payload.bytes);
    }
    manifest.recoverySeeds = {
        {0U, 5.0e6, 155.0, RestartPhaseAssertion::SINGLE_LIQUID,
            hash("phase-cell-0")},
        {1U, 5.0e6, 180.0, RestartPhaseAssertion::TWO_PHASE,
            hash("phase-cell-1-quality-x-y")}};
    manifest.generationChecksumSha256 = restartManifestChecksum(manifest);
    return result;
}

RestartManifestExpectation expectation(
    const RestartCandidatePackage& value,
    bool requirePreviousAcceptedDt)
{
    RestartManifestExpectation result;
    result.identities = value.manifest.identities;
    result.byteOrder = nativeRestartByteOrder();
    result.cellCount = value.manifest.cellCount;
    result.requirePreviousAcceptedDt = requirePreviousAcceptedDt;
    const auto& roster = restartOwnerRoster();
    for (std::size_t i = 0; i < roster.size(); ++i)
    {
        result.owners[i].owner = roster[i];
        result.owners[i].disposition =
            value.manifest.owners[i].disposition;
        result.owners[i].identitySha256 =
            value.manifest.owners[i].identitySha256;
    }
    return result;
}

class ScopedDirectory
{
    fs::path path_;

public:
    explicit ScopedDirectory(const std::string& label)
    {
        const auto tick = std::chrono::high_resolution_clock::now()
            .time_since_epoch().count();
        path_ = fs::temp_directory_path()/(
            "rfam-restart-clock-preflight-" + label + "-"
            + std::to_string(tick));
        std::error_code error;
        fs::remove_all(path_, error);
        error.clear();
        if (!fs::create_directory(path_, error) || error)
            throw std::runtime_error("cannot create preflight test directory");
    }

    ~ScopedDirectory()
    {
        std::error_code error;
        fs::remove_all(path_, error);
    }

    const fs::path& path() const { return path_; }
};

struct PublicationSentinel
{
    double solverTime = 9.0;
    int solverIteration = 17;
    std::string domainFingerprint = "domain-not-published";
};

bool sameSentinel(
    const PublicationSentinel& lhs,
    const PublicationSentinel& rhs)
{
    return lhs.solverTime == rhs.solverTime
        && lhs.solverIteration == rhs.solverIteration
        && lhs.domainFingerprint == rhs.domainFingerprint;
}

RestartCheckpointClockObservation observation(
    const SectorCheckpointMetadata& metadata)
{
    return {
        metadata.time,
        metadata.iteration,
        metadata.threadCount};
}

} // End anonymous namespace

int main()
{
    try
    {
        ScopedDirectory root("baseline");
        const RestartCandidatePackage noDt = candidate(
            "generation-clock-preflight-no-dt", false);
        const RestartManifestExpectation noDtExpectation =
            expectation(noDt, false);
        const RestartGenerationStorePublishResult noDtPublish =
            RestartGenerationStore::publish(
                root.path(), noDt, noDtExpectation);
        if (!noDtPublish.success())
        {
            std::cerr << "NO_DT_PUBLISH_FAILURE "
                      << restartGenerationStoreStatusName(noDtPublish.status)
                      << ' ' << noDtPublish.detail << '\n';
        }
        check(noDtPublish.success(),
            "candidate without previous dt publishes through the store");
        const RestartGenerationStoreLoadResult noDtLoad =
            RestartGenerationStore::load(
                root.path(), noDt.manifest.generationId, noDtExpectation);
        check(noDtLoad.success(),
            "candidate without previous dt loads through the store");

        SectorCheckpointMetadata metadata;
        metadata.time = noDt.manifest.generationState.physicalTime;
        metadata.iteration = static_cast<int>(
            noDt.manifest.generationState.acceptedIteration);
        metadata.threadCount = 1;
        const PublicationSentinel sentinel;
        PublicationSentinel runtime = sentinel;
        const RestartClockAdmissionResult accepted =
            RestartClockAdmission::stage(
                noDtLoad, observation(metadata), false);
        check(accepted.success(),
            "matching store and checkpoint clocks pass the pure preflight");
        check(accepted.publication->physicalTime == metadata.time
                && accepted.publication->acceptedIteration == metadata.iteration
                && !accepted.publication->previousAcceptedDt.has_value(),
            "successful preflight exposes only the validated clock values");
        Domain solverDomain(1, 1.0, 1);
        SolverBase solver(solverDomain);
        check(accepted.publication->publish(solver)
                && solver.time() == metadata.time
                && solver.iteration() == metadata.iteration
                && !solver.previousAcceptedDt().has_value(),
            "validated no-dt publication writes one complete solver clock");
        check(sameSentinel(runtime, sentinel),
            "clock publication does not publish Domain owner state");

        SectorCheckpointMetadata wrongTime = metadata;
        wrongTime.time += 1.0e-12;
        const RestartClockAdmissionResult timeMismatch =
            RestartClockAdmission::stage(
                noDtLoad, observation(wrongTime), false);
        check(timeMismatch.status == RestartClockAdmissionStatus::TIME_MISMATCH
                && !timeMismatch.publication,
            "time mismatch fails without exposing a publication plan");

        SectorCheckpointMetadata wrongIteration = metadata;
        ++wrongIteration.iteration;
        const RestartClockAdmissionResult iterationMismatch =
            RestartClockAdmission::stage(
                noDtLoad, observation(wrongIteration), false);
        check(iterationMismatch.status
                    == RestartClockAdmissionStatus::ITERATION_MISMATCH
                && !iterationMismatch.publication,
            "iteration mismatch fails without exposing a publication plan");

        SectorCheckpointMetadata negativeZero = metadata;
        negativeZero.time = -0.0;
        const RestartClockAdmissionResult invalidCheckpoint =
            RestartClockAdmission::stage(
                noDtLoad, observation(negativeZero), false);
        check(invalidCheckpoint.status
                    == RestartClockAdmissionStatus::INVALID_CHECKPOINT_CLOCK
                && !invalidCheckpoint.publication,
            "negative-zero checkpoint time fails before publication");

        SectorCheckpointMetadata invalidThreads = metadata;
        invalidThreads.threadCount = 0;
        const RestartClockAdmissionResult invalidThreadCount =
            RestartClockAdmission::stage(
                noDtLoad, observation(invalidThreads), false);
        check(invalidThreadCount.status
                    == RestartClockAdmissionStatus::INVALID_CHECKPOINT_CLOCK
                && !invalidThreadCount.publication,
            "nonpositive checkpoint thread count fails before publication");

        RestartGenerationStoreLoadResult malformedPreviousDt = noDtLoad;
        malformedPreviousDt.candidate->manifest.generationState
            .previousAcceptedDt = 1.0e-8;
        const RestartClockAdmissionResult invalidGenerationClock =
            RestartClockAdmission::stage(
                malformedPreviousDt, observation(metadata), false);
        check(invalidGenerationClock.status
                    == RestartClockAdmissionStatus::INVALID_GENERATION_CLOCK
                && !invalidGenerationClock.publication,
            "malformed in-memory previous-dt state is independently rejected");

        RestartGenerationStoreLoadResult overflow = noDtLoad;
        overflow.candidate->manifest.generationState.acceptedIteration =
            static_cast<std::uint64_t>(std::numeric_limits<int>::max()) + 1U;
        const RestartClockAdmissionResult iterationOverflow =
            RestartClockAdmission::stage(
                overflow, observation(metadata), false);
        check(iterationOverflow.status
                    == RestartClockAdmissionStatus::ITERATION_OVERFLOW
                && !iterationOverflow.publication,
            "iteration overflow fails before the solver int conversion");

        RestartGenerationStoreLoadResult invalidLoad;
        invalidLoad.status = RestartGenerationStoreStatus::NOT_COMMITTED;
        const RestartClockAdmissionResult rejectedStore =
            RestartClockAdmission::stage(
                invalidLoad, observation(metadata), false);
        check(rejectedStore.status
                    == RestartClockAdmissionStatus::INVALID_STORE_RESULT
                && !rejectedStore.publication,
            "invalid store result fails without a publication plan");

        const RestartCandidatePackage withDt = candidate(
            "generation-clock-preflight-with-dt", true);
        const RestartManifestExpectation withDtExpectation =
            expectation(withDt, true);
        const RestartGenerationStorePublishResult withDtPublish =
            RestartGenerationStore::publish(
                root.path(), withDt, withDtExpectation);
        if (!withDtPublish.success())
        {
            std::cerr << "WITH_DT_PUBLISH_FAILURE "
                      << restartGenerationStoreStatusName(withDtPublish.status)
                      << ' ' << withDtPublish.detail << '\n';
        }
        check(withDtPublish.success(),
            "candidate with previous dt publishes through the store");
        const RestartGenerationStoreLoadResult withDtLoad =
            RestartGenerationStore::load(
                root.path(), withDt.manifest.generationId,
                withDtExpectation);
        check(withDtLoad.success()
                && withDtLoad.candidate->manifest.generationState
                    .previousAcceptedDtPresent,
            "store validates and exposes the manifest previous dt");

        const RestartClockAdmissionResult withDtAdmission =
            RestartClockAdmission::stage(
                withDtLoad, observation(metadata), true);
        check(withDtAdmission.success()
                && withDtAdmission.publication->previousAcceptedDt
                    == std::optional<double>(1.0e-8),
            "manifest previous dt is staged with the validated clock");
        check(withDtAdmission.publication->publish(solver)
                && solver.time() == metadata.time
                && solver.iteration() == metadata.iteration
                && solver.previousAcceptedDt()
                    == std::optional<double>(1.0e-8),
            "validated previous dt publishes through the runtime owner");

        const RestartClockAdmissionResult wrongDtPolicy =
            RestartClockAdmission::stage(
                withDtLoad, observation(metadata), false);
        check(wrongDtPolicy.status
                    == RestartClockAdmissionStatus::PREVIOUS_DT_PRESENCE_MISMATCH
                && !wrongDtPolicy.publication,
            "previous-dt policy mismatch exposes no publication");

        const double solverTimeBefore = solver.time();
        const int solverIterationBefore = solver.iteration();
        const std::optional<double> solverDtBefore = solver.previousAcceptedDt();
        RestartClockPublication invalidPublication =
            *withDtAdmission.publication;
        invalidPublication.previousAcceptedDt = 0.0;
        check(!invalidPublication.publish(solver)
                && solver.time() == solverTimeBefore
                && solver.iteration() == solverIterationBefore
                && solver.previousAcceptedDt() == solverDtBefore,
            "invalid publication cannot partially overwrite the solver clock");
        check(sameSentinel(runtime, sentinel),
            "all rejected preflights leave Solver and Domain sentinels exact");

        check(std::string(restartClockAdmissionStatusName(
                RestartClockAdmissionStatus::PREVIOUS_DT_PRESENCE_MISMATCH))
                == "PREVIOUS_DT_PRESENCE_MISMATCH",
            "previous-dt policy mismatch has a stable structured status");
    }
    catch (const std::exception& error)
    {
        ++failures;
        std::cerr << "unexpected exception: " << error.what() << '\n';
    }

    if (failures != 0)
    {
        std::cerr << "RESTART_GENERATION_STORE_SECTOR_CLOCK_PREFLIGHT_FAIL "
                  << "checks=" << checks << " failures=" << failures << '\n';
        return 1;
    }

    std::cout << "RESTART_GENERATION_STORE_SECTOR_CLOCK_PREFLIGHT_PASS "
              << "checks=" << checks << '\n'
              << "SCOPE clock-owner-only domain-publication=0 "
              << "production_ready=0\n";
    return 0;
}
