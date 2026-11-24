#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

#include "CompleteCellState.H"
#include "ModelCapabilities.H"
#include "PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "PhaseEquilibrium/PRMixtureKernel.H"
#include "PhaseEquilibrium/PRPhaseCaloric.H"
#include "PhaseEquilibrium/PRPhaseDisplacementAuthority.H"
#include "PhaseEquilibrium/PRPhaseRhoERecovery.H"
#include "PhaseEquilibrium/PRPhaseThermoOwnerAdapter.H"
#include "PhaseEquilibrium/PRTPFlash.H"
#include "StateContract.H"

namespace
{

using Foam::BinaryInteractionDB;
using Foam::CompleteCellConversionResult;
using Foam::CompleteCellState;
using Foam::CompleteStateConversion;
using Foam::ConservativeStateContract;
using Foam::ModelCapabilities;
using Foam::PengRobinsonEOS;
using Foam::PrimitiveVars;
using Foam::SpeciesCriticalDB;
using Foam::phaseEquilibrium::EquilibriumRegion;
using Foam::phaseEquilibrium::PRMixtureKernel;
using Foam::phaseEquilibrium::PRPhaseCaloric;
using Foam::phaseEquilibrium::PRPhaseCaloricProperties;
using Foam::phaseEquilibrium::PRPhaseDisplacementAuthority;
using Foam::phaseEquilibrium::PRPhaseDisplacementBounds;
using Foam::phaseEquilibrium::PRPhaseDisplacementStatus;
using Foam::phaseEquilibrium::PRPhaseIdentity;
using Foam::phaseEquilibrium::PRPhaseRhoEOptions;
using Foam::phaseEquilibrium::PRPhaseRhoERecovery;
using Foam::phaseEquilibrium::PRPhaseRootIdentity;
using Foam::phaseEquilibrium::PRPhaseThermoLiveToken;
using Foam::phaseEquilibrium::PRPhaseThermoOwnerAdapter;
using Foam::phaseEquilibrium::PRPhaseThermoOwnerRequest;
using Foam::phaseEquilibrium::PRPhaseThermoOwnerResult;
using Foam::phaseEquilibrium::PRPhaseThermoOwnerStatus;
using Foam::phaseEquilibrium::PRPhaseThermoProbeResult;
using Foam::phaseEquilibrium::PRPhaseThermoParentAdmission;
using Foam::phaseEquilibrium::PRPhaseThermoParentView;
using Foam::phaseEquilibrium::PRTPFlash;
using Foam::phaseEquilibrium::RhoEFlashOptions;
using Foam::phaseEquilibrium::TPFlashResult;
using Foam::thermo::EquilibriumResult;
using Foam::thermo::EquilibriumState;
using Foam::thermo::PREquilibriumThermoModel;

static_assert(!ModelCapabilities::phaseAwareChemistry,
    "token-bound displacement must not enable phase-aware chemistry");
static_assert(!std::is_copy_constructible<PRPhaseThermoOwnerAdapter>::value,
    "a live displacement token authority must remain unique");

int checks = 0;
int failures = 0;

void check(bool condition, const std::string& label)
{
    ++checks;
    std::cout << (condition ? "  PASS " : "  FAIL ") << label << '\n';
    if (!condition) ++failures;
}

bool noPublication(const PRPhaseThermoOwnerResult& result)
{
    return !result.phaseState && !result.derivatives && !result.acceptedToken;
}

bool close(double lhs, double rhs, double relative = 1.0e-9)
{
    return std::isfinite(lhs) && std::isfinite(rhs)
        && std::abs(lhs - rhs) <= relative
            *std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

RhoEFlashOptions parentOptions()
{
    RhoEFlashOptions options;
    options.minimumPressure = 1.0e5;
    options.maximumPressure = 2.0e7;
    options.minimumTemperature = 160.0;
    options.maximumTemperature = 1500.0;
    options.initialPressure = 2.0e6;
    options.initialTemperature = 300.0;
    options.pressureScanIntervals = 64;
    options.maximumIterations = 120;
    options.finiteDifferenceStep = 1.0e-4;
    options.logStepLimit = 0.75;
    options.densityAbsoluteTolerance = 1.0e-7;
    options.densityRelativeTolerance = 2.0e-8;
    options.energyAbsoluteTolerance = 1.0e-4;
    options.energyRelativeTolerance = 2.0e-10;
    return options;
}

PRPhaseRhoEOptions phaseOptions()
{
    PRPhaseRhoEOptions options;
    options.minimumTemperature = 250.0;
    options.maximumTemperature = 500.0;
    options.temperatureScanIntervals = 128;
    options.maximumIterations = 80;
    options.derivativeRelativeStep = 3.0e-8;
    options.energyAbsoluteTolerance = 1.0e-5;
    options.energyRelativeTolerance = 2.0e-10;
    options.densityRelativeTolerance = 2.0e-9;
    options.rootRelativeTolerance = 2.0e-12;
    options.derivativeIdentityRelativeTolerance = 2.0e-6;
    options.derivativeIdentityAbsoluteTolerance = 1.0e-4;
    return options;
}

PRPhaseDisplacementBounds displacementBounds()
{
    return PRPhaseDisplacementBounds(
        0.02,
        0.05,
        1.0e-3,
        3,
        1.0e-5,
        1.0e-4,
        0.05,
        0.05,
        0.05,
        2.0e-9,
        1.0e-5,
        2.0e-10,
        2.0e-12);
}

PRPhaseDisplacementBounds narrowerDisplacementBounds()
{
    return PRPhaseDisplacementBounds(
        0.01,
        0.025,
        5.0e-4,
        4,
        2.0e-5,
        2.0e-4,
        0.025,
        0.025,
        0.025,
        1.0e-9,
        5.0e-6,
        1.0e-10,
        1.0e-12);
}

std::optional<CompleteCellState> completeState(
    const EquilibriumState& equilibrium,
    const PREquilibriumThermoModel& model)
{
    PrimitiveVars primitive(
        static_cast<int>(equilibrium.overallMassFractions.size()));
    primitive.p = equilibrium.p;
    primitive.T = equilibrium.T;
    primitive.U = {3.0, -2.0, 1.0};
    primitive.Y = equilibrium.overallMassFractions;
    primitive.Z = 0.25;
    primitive.C = 0.5;
    const CompleteCellConversionResult converted =
        CompleteStateConversion::fromPrimitiveAndEquilibrium(
            primitive, equilibrium, model);
    return converted.success() ? converted.state : std::nullopt;
}

std::optional<PRPhaseCaloricProperties> fixedDensityVaporState(
    const PRMixtureKernel& kernel,
    const PRPhaseCaloric& caloric,
    const PengRobinsonEOS& eos,
    double density,
    double temperature,
    const std::vector<double>& massFractions)
{
    std::vector<double> moleFractions;
    if (kernel.massToMoleFractions(massFractions, moleFractions)
        != Foam::phaseEquilibrium::PhaseEquilibriumStatus::SUCCESS)
    {
        return std::nullopt;
    }
    const double pressure = eos.pressure(
        density, temperature, massFractions);
    if (!(std::isfinite(pressure) && pressure > 0.0)) return std::nullopt;
    const auto evaluated = caloric.evaluate(
        pressure,
        temperature,
        moleFractions,
        Foam::phaseEquilibrium::PhaseRoot::VAPOR);
    if (!evaluated.success() || !close(evaluated.state->rho, density, 2.0e-9))
        return std::nullopt;
    return *evaluated.state;
}

PRPhaseThermoOwnerRequest requestFor(
    const PRPhaseCaloricProperties& phase,
    double pressureSeed)
{
    PRPhaseThermoOwnerRequest request;
    request.phaseDensity = phase.rho;
    request.phaseSpecificInternalEnergy = phase.e;
    request.phaseMassFractions = phase.massFractions;
    request.pressureSeed = pressureSeed;
    request.temperatureSeed = phase.temperature;
    return request;
}

bool parentUnchanged(
    const CompleteCellState& before,
    const CompleteCellState& after)
{
    return before.U.rho == after.U.rho
        && before.U.rhoE == after.U.rhoE
        && before.U.rhoY == after.U.rhoY
        && before.Q.p == after.Q.p
        && before.Q.T == after.Q.T
        && before.Q.rho == after.Q.rho
        && before.Q.e == after.Q.e
        && before.Q.region == after.Q.region
        && before.Q.vapor.has_value() == after.Q.vapor.has_value()
        && (!before.Q.vapor
            || (before.Q.vapor->rho == after.Q.vapor->rho
                && before.Q.vapor->e == after.Q.vapor->e
                && before.Q.vapor->massFractions
                    == after.Q.vapor->massFractions));
}

} // End anonymous namespace

int main()
{
    std::cout << "=== Token-bound PR phase displacement owner ===\n";

    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB kij = BinaryInteractionDB::createDefault();
    const std::vector<std::string> species = {"O2", "C3H8"};

    PRTPFlash flash(critical, kij, species);
    const TPFlashResult feed = flash.flashMoleFractions(
        2.0e6, 300.0, {0.2, 0.8});
    check(feed.valid() && feed.region == EquilibriumRegion::TWO_PHASE,
        "fixture supplies a two-phase parent");
    if (!feed.valid() || feed.region != EquilibriumRegion::TWO_PHASE) return 1;

    const RhoEFlashOptions mixtureOptions = parentOptions();
    PREquilibriumThermoModel model(
        critical, kij, species, mixtureOptions);
    const EquilibriumResult initialized = model.initializePTY(
        feed.pressure, feed.temperature, feed.feedMassFractions);
    check(initialized.success() && initialized.state
            && initialized.state->vapor,
        "whole-mixture owner publishes a vapor parent payload");
    if (!initialized.success() || !initialized.state
        || !initialized.state->vapor)
    {
        return 1;
    }

    const std::optional<CompleteCellState> committed = completeState(
        *initialized.state, model);
    check(committed && committed->valid(),
        "fixture closes a complete parent U-W-Q state");
    if (!committed) return 1;
    const CompleteCellState immutableParent = *committed;

    ConservativeStateContract stateContract;
    stateContract.setSpeciesLayout(species, species.size());
    const PRPhaseRhoEOptions localOptions = phaseOptions();
    PRPhaseRhoERecovery recovery(critical, kij, species, localOptions);
    PRPhaseDisplacementAuthority displacement(
        critical, kij, species, localOptions, displacementBounds());
    check(recovery.configured() && displacement.configured()
            && recovery.authorityFingerprint()
                == displacement.recoveryAuthorityFingerprint()
            && recovery.physicalModelAuthorityFingerprint()
                == displacement.physicalModelAuthorityFingerprint()
            && recovery.physicalModelAuthorityFingerprint()
                == model.phaseLocalPhysicalModelFingerprint(),
        "recovery, displacement, and model share one frozen authority");

    const PRPhaseThermoParentAdmission admission{
        "box-c/cell-17", 7U, 11U,
        PRPhaseIdentity::VAPOR,
        PRPhaseRootIdentity::PR_VAPOR_ROOT};
    const PRPhaseThermoParentView parent(
        stateContract, *committed, model, admission);

    PRMixtureKernel kernel(critical, kij, species);
    PRPhaseCaloric caloric(critical, kij, species);
    PengRobinsonEOS eos(critical, kij);
    eos.setSpeciesNames(species);
    std::vector<double> displacedY = committed->Q.vapor->massFractions;
    displacedY[0] += 1.0e-4;
    displacedY[1] -= 1.0e-4;
    const auto displacedState = fixedDensityVaporState(
        kernel,
        caloric,
        eos,
        committed->Q.vapor->rho,
        committed->Q.T + 1.0,
        displacedY);
    check(displacedState.has_value(),
        "nearby rho-e-Y request is independently materialized");
    if (!displacedState) return 1;

    PRPhaseThermoOwnerAdapter adapter(recovery, displacement, model);
    const PRPhaseThermoOwnerResult opened = adapter.open(parent, 310.0);
    check(opened.success() && opened.acceptedToken
            && opened.acceptedToken->transactionSequence == 1U
            && !opened.diagnostics.displacementAuthorityInvoked,
        "parent replay opens sequence one before displacement");
    if (!opened.acceptedToken) return 1;

    PRPhaseDisplacementAuthority narrowerDisplacement(
        critical, kij, species, localOptions, narrowerDisplacementBounds());
    PRPhaseThermoOwnerAdapter narrowerAdapter(
        recovery, narrowerDisplacement, model);
    const PRPhaseThermoOwnerResult narrowerOpen = narrowerAdapter.open(
        parent, 310.0);
    check(narrowerOpen.success() && narrowerOpen.acceptedToken
            && displacement.authorityFingerprint()
                != narrowerDisplacement.authorityFingerprint()
            && opened.diagnostics.parentBindingDigest
                != narrowerOpen.diagnostics.parentBindingDigest
            && opened.acceptedToken->transactionId
                != narrowerOpen.acceptedToken->transactionId
            && opened.acceptedToken->tokenDigest
                != narrowerOpen.acceptedToken->tokenDigest,
        "displacement bounds are part of parent, transaction, and token identity");
    narrowerAdapter.abort();

    const PRPhaseThermoOwnerRequest displacedRequest = requestFor(
        *displacedState, committed->Q.p);
    const PRPhaseThermoOwnerRequest parentReplayRequest = requestFor(
        *opened.phaseState, committed->Q.p);
    const PRPhaseThermoProbeResult firstProbe = adapter.probe(
        parent, displacedRequest, *opened.acceptedToken);
    const PRPhaseThermoProbeResult replayProbe = adapter.probe(
        parent, parentReplayRequest, *opened.acceptedToken);
    const PRPhaseThermoProbeResult repeatedProbe = adapter.probe(
        parent, displacedRequest, *opened.acceptedToken);
    PRPhaseThermoOwnerRequest roundoffProbeRequest = parentReplayRequest;
    roundoffProbeRequest.phaseSpecificInternalEnergy = std::nextafter(
        roundoffProbeRequest.phaseSpecificInternalEnergy,
        std::numeric_limits<double>::infinity());
    roundoffProbeRequest.phaseMassFractions[0] = std::nextafter(
        roundoffProbeRequest.phaseMassFractions[0],
        std::numeric_limits<double>::infinity());
    const PRPhaseThermoProbeResult roundoffProbe = adapter.probe(
        parent, roundoffProbeRequest, *opened.acceptedToken);
    check(firstProbe.success() && replayProbe.success()
            && repeatedProbe.success()
            && firstProbe.phaseState && repeatedProbe.phaseState
            && close(
                firstProbe.phaseState->rho,
                repeatedProbe.phaseState->rho, 0.0)
            && close(
                firstProbe.phaseState->e,
                repeatedProbe.phaseState->e, 0.0)
            && firstProbe.phaseState->massFractions
                == repeatedProbe.phaseState->massFractions
            && !firstProbe.diagnostics.resultTokenIssued
            && !replayProbe.diagnostics.resultTokenIssued
            && !repeatedProbe.diagnostics.resultTokenIssued
            && adapter.currentSequence() == 1U
            && adapter.transactionActive(),
        "repeated and out-of-order probes are pure and sequence independent");
    check(roundoffProbe.success()
            && roundoffProbe.diagnostics
                .requestCanonicalizedToAdmittedState
            && !roundoffProbe.diagnostics.resultTokenIssued
            && adapter.currentSequence() == 1U
            && adapter.transactionActive(),
        "sub-closure numerical noise probes canonicalize without token mutation");

    PRPhaseThermoOwnerRequest rejectedProbeRequest = displacedRequest;
    rejectedProbeRequest.phaseSpecificInternalEnergy += 1.0e9;
    const PRPhaseThermoProbeResult rejectedProbe = adapter.probe(
        parent, rejectedProbeRequest, *opened.acceptedToken);
    check(rejectedProbe.status
                == PRPhaseThermoOwnerStatus::DISPLACEMENT_AUTHORITY_FAILED
            && !rejectedProbe.phaseState && !rejectedProbe.derivatives
            && !rejectedProbe.diagnostics.resultTokenIssued
            && adapter.currentSequence() == 1U
            && adapter.transactionActive(),
        "failed numerical probe publishes nothing and preserves the live token");

    const PRPhaseThermoOwnerResult moved = adapter.recover(
        parent,
        displacedRequest,
        *opened.acceptedToken);
    check(moved.success() && moved.acceptedToken && moved.phaseState
            && moved.acceptedToken->transactionSequence == 2U
            && moved.acceptedToken->rootIdentity
                == PRPhaseRootIdentity::PR_VAPOR_ROOT
            && moved.diagnostics.inputTokenAdmitted
            && moved.diagnostics.displacementAuthorityInvoked
            && moved.diagnostics.displacementStatus
                == PRPhaseDisplacementStatus::SUCCESS,
        "current token authorizes bounded displacement on the parent-derived root");
    check(moved.phaseState
            && close(moved.phaseState->rho, displacedState->rho, 2.0e-9)
            && close(moved.phaseState->e, displacedState->e, 2.0e-10)
            && moved.phaseState->massFractions == displacedState->massFractions,
        "successful publication closes against requested rho-e-Y");
    check(parentUnchanged(immutableParent, *committed),
        "phase-local success cannot mutate the immutable outer parent");
    if (!moved.acceptedToken || !moved.phaseState) return 1;

    const PRPhaseThermoOwnerResult replayMoved = adapter.recover(
        parent,
        requestFor(*moved.phaseState, moved.phaseState->pressure),
        *moved.acceptedToken);
    check(replayMoved.success() && replayMoved.acceptedToken
            && replayMoved.acceptedToken->transactionSequence == 3U
            && replayMoved.diagnostics.displacementAuthorityInvoked,
        "exact replay of displaced published bytes replaces the live token");
    if (!replayMoved.acceptedToken) return 1;

    const PRPhaseThermoOwnerResult stale = adapter.recover(
        parent,
        requestFor(*moved.phaseState, moved.phaseState->pressure),
        *moved.acceptedToken);
    check(stale.status == PRPhaseThermoOwnerStatus::STALE_TOKEN
            && noPublication(stale)
            && !stale.diagnostics.displacementAuthorityInvoked
            && !adapter.transactionActive(),
        "superseded token closes before displacement work");

    PRPhaseThermoOwnerAdapter driftAdapter(recovery, displacement, model);
    const PRPhaseThermoOwnerResult driftOpen = driftAdapter.open(parent, 310.0);
    check(driftOpen.success() && driftOpen.acceptedToken,
        "generation-drift regression opens a live transaction");
    if (!driftOpen.acceptedToken) return 1;
    PRPhaseThermoParentAdmission nextGeneration = admission;
    ++nextGeneration.parentGeneration;
    const PRPhaseThermoParentView driftedParent(
        stateContract, *committed, model, nextGeneration);
    const PRPhaseThermoOwnerResult generationDrift = driftAdapter.recover(
        driftedParent,
        requestFor(*displacedState, committed->Q.p),
        *driftOpen.acceptedToken);
    check(generationDrift.status
                == PRPhaseThermoOwnerStatus::PARENT_BINDING_DRIFT
            && noPublication(generationDrift)
            && !generationDrift.diagnostics.displacementAuthorityInvoked
            && !driftAdapter.transactionActive(),
        "generation drift closes before token-bound displacement");

    const EquilibriumResult driftedInitialization = model.initializePTY(
        feed.pressure, feed.temperature + 1.0, feed.feedMassFractions);
    const std::optional<CompleteCellState> driftedCommitted =
        driftedInitialization.success() && driftedInitialization.state
        ? completeState(*driftedInitialization.state, model)
        : std::nullopt;
    check(driftedCommitted && driftedCommitted->valid(),
        "parent-byte drift fixture remains a complete model-valid state");
    if (!driftedCommitted) return 1;
    PRPhaseThermoOwnerAdapter byteDriftAdapter(
        recovery, displacement, model);
    const PRPhaseThermoOwnerResult byteDriftOpen = byteDriftAdapter.open(
        parent, 310.0);
    check(byteDriftOpen.success() && byteDriftOpen.acceptedToken,
        "parent-byte drift regression opens a live transaction");
    if (!byteDriftOpen.acceptedToken) return 1;
    const PRPhaseThermoParentView byteDriftedParent(
        stateContract, *driftedCommitted, model, admission);
    const PRPhaseThermoOwnerResult byteDriftRejected =
        byteDriftAdapter.recover(
            byteDriftedParent,
            requestFor(*displacedState, committed->Q.p),
            *byteDriftOpen.acceptedToken);
    check(byteDriftRejected.status
                == PRPhaseThermoOwnerStatus::PARENT_BINDING_DRIFT
            && noPublication(byteDriftRejected)
            && !byteDriftRejected.diagnostics.displacementAuthorityInvoked
            && !byteDriftAdapter.transactionActive(),
        "valid changed parent bytes close before token-bound displacement");

    PREquilibriumThermoModel modelClone(
        critical, kij, species, mixtureOptions);
    PRPhaseThermoOwnerAdapter cloneAdapter(recovery, displacement, model);
    const PRPhaseThermoParentView clonedModelParent(
        stateContract, *committed, modelClone, admission);
    const PRPhaseThermoOwnerResult cloneRejected = cloneAdapter.open(
        clonedModelParent, 310.0);
    check(cloneRejected.status
                == PRPhaseThermoOwnerStatus::INVALID_PARENT_ADMISSION
            && noPublication(cloneRejected),
        "same-configuration model clone cannot replace the bound model instance");

    PRPhaseRhoEOptions mismatchedOptions = localOptions;
    mismatchedOptions.maximumIterations += 1;
    PRPhaseDisplacementAuthority mismatchedDisplacement(
        critical,
        kij,
        species,
        mismatchedOptions,
        displacementBounds());
    PRPhaseThermoOwnerAdapter mismatchedAdapter(
        recovery, mismatchedDisplacement, model);
    const PRPhaseThermoOwnerResult mismatchRejected =
        mismatchedAdapter.open(parent, 310.0);
    check(mismatchRejected.status
                == PRPhaseThermoOwnerStatus::INVALID_CONFIGURATION
            && noPublication(mismatchRejected),
        "recovery-option rebind fails before parent or token work");

    BinaryInteractionDB shiftedKij = kij;
    shiftedKij.add("O2", "C3H8", 0.04);
    PRPhaseDisplacementAuthority shiftedDisplacement(
        critical,
        shiftedKij,
        species,
        localOptions,
        displacementBounds());
    PRPhaseThermoOwnerAdapter shiftedAdapter(
        recovery, shiftedDisplacement, model);
    const PRPhaseThermoOwnerResult shiftedRejected = shiftedAdapter.open(
        parent, 310.0);
    check(shiftedRejected.status
                == PRPhaseThermoOwnerStatus::INVALID_CONFIGURATION
            && noPublication(shiftedRejected),
        "critical-kij-caloric physical authority rebind fails closed");

    PREquilibriumThermoModel shiftedModel(
        critical, shiftedKij, species, mixtureOptions);
    PRPhaseThermoOwnerAdapter shiftedModelAdapter(
        recovery, displacement, shiftedModel);
    const PRPhaseThermoOwnerResult shiftedModelRejected =
        shiftedModelAdapter.open(parent, 310.0);
    check(shiftedModelRejected.status
                == PRPhaseThermoOwnerStatus::INVALID_CONFIGURATION
            && noPublication(shiftedModelRejected),
        "whole-mixture model physical authority must match phase kernels");

    PRPhaseThermoOwnerAdapter boundsAdapter(recovery, displacement, model);
    const PRPhaseThermoOwnerResult boundsOpen = boundsAdapter.open(parent, 310.0);
    check(boundsOpen.success() && boundsOpen.acceptedToken,
        "rollback regression opens a live transaction");
    if (!boundsOpen.acceptedToken) return 1;
    PRPhaseThermoOwnerRequest outOfBounds = requestFor(
        *displacedState, committed->Q.p);
    outOfBounds.phaseSpecificInternalEnergy += 1.0e9;
    const PRPhaseThermoOwnerResult rejected = boundsAdapter.recover(
        parent, outOfBounds, *boundsOpen.acceptedToken);
    check(rejected.status
                == PRPhaseThermoOwnerStatus::DISPLACEMENT_AUTHORITY_FAILED
            && rejected.diagnostics.displacementAuthorityInvoked
            && rejected.diagnostics.displacementStatus
                == PRPhaseDisplacementStatus::DISPLACEMENT_OUT_OF_BOUNDS
            && noPublication(rejected)
            && !boundsAdapter.transactionActive()
            && parentUnchanged(immutableParent, *committed),
        "failed displacement publishes nothing, closes token, and rolls back locally");

    PRPhaseThermoParentAdmission forgedAdmission = admission;
    forgedAdmission.phaseIdentity = PRPhaseIdentity::SUPERCRITICAL;
    forgedAdmission.rootIdentity =
        PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT;
    const PRPhaseThermoParentView forgedParent(
        stateContract, *committed, model, forgedAdmission);
    PRPhaseThermoOwnerAdapter forgedAdapter(recovery, displacement, model);
    const PRPhaseThermoOwnerResult forgedRejected = forgedAdapter.open(
        forgedParent, 310.0);
    check(forgedRejected.status
                == PRPhaseThermoOwnerStatus::INVALID_PARENT_PHASE_ADMISSION
            && noPublication(forgedRejected),
        "raw phase/root assertions cannot override the validated parent region");

    if (failures == 0)
    {
        std::cout
            << "PR_PHASE_THERMO_OWNER_DISPLACEMENT_PASS checks=" << checks
            << " token_bound=1 parent_lineage=1 exact_replay=1 rollback=1"
            << " model_authority=1 raw_identity_input=0 capability=0"
            << " production_ready=0\n";
        return 0;
    }
    std::cout << "PR_PHASE_THERMO_OWNER_DISPLACEMENT_FAIL checks=" << checks
              << " failures=" << failures << '\n';
    return 1;
}
