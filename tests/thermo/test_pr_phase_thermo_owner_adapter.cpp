#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

#include "CompleteCellState.H"
#include "MixingRules/BinaryInteractionDB.H"
#include "MixingRules/SpeciesCriticalDB.H"
#include "PhaseEquilibrium/PREquilibriumThermoModel.H"
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
using Foam::PrimitiveVars;
using Foam::SpeciesCriticalDB;
using Foam::phaseEquilibrium::PRPhaseIdentity;
using Foam::phaseEquilibrium::PRPhaseRhoEOptions;
using Foam::phaseEquilibrium::PRPhaseRhoERecovery;
using Foam::phaseEquilibrium::PRPhaseRhoEStatus;
using Foam::phaseEquilibrium::PRPhaseRootIdentity;
using Foam::phaseEquilibrium::PRPhaseThermoLiveToken;
using Foam::phaseEquilibrium::PRPhaseThermoOwnerAdapter;
using Foam::phaseEquilibrium::PRPhaseThermoOwnerRequest;
using Foam::phaseEquilibrium::PRPhaseThermoOwnerResult;
using Foam::phaseEquilibrium::PRPhaseThermoOwnerStatus;
using Foam::phaseEquilibrium::PRPhaseThermoParentAdmission;
using Foam::phaseEquilibrium::PRPhaseThermoParentView;
using Foam::phaseEquilibrium::PRTPFlash;
using Foam::phaseEquilibrium::RhoEFlashOptions;
using Foam::phaseEquilibrium::TPFlashResult;
using Foam::thermo::EquilibriumRegion;
using Foam::thermo::EquilibriumResult;
using Foam::thermo::EquilibriumState;
using Foam::thermo::PREquilibriumThermoModel;

static_assert(!std::is_copy_constructible<PRPhaseThermoOwnerAdapter>::value,
    "a phase owner must not duplicate live token authority");
static_assert(!std::is_move_constructible<PRPhaseThermoOwnerAdapter>::value,
    "a phase owner must not leave moved-from token authority");
static_assert(!std::is_constructible<
        PRPhaseThermoOwnerAdapter, PRPhaseRhoERecovery&&>::value,
    "the non-owning adapter must reject a temporary recovery kernel");

int checks = 0;
int failures = 0;

void check(bool condition, const std::string& label)
{
    ++checks;
    std::cout << (condition ? "  PASS " : "  FAIL ") << label << '\n';
    if (!condition) ++failures;
}

bool close(double lhs, double rhs, double relative = 1.0e-8)
{
    return std::isfinite(lhs) && std::isfinite(rhs)
        && std::abs(lhs - rhs) <= relative
            *std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

bool digest(const std::string& value)
{
    if (value.size() != 64U || value == std::string(64U, '0')) return false;
    for (char character : value)
    {
        if (!((character >= '0' && character <= '9')
            || (character >= 'a' && character <= 'f')))
            return false;
    }
    return true;
}

RhoEFlashOptions parentOptions(
    double minimumTemperature = 160.0,
    double maximumTemperature = 1500.0)
{
    RhoEFlashOptions options;
    options.minimumPressure = 1.0e5;
    options.maximumPressure = 2.0e7;
    options.minimumTemperature = minimumTemperature;
    options.maximumTemperature = maximumTemperature;
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

PRPhaseRhoEOptions phaseOptions(
    double minimumTemperature = 250.0,
    double maximumTemperature = 500.0)
{
    PRPhaseRhoEOptions options;
    options.minimumTemperature = minimumTemperature;
    options.maximumTemperature = maximumTemperature;
    options.temperatureScanIntervals = 128;
    options.maximumIterations = 80;
    options.derivativeRelativeStep = 3.0e-8;
    return options;
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

PRPhaseThermoOwnerRequest requestFor(
    const Foam::thermo::EquilibriumPhaseProperties& phase,
    double pressureSeed,
    double temperatureSeed)
{
    PRPhaseThermoOwnerRequest request;
    request.phaseDensity = phase.rho;
    request.phaseSpecificInternalEnergy = phase.e;
    request.phaseMassFractions = phase.massFractions;
    request.pressureSeed = pressureSeed;
    request.temperatureSeed = temperatureSeed;
    return request;
}

bool noPublication(const PRPhaseThermoOwnerResult& result)
{
    return !result.phaseState.has_value()
        && !result.derivatives.has_value()
        && !result.acceptedToken.has_value();
}

} // namespace

int main()
{
    std::cout << "=== Disconnected token-bound PR phase owner adapter ===\n";

    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB kij = BinaryInteractionDB::createDefault();
    const std::vector<std::string> species = {"O2", "C3H8"};

    PRTPFlash feedFlash(critical, kij, species);
    const TPFlashResult feed = feedFlash.flashMoleFractions(
        2.0e6, 300.0, {0.2, 0.8});
    check(feed.valid() && feed.region
            == Foam::phaseEquilibrium::EquilibriumRegion::TWO_PHASE,
        "fixture supplies a two-phase O2-C3H8 parent");
    if (!feed.valid()
        || feed.region != Foam::phaseEquilibrium::EquilibriumRegion::TWO_PHASE)
        return 1;

    const RhoEFlashOptions mixtureOptions = parentOptions();
    PREquilibriumThermoModel model(
        critical, kij, species, mixtureOptions);
    const EquilibriumResult initialized = model.initializePTY(
        2.0e6, 300.0, feed.feedMassFractions);
    check(initialized.success() && initialized.state
            && initialized.state->region == EquilibriumRegion::TWO_PHASE
            && initialized.state->vapor.has_value(),
        "parent model publishes a complete model-valid two-phase state");
    if (!initialized.success() || !initialized.state
        || !initialized.state->vapor)
        return 1;

    const std::optional<CompleteCellState> committed = completeState(
        *initialized.state, model);
    check(committed.has_value() && committed->valid()
            && model.stateValidForModel(committed->Q),
        "fixture closes the full U-W-Q parent state");
    if (!committed) return 1;

    ConservativeStateContract stateContract;
    stateContract.setSpeciesLayout(species, species.size());
    const PRPhaseRhoEOptions localOptions = phaseOptions();
    PRPhaseRhoERecovery recovery(
        critical, kij, species, localOptions);
    check(recovery.configured(), "phase-local recovery kernel configures");

    const PRPhaseThermoParentAdmission admission{
        "box-c/cell-17", 7U, 11U,
        PRPhaseIdentity::VAPOR,
        PRPhaseRootIdentity::PR_VAPOR_ROOT};
    const PRPhaseThermoParentView parent(
        stateContract, *committed, model, admission);

    PRPhaseThermoOwnerAdapter adapter(recovery);
    const PRPhaseThermoOwnerResult opened = adapter.open(parent, 310.0);
    std::cout << "    vapor_open_status="
              << Foam::phaseEquilibrium::phaseThermoOwnerStatusName(opened.status)
              << " kernel="
              << Foam::phaseEquilibrium::phaseRhoEStatusName(
                    opened.diagnostics.kernelStatus)
              << " residual_ok=" << opened.residuals.accepted()
              << " residuals=" << opened.residuals.pressureRelative << ','
              << opened.residuals.temperatureRelative << ','
              << opened.residuals.compressibilityRelative << ','
              << opened.residuals.densityRelative << ','
              << opened.residuals.specificInternalEnergyRelative << ','
              << opened.residuals.specificEnthalpyRelative << ','
              << opened.residuals.massFractionL1 << ','
              << opened.residuals.eosRoot << ','
              << opened.residuals.departureCaloric << ','
              << opened.residuals.cvMinusDedTRho << ','
              << opened.residuals.rootContinuity << '\n';
    check(opened.success() && opened.acceptedToken
            && opened.acceptedToken->transactionSequence == 1U,
        "validated parent replay opens one sequence-one transaction");
    check(opened.diagnostics.parentStateContractValidated
            && opened.diagnostics.parentFingerprintsComputed
            && !opened.diagnostics.inputTokenAdmitted
            && opened.diagnostics.kernelInvoked
            && opened.diagnostics.resultTokenIssued,
        "open ordering records parent validation before fingerprint and kernel");
    check(digest(opened.diagnostics.committedParentQFingerprint)
            && digest(opened.diagnostics.parentPhasePayloadFingerprint)
            && digest(opened.diagnostics.parentBindingDigest)
            && digest(opened.acceptedToken->transactionId)
            && digest(opened.acceptedToken->tokenDigest),
        "successful replay emits closed nonzero SHA-256 identities");
    check(opened.diagnostics.committedParentQFingerprint
                == "179b898fc0ac3b5b400bdc9c76e311a8c6373829447db72b34883264b7b033ea"
            && opened.diagnostics.parentPhasePayloadFingerprint
                == "09b088b92c61079bd4b6b483a527b3ed417031ecc142644ca291738cb248af87"
            && opened.diagnostics.parentBindingDigest
                == "d8e9e191ae8ac59f9b087374d985f4e72a307921094bb4f1aff4e1ac4db07c93"
            && opened.acceptedToken->transactionId
                == "b83c5ba0d7bb9fe83eb78e747bda344ce1e4d90d29347501e2d9b2906ccb6f0b"
            && opened.acceptedToken->tokenDigest
                == "6e29ce1d42c0692b39b2984ab8933813a012271e63dacff40717e495c128e3f4",
        "C++ parent and token identities match the independent v4 oracle");

    PRPhaseThermoOwnerAdapter tamperedTokenAdapter(recovery);
    const PRPhaseThermoOwnerResult tamperedTokenOpen =
        tamperedTokenAdapter.open(parent, 310.0);
    check(tamperedTokenOpen.success() && tamperedTokenOpen.acceptedToken,
        "token-tampering regression starts from a current live token");
    if (!tamperedTokenOpen.acceptedToken) return 1;
    PRPhaseThermoLiveToken tamperedToken = *tamperedTokenOpen.acceptedToken;
    tamperedToken.tokenDigest[0] = tamperedToken.tokenDigest[0] == '0'
        ? '1' : '0';
    const PRPhaseThermoOwnerResult tamperedTokenRejected =
        tamperedTokenAdapter.recover(
            parent,
            requestFor(*committed->Q.vapor, committed->Q.p, committed->Q.T),
            tamperedToken);
    check(tamperedTokenRejected.status
                == PRPhaseThermoOwnerStatus::STALE_TOKEN
            && noPublication(tamperedTokenRejected)
            && !tamperedTokenRejected.diagnostics.inputTokenAdmitted
            && !tamperedTokenRejected.diagnostics.kernelInvoked
            && !tamperedTokenAdapter.transactionActive(),
        "tampered current token fails before kernel work and closes locally");

    PRPhaseThermoOwnerAdapter perParentOrdinalAdapter(recovery);
    const PRPhaseThermoOwnerResult firstParentOrdinal =
        perParentOrdinalAdapter.open(parent, 310.0);
    check(firstParentOrdinal.success(),
        "per-parent ordinal regression opens the first parent");
    perParentOrdinalAdapter.abort();
    PRPhaseThermoParentAdmission secondAdmission = admission;
    secondAdmission.parentOwnerId = "box-c/cell-18";
    const PRPhaseThermoParentView secondParent(
        stateContract, *committed, model, secondAdmission);
    const PRPhaseThermoOwnerResult reusedSecondParent =
        perParentOrdinalAdapter.open(secondParent, 310.0);
    PRPhaseThermoOwnerAdapter freshSecondParentAdapter(recovery);
    const PRPhaseThermoOwnerResult freshSecondParent =
        freshSecondParentAdapter.open(secondParent, 310.0);
    check(reusedSecondParent.success() && reusedSecondParent.acceptedToken
            && freshSecondParent.success() && freshSecondParent.acceptedToken
            && reusedSecondParent.acceptedToken->transactionId
                == freshSecondParent.acceptedToken->transactionId
            && reusedSecondParent.acceptedToken->tokenDigest
                == freshSecondParent.acceptedToken->tokenDigest,
        "transaction ordinals are independent for distinct parent owners");

    PRPhaseThermoOwnerAdapter explicitCloseAdapter(recovery);
    const PRPhaseThermoOwnerResult explicitCloseOpen =
        explicitCloseAdapter.open(parent, 310.0);
    check(explicitCloseOpen.success() && explicitCloseOpen.acceptedToken,
        "explicit-close regression starts from a current live token");
    if (!explicitCloseOpen.acceptedToken) return 1;
    explicitCloseAdapter.closeTransaction();
    const PRPhaseThermoOwnerResult afterExplicitClose =
        explicitCloseAdapter.recover(
            parent,
            requestFor(*committed->Q.vapor, committed->Q.p, committed->Q.T),
            *explicitCloseOpen.acceptedToken);
    check(afterExplicitClose.status
                == PRPhaseThermoOwnerStatus::NO_ACTIVE_TRANSACTION
            && noPublication(afterExplicitClose)
            && afterExplicitClose.diagnostics.parentStateContractValidated
            && afterExplicitClose.diagnostics.parentFingerprintsComputed
            && !afterExplicitClose.diagnostics.inputTokenAdmitted
            && !afterExplicitClose.diagnostics.kernelInvoked
            && !explicitCloseAdapter.transactionActive(),
        "explicit close invalidates the token without kernel or publication");
    check(opened.phaseState
            && close(opened.phaseState->pressure, committed->Q.p)
            && close(opened.phaseState->temperature, committed->Q.T)
            && close(opened.phaseState->rho, committed->Q.vapor->rho)
            && close(opened.phaseState->e, committed->Q.vapor->e),
        "open replays parent p-T-rho-e on the admitted vapor branch");
    check(opened.residuals.accepted(),
        "open validates every v4 replay residual before token issuance");
    check(adapter.transactionActive() && adapter.currentSequence() == 1U,
        "adapter publishes only the completed local transaction");
    if (!opened.success() || !opened.acceptedToken) return 1;

    PRPhaseRhoEOptions defaultDerivativeOptions = phaseOptions();
    defaultDerivativeOptions.derivativeRelativeStep = 2.0e-7;
    PRPhaseRhoERecovery defaultDerivativeRecovery(
        critical, kij, species, defaultDerivativeOptions);
    PRPhaseThermoOwnerAdapter defaultDerivativeAdapter(
        defaultDerivativeRecovery);
    const PRPhaseThermoOwnerResult defaultDerivativeRejected =
        defaultDerivativeAdapter.open(parent, 310.0);
    check(defaultDerivativeRejected.status
            == PRPhaseThermoOwnerStatus::RECOVERY_CONTRACT_FAILED
            && defaultDerivativeRejected.diagnostics.kernelStatus
                == PRPhaseRhoEStatus::SUCCESS
            && defaultDerivativeRejected.residuals.cvMinusDedTRho > 1.0e-8
            && noPublication(defaultDerivativeRejected)
            && !defaultDerivativeAdapter.transactionActive(),
        "v4 residual gate rejects a kernel success with coarse Cv differencing");

    const PRPhaseThermoOwnerResult duplicateOpen = adapter.open(parent, 300.0);
    check(duplicateOpen.status
            == PRPhaseThermoOwnerStatus::TRANSACTION_ALREADY_ACTIVE
            && noPublication(duplicateOpen)
            && adapter.transactionActive()
            && adapter.currentSequence() == 1U,
        "a second transaction for the same adapter fails without disturbing it");

    const PRPhaseThermoOwnerRequest replayRequest = requestFor(
        *committed->Q.vapor, committed->Q.p, 470.0);
    const PRPhaseThermoOwnerResult continued = adapter.recover(
        parent, replayRequest, *opened.acceptedToken);
    check(continued.success() && continued.acceptedToken
            && continued.acceptedToken->transactionSequence == 2U
            && continued.acceptedToken->tokenDigest
                != opened.acceptedToken->tokenDigest,
        "continued exact replay replaces the live token at sequence two");
    check(continued.diagnostics.parentStateContractValidated
            && continued.diagnostics.parentFingerprintsComputed
            && continued.diagnostics.inputTokenAdmitted
            && continued.diagnostics.kernelInvoked
            && continued.diagnostics.resultTokenIssued,
        "continued recovery follows validate-fingerprint-token-kernel order");

    const PRPhaseThermoOwnerResult stale = adapter.recover(
        parent, replayRequest, *opened.acceptedToken);
    check(stale.status == PRPhaseThermoOwnerStatus::STALE_TOKEN
            && noPublication(stale)
            && stale.diagnostics.parentStateContractValidated
            && stale.diagnostics.parentFingerprintsComputed
            && !stale.diagnostics.inputTokenAdmitted
            && !stale.diagnostics.kernelInvoked
            && !adapter.transactionActive(),
        "superseded token fails before the kernel and closes the transaction");

    const PRPhaseThermoOwnerResult reopened = adapter.open(parent, 305.0);
    check(reopened.success() && reopened.acceptedToken
            && reopened.acceptedToken->transactionSequence == 1U
            && reopened.acceptedToken->transactionId
                != opened.acceptedToken->transactionId,
        "fresh replay after abort creates a distinct sequence-one transaction");
    if (!reopened.acceptedToken) return 1;

    CompleteCellState malformedParent = *committed;
    malformedParent.U.rhoY[0] += 0.01*malformedParent.U.rho;
    PRPhaseThermoLiveToken forged = *reopened.acceptedToken;
    forged.tokenDigest.assign(64U, 'a');
    const PRPhaseThermoParentView malformedView(
        stateContract, malformedParent, model, admission);
    const PRPhaseThermoOwnerResult malformedFirst = adapter.recover(
        malformedView, replayRequest, forged);
    check(malformedFirst.status
            == PRPhaseThermoOwnerStatus::INVALID_PARENT_STATE_CONTRACT
            && noPublication(malformedFirst)
            && !malformedFirst.diagnostics.parentStateContractValidated
            && !malformedFirst.diagnostics.parentFingerprintsComputed
            && !malformedFirst.diagnostics.inputTokenAdmitted
            && !malformedFirst.diagnostics.kernelInvoked
            && !adapter.transactionActive(),
        "invalid full parent state wins over forged fingerprint and token");

    const PRPhaseThermoOwnerResult displacementOpen = adapter.open(parent, 300.0);
    check(displacementOpen.success() && displacementOpen.acceptedToken,
        "displacement rejection starts from a fresh bound transaction");
    if (!displacementOpen.acceptedToken) return 1;
    PRPhaseThermoOwnerRequest displacedRequest = replayRequest;
    displacedRequest.phaseSpecificInternalEnergy = std::nextafter(
        displacedRequest.phaseSpecificInternalEnergy,
        std::numeric_limits<double>::infinity());
    const PRPhaseThermoOwnerResult displaced = adapter.recover(
        parent, displacedRequest, *displacementOpen.acceptedToken);
    check(displaced.status
            == PRPhaseThermoOwnerStatus::DISPLACEMENT_NOT_AUTHORIZED
            && noPublication(displaced)
            && displaced.diagnostics.inputTokenAdmitted
            && !displaced.diagnostics.kernelInvoked
            && !adapter.transactionActive(),
        "unvalidated phase displacement fails closed before recovery");

    const EquilibriumResult driftedInitialization = model.initializePTY(
        2.0e6, 301.0, feed.feedMassFractions);
    const std::optional<CompleteCellState> driftedCommitted =
        driftedInitialization.success() && driftedInitialization.state
        ? completeState(*driftedInitialization.state, model) : std::nullopt;
    check(driftedCommitted.has_value()
            && driftedCommitted->Q.region == EquilibriumRegion::TWO_PHASE,
        "binding-drift fixture remains a valid two-phase parent");
    if (!driftedCommitted) return 1;
    const PRPhaseThermoOwnerResult driftOpen = adapter.open(parent, 300.0);
    check(driftOpen.success() && driftOpen.acceptedToken,
        "binding drift starts from a fresh bound transaction");
    if (!driftOpen.acceptedToken) return 1;
    const PRPhaseThermoParentView driftedView(
        stateContract, *driftedCommitted, model, admission);
    const PRPhaseThermoOwnerRequest driftedRequest = requestFor(
        *driftedCommitted->Q.vapor,
        driftedCommitted->Q.p,
        driftedCommitted->Q.T);
    const PRPhaseThermoOwnerResult drifted = adapter.recover(
        driftedView, driftedRequest, *driftOpen.acceptedToken);
    check(drifted.status == PRPhaseThermoOwnerStatus::PARENT_BINDING_DRIFT
            && noPublication(drifted)
            && drifted.diagnostics.parentStateContractValidated
            && drifted.diagnostics.parentFingerprintsComputed
            && !drifted.diagnostics.inputTokenAdmitted
            && !drifted.diagnostics.kernelInvoked
            && !adapter.transactionActive(),
        "valid-but-different parent bytes invalidate the bound transaction");

    const std::vector<std::string> propaneSpecies = {"C3H8"};
    PREquilibriumThermoModel propaneModel(
        critical, kij, propaneSpecies, parentOptions());
    const EquilibriumResult liquidInitialization = propaneModel.initializePTY(
        5.0e6, 250.0, {1.0});
    const std::optional<CompleteCellState> liquidCommitted =
        liquidInitialization.success() && liquidInitialization.state
        ? completeState(*liquidInitialization.state, propaneModel)
        : std::nullopt;
    check(liquidCommitted.has_value()
            && liquidCommitted->Q.region == EquilibriumRegion::SINGLE_LIQUID,
        "forged-identity fixture is a parent-admitted liquid state");
    if (!liquidCommitted) return 1;
    ConservativeStateContract propaneContract;
    propaneContract.setSpeciesLayout(propaneSpecies, propaneSpecies.size());
    PRPhaseRhoERecovery propaneRecovery(
        critical, kij, propaneSpecies, phaseOptions(240.0, 500.0));
    PRPhaseThermoOwnerAdapter propaneAdapter(propaneRecovery);
    const PRPhaseThermoParentAdmission forgedSupercriticalAdmission{
        "cell-liquid", 7U, 3U,
        PRPhaseIdentity::SUPERCRITICAL,
        PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT};
    const PRPhaseThermoParentView forgedSupercriticalParent(
        propaneContract, *liquidCommitted, propaneModel,
        forgedSupercriticalAdmission);
    const PRPhaseThermoOwnerResult forgedSupercritical =
        propaneAdapter.open(forgedSupercriticalParent, 250.0);
    check(forgedSupercritical.status
            == PRPhaseThermoOwnerStatus::INVALID_PARENT_PHASE_ADMISSION
            && noPublication(forgedSupercritical)
            && forgedSupercritical.diagnostics.parentStateContractValidated
            && !forgedSupercritical.diagnostics.parentFingerprintsComputed
            && !forgedSupercritical.diagnostics.kernelInvoked
            && !propaneAdapter.transactionActive(),
        "raw liquid-as-supercritical enums cannot cross parent admission");

    const std::vector<std::string> oxygenSpecies = {"O2"};
    PREquilibriumThermoModel oxygenModel(
        critical, kij, oxygenSpecies, parentOptions());
    const EquilibriumResult oxygenInitialization = oxygenModel.initializePTY(
        6.0e6, 180.0, {1.0});
    const std::optional<CompleteCellState> oxygenCommitted =
        oxygenInitialization.success() && oxygenInitialization.state
        ? completeState(*oxygenInitialization.state, oxygenModel)
        : std::nullopt;
    check(oxygenCommitted.has_value()
            && oxygenCommitted->Q.region == EquilibriumRegion::SUPERCRITICAL,
        "supercritical fixture carries an authenticated parent region");
    if (!oxygenCommitted) return 1;
    ConservativeStateContract oxygenContract;
    oxygenContract.setSpeciesLayout(oxygenSpecies, oxygenSpecies.size());
    PRPhaseRhoERecovery oxygenRecovery(
        critical, kij, oxygenSpecies, phaseOptions(160.0, 500.0));
    PRPhaseThermoOwnerAdapter oxygenAdapter(oxygenRecovery);
    const PRPhaseThermoParentAdmission oxygenAdmission{
        "cell-supercritical", 9U, 21U,
        PRPhaseIdentity::SUPERCRITICAL,
        PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT};
    const PRPhaseThermoParentView oxygenParent(
        oxygenContract, *oxygenCommitted, oxygenModel, oxygenAdmission);
    const PRPhaseThermoOwnerResult oxygenOpen = oxygenAdapter.open(
        oxygenParent, 190.0);
    std::cout << "    supercritical_open_status="
              << Foam::phaseEquilibrium::phaseThermoOwnerStatusName(
                    oxygenOpen.status)
              << " kernel="
              << Foam::phaseEquilibrium::phaseRhoEStatusName(
                    oxygenOpen.diagnostics.kernelStatus)
              << " residual_ok=" << oxygenOpen.residuals.accepted()
              << " residuals=" << oxygenOpen.residuals.pressureRelative << ','
              << oxygenOpen.residuals.temperatureRelative << ','
              << oxygenOpen.residuals.compressibilityRelative << ','
              << oxygenOpen.residuals.densityRelative << ','
              << oxygenOpen.residuals.specificInternalEnergyRelative << ','
              << oxygenOpen.residuals.specificEnthalpyRelative << ','
              << oxygenOpen.residuals.massFractionL1 << ','
              << oxygenOpen.residuals.eosRoot << ','
              << oxygenOpen.residuals.departureCaloric << ','
              << oxygenOpen.residuals.cvMinusDedTRho << ','
              << oxygenOpen.residuals.rootContinuity << '\n';
    check(oxygenOpen.success() && oxygenOpen.phaseState
            && oxygenOpen.acceptedToken
            && oxygenOpen.phaseState->root
                == Foam::phaseEquilibrium::PhaseRoot::VAPOR,
        "parent-authenticated unique-supercritical replay reaches the kernel");

    PRPhaseRhoERecovery excludedRangeRecovery(
        critical, kij, oxygenSpecies, phaseOptions(250.0, 500.0));
    PRPhaseThermoOwnerAdapter excludedRangeAdapter(excludedRangeRecovery);
    const PRPhaseThermoOwnerResult kernelFailure = excludedRangeAdapter.open(
        oxygenParent, 300.0);
    std::cout << "    excluded_open_status="
              << Foam::phaseEquilibrium::phaseThermoOwnerStatusName(
                    kernelFailure.status)
              << " kernel="
              << Foam::phaseEquilibrium::phaseRhoEStatusName(
                    kernelFailure.diagnostics.kernelStatus) << '\n';
    check(kernelFailure.status == PRPhaseThermoOwnerStatus::KERNEL_FAILURE
            && kernelFailure.diagnostics.kernelStatus
                == PRPhaseRhoEStatus::NOT_BRACKETED
            && kernelFailure.diagnostics.parentStateContractValidated
            && kernelFailure.diagnostics.parentFingerprintsComputed
            && kernelFailure.diagnostics.kernelInvoked
            && !kernelFailure.diagnostics.resultTokenIssued
            && noPublication(kernelFailure)
            && !excludedRangeAdapter.transactionActive(),
        "kernel failure rolls back without phase state or live token");

    ConservativeStateContract reorderedContract;
    reorderedContract.setSpeciesLayout(
        {"C3H8", "O2"}, species.size());
    PRPhaseThermoOwnerAdapter reorderedAdapter(recovery);
    const PRPhaseThermoParentView reorderedParent(
        reorderedContract, *committed, model, admission);
    const PRPhaseThermoOwnerResult reordered = reorderedAdapter.open(
        reorderedParent, 300.0);
    check(reordered.status
            == PRPhaseThermoOwnerStatus::INVALID_PARENT_STATE_CONTRACT
            && noPublication(reordered)
            && !reordered.diagnostics.parentFingerprintsComputed
            && !reordered.diagnostics.kernelInvoked,
        "ordered StateContract species drift fails before fingerprinting");

    if (failures != 0)
    {
        std::cerr << "PR_PHASE_THERMO_OWNER_ADAPTER_FAIL checks=" << checks
                  << " failures=" << failures << '\n';
        return 1;
    }
    std::cout << "PR_PHASE_THERMO_OWNER_ADAPTER_PASS checks=" << checks
              << " scope=disconnected-exact-parent-replay-only"
              << " Domain=0 Fluent=0 capability=0 production_ready=0\n";
    return 0;
}
