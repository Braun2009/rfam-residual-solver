#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "CompleteCellState.H"
#include "FiniteRate/RealFluidCVODEWrapper.H"
#include "FiniteRateReactor.H"
#include "MixingRules/BinaryInteractionDB.H"
#include "MixingRules/SpeciesCriticalDB.H"
#include "PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "PhaseEquilibrium/PRPhaseRhoERecovery.H"
#include "PhaseEquilibrium/PRTPFlash.H"
#include "SpeciesRegistry.H"
#include "StateContract.H"

namespace
{

using namespace Foam;
using namespace Foam::phaseEquilibrium;
using namespace Foam::thermo;

int checks = 0;
int failures = 0;

void check(bool condition, const std::string& label)
{
    ++checks;
    if (!condition)
    {
        ++failures;
        std::cerr << "FAIL: " << label << '\n';
    }
}

class ControlledChemistry final : public ChemistryModel
{
    bool active_;

public:
    mutable int calls = 0;

    ControlledChemistry(const SpeciesRegistry& registry, bool active)
        : active_(active)
    {
        nSpecies_ = 1;
        nReactions_ = active ? 1 : 0;
        speciesNames_ = {"O2"};
        molecularWeights_ = {
            registry.molecularWeight(registry.index("O2"))};
    }

    void computeSourceTerms(
        double,
        double,
        double,
        const std::vector<double>&,
        std::vector<double>& omega,
        double& heatRelease) const override
    {
        ++calls;
        omega.assign(1, 0.0);
        heatRelease = active_ ? 1.0e8 : 0.0;
    }
};

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
    options.minimumTemperature = 160.0;
    options.maximumTemperature = 500.0;
    options.temperatureScanIntervals = 128;
    options.maximumIterations = 80;
    options.derivativeRelativeStep = 3.0e-8;
    return options;
}

PRPhaseDisplacementBounds phaseDisplacementBounds()
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

std::optional<CompleteCellState> completeState(
    const EquilibriumState& equilibrium,
    const PREquilibriumThermoModel& model)
{
    PrimitiveVars primitive(
        static_cast<int>(equilibrium.overallMassFractions.size()));
    primitive.p = equilibrium.p;
    primitive.T = equilibrium.T;
    primitive.U = {0.0, 0.0, 0.0};
    primitive.Y = equilibrium.overallMassFractions;
    const CompleteCellConversionResult converted =
        CompleteStateConversion::fromPrimitiveAndEquilibrium(
            primitive, equilibrium, model);
    return converted.success() ? converted.state : std::nullopt;
}

} // namespace

int main()
{
    SpeciesRegistry registry;
    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    BinaryInteractionDB kij = BinaryInteractionDB::createDefault();
    const std::vector<std::string> species = registry.names();
    for (std::size_t first = 0; first < species.size(); ++first)
    {
        for (std::size_t second = first + 1;
             second < species.size(); ++second)
        {
            if (!kij.has(species[first], species[second]))
                kij.add(species[first], species[second], 0.0);
        }
    }
    ConservativeStateContract stateContract;
    stateContract.setSpeciesLayout(species, species.size());

    std::vector<double> feedMoleFractions(species.size(), 0.0);
    feedMoleFractions[registry.index("O2")] = 0.2;
    feedMoleFractions[registry.index("C3H8")] = 0.8;
    PRTPFlash flash(critical, kij, species);
    const TPFlashResult feed = flash.flashMoleFractions(
        2.0e6, 300.0, feedMoleFractions);
    check(feed.valid()
            && feed.region
                == Foam::phaseEquilibrium::EquilibriumRegion::TWO_PHASE,
        "two-phase vapor replay fixture is available");
    if (!feed.valid()
        || feed.region
            != Foam::phaseEquilibrium::EquilibriumRegion::TWO_PHASE)
        return 1;

    auto twoPhaseModel = std::make_shared<PREquilibriumThermoModel>(
        critical, kij, species, parentOptions());
    const EquilibriumResult initialized = twoPhaseModel->initializePTY(
        2.0e6, 300.0, feed.feedMassFractions);
    const std::optional<CompleteCellState> committed =
        initialized.success() && initialized.state
        ? completeState(*initialized.state, *twoPhaseModel)
        : std::nullopt;
    check(committed.has_value()
            && committed->Q.region
                == Foam::thermo::EquilibriumRegion::TWO_PHASE
            && committed->Q.vapor.has_value(),
        "complete parent carries an authoritative vapor payload");
    if (!committed || !committed->Q.vapor) return 1;

    PRPhaseRhoERecovery twoPhaseRecovery(
        critical, kij, species, phaseOptions());
    const PRPhaseThermoParentAdmission vaporAdmission{
        "phase-cvode/vapor-cell", 1U, 1U,
        PRPhaseIdentity::VAPOR, PRPhaseRootIdentity::PR_VAPOR_ROOT};
    const PRPhaseThermoParentView vaporParent(
        stateContract, *committed, *twoPhaseModel, vaporAdmission);

    ControlledChemistry inert(registry, false);
    SpeciesMapper inertMapper(registry, inert.speciesNames());
    RealFluidCVODEWrapper replay(
        &inert, registry, inertMapper, twoPhaseModel, twoPhaseRecovery);
    RealFluidCVODEPhaseState replayInitial;
    const bool replayOpened = replay.initPhase(
        0.0, vaporParent, 300.0, replayInitial);
    check(replayOpened && replay.phaseTransactionActive(),
        "wrapper opens one token-bound vapor transaction");
    if (!replayOpened) return 1;
    std::vector<double> replayY(species.size(), -1.0);
    double replayEnergy = -1.0;
    double replayTemperature = -1.0;
    double replayPressure = -1.0;
    check(replay.advanceTo(
            0.0, replayY, replayEnergy,
            replayTemperature, replayPressure),
        "zero-step phase replay succeeds without a CVODE RHS");
    check(replayY == replayInitial.massFractions
            && replayEnergy == replayInitial.specificEnergy
            && replayTemperature == replayInitial.temperature
            && replayPressure == replayInitial.pressure
            && inert.calls == 0
            && !replay.phaseTransactionActive()
            && replay.lastFailureClass()
                == RealFluidCVODEFailureClass::NONE,
        "zero-step replay publishes the exact local phase and closes token");

    CompleteCellState driftedCommitted = *committed;
    const PRPhaseThermoParentView driftedParent(
        stateContract, driftedCommitted, *twoPhaseModel, vaporAdmission);
    RealFluidCVODEPhaseState driftedInitial;
    check(replay.initPhase(0.0, driftedParent, 300.0, driftedInitial),
        "zero-step drift regression opens from an authoritative parent");
    driftedCommitted.Q.p *= 1.001;
    std::vector<double> driftedY(species.size(), -2.0);
    const std::vector<double> driftedYBefore = driftedY;
    double driftedEnergy = -2.0;
    double driftedTemperature = -2.0;
    double driftedPressure = -2.0;
    check(!replay.advanceTo(
            0.0, driftedY, driftedEnergy,
            driftedTemperature, driftedPressure)
            && replay.lastFailureIsPermanent()
            && !replay.phaseTransactionActive(),
        "zero-step publication rejects parent drift after token issue");
    check(driftedY == driftedYBefore && driftedEnergy == -2.0
            && driftedTemperature == -2.0 && driftedPressure == -2.0
            && inert.calls == 0,
        "zero-step parent drift publishes no output and executes no chemistry");

    RealFluidCVODEPhaseState reopenedInitial;
    check(replay.initPhase(0.0, vaporParent, 300.0, reopenedInitial),
        "wrapper can open a fresh transaction after permanent drift failure");
    std::vector<double> reopenedY(species.size(), -3.0);
    double reopenedEnergy = -3.0;
    double reopenedTemperature = -3.0;
    double reopenedPressure = -3.0;
    check(replay.advanceTo(
            0.0, reopenedY, reopenedEnergy,
            reopenedTemperature, reopenedPressure)
            && reopenedY == reopenedInitial.massFractions
            && reopenedEnergy == reopenedInitial.specificEnergy
            && reopenedTemperature == reopenedInitial.temperature
            && reopenedPressure == reopenedInitial.pressure
            && !replay.phaseTransactionActive()
            && replay.lastFailureClass()
                == RealFluidCVODEFailureClass::NONE,
        "fresh transaction remains usable after prior fail-closed cleanup");

    auto cloneModel = std::make_shared<PREquilibriumThermoModel>(
        critical, kij, species, parentOptions());
    const PRPhaseThermoParentView cloneParent(
        stateContract, *committed, *cloneModel, vaporAdmission);
    RealFluidCVODEPhaseState rejectedInitial;
    check(!replay.initPhase(0.0, cloneParent, 300.0, rejectedInitial)
            && replay.lastFailureIsPermanent()
            && !replay.phaseTransactionActive(),
        "same-configuration model clone cannot replace the shared owner");

    ControlledChemistry active(registry, true);
    SpeciesMapper activeMapper(registry, active.speciesNames());
    RealFluidCVODEWrapper displaced(
        &active, registry, activeMapper, twoPhaseModel, twoPhaseRecovery);
    RealFluidCVODEPhaseState displacedInitial;
    const bool displacementOpened = displaced.initPhase(
        0.0, vaporParent, 300.0, displacedInitial);
    check(displacementOpened,
        "active-source regression opens from exact vapor replay");
    if (!displacementOpened) return 1;
    std::vector<double> displacedY = displacedInitial.massFractions;
    const std::vector<double> displacedYBefore = displacedY;
    double displacedEnergy = displacedInitial.specificEnergy;
    const double displacedEnergyBefore = displacedEnergy;
    double displacedTemperature = displacedInitial.temperature;
    const double displacedTemperatureBefore = displacedTemperature;
    double displacedPressure = displacedInitial.pressure;
    const double displacedPressureBefore = displacedPressure;
    check(!displaced.advanceTo(
            1.0e-8, displacedY, displacedEnergy,
            displacedTemperature, displacedPressure)
            && displaced.lastFailureIsPermanent()
            && !displaced.phaseTransactionActive()
            && active.calls > 0,
        "first nonzero energy displacement is an unrecoverable CVODE failure");
    check(displacedY == displacedYBefore
            && displacedEnergy == displacedEnergyBefore
            && displacedTemperature == displacedTemperatureBefore
            && displacedPressure == displacedPressureBefore,
        "failed phase CVODE attempt publishes no caller output");

    PRPhaseDisplacementAuthority displacementAuthority(
        critical, kij, species, phaseOptions(), phaseDisplacementBounds());
    check(displacementAuthority.configured()
            && displacementAuthority.physicalModelAuthorityFingerprint()
                == twoPhaseRecovery.physicalModelAuthorityFingerprint(),
        "authorized phase CVODE uses the same physical model authority");
    RealFluidCVODEWrapper authorized(
        &active, registry, activeMapper, twoPhaseModel, twoPhaseRecovery,
        displacementAuthority);
    RealFluidCVODEPhaseState authorizedInitial;
    const bool authorizedOpened = authorized.initPhase(
        0.0, vaporParent, 300.0, authorizedInitial);
    check(authorizedOpened,
        "bounded phase CVODE opens one token-bound transaction");
    if (!authorizedOpened) return 1;
    std::vector<double> authorizedY = authorizedInitial.massFractions;
    double authorizedEnergy = authorizedInitial.specificEnergy;
    double authorizedTemperature = authorizedInitial.temperature;
    double authorizedPressure = authorizedInitial.pressure;
    const bool authorizedAdvanced = authorized.advanceTo(
        1.0e-12, authorizedY, authorizedEnergy,
        authorizedTemperature, authorizedPressure);
    if (!authorizedAdvanced)
    {
        std::cerr << "authorized phase CVODE failed owner="
                  << phaseThermoOwnerStatusName(
                        authorized.lastPhaseOwnerStatus())
                  << " displacement="
                  << phaseDisplacementStatusName(
                        authorized.lastPhaseDisplacementStatus())
                  << " probes=" << authorized.phaseProbeCount()
                  << " accepts=" << authorized.phaseAcceptCount() << '\n';
    }
    check(authorizedAdvanced
            && authorized.phaseProbeCount() > 1U
            && authorized.phaseAcceptCount() == 1U
            && !authorized.phaseTransactionActive()
            && authorized.lastFailureClass()
                == RealFluidCVODEFailureClass::NONE
            && authorizedEnergy > authorizedInitial.specificEnergy,
        "CVODE probes remain pure and one final phase state consumes the token");

    FiniteRateReactor reactor(
        &active, registry, twoPhaseModel, twoPhaseRecovery);
    std::vector<double> reactorY =
        committed->Q.vapor->massFractions;
    const std::vector<double> reactorYBefore = reactorY;
    double reactorDelta = 123.0;
    double reactorTemperature = committed->Q.T;
    const double reactorTemperatureBefore = reactorTemperature;
    double reactorPressure = committed->Q.p;
    const double reactorPressureBefore = reactorPressure;
    const int callsBeforeReactor = active.calls;
    check(!reactor.reactPhaseConservative(
            vaporParent, reactorY, 1.0e-8, reactorDelta,
            reactorTemperature, reactorPressure)
            && reactor.hasPhaseThermoOwner()
            && reactor.lastFailureCode() == 6
            && reactor.lastSubcycles() == -1,
        "reactor classifies closed-token displacement as permanent admission");
    check(reactorY == reactorYBefore && reactorDelta == 0.0
            && reactorTemperature == reactorTemperatureBefore
            && reactorPressure == reactorPressureBefore
            && active.calls > callsBeforeReactor
            && active.calls - callsBeforeReactor < 10,
        "permanent admission failure rolls back without 40 refinements");

    std::vector<double> zeroY = reactorYBefore;
    double zeroDelta = 456.0;
    double zeroTemperature = reactorTemperatureBefore;
    double zeroPressure = reactorPressureBefore;
    const int callsBeforeZero = active.calls;
    check(reactor.reactPhaseConservative(
            vaporParent, zeroY, 0.0, zeroDelta,
            zeroTemperature, zeroPressure)
            && zeroY == reactorYBefore && zeroDelta == 0.0
            && active.calls == callsBeforeZero,
        "reactor zero-step replay closes without evaluating chemistry");

    if (failures != 0)
    {
        std::cerr << "REAL_FLUID_PHASE_OWNER_CVODE_FAIL checks=" << checks
                  << " failures=" << failures << '\n';
        return 1;
    }
    std::cout << "REAL_FLUID_PHASE_OWNER_CVODE_PASS checks=" << checks
              << " scope=pure-probe-final-accept"
              << " displacement=1 Domain=0 capability=0 production_ready=0\n";
    return 0;
}
