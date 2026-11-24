#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "PhaseEquilibrium/HEMStateOwnershipContract.H"
#include "PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "PhaseEquilibrium/PRTPFlash.H"
#include "MixingRules/BinaryInteractionDB.H"
#include "MixingRules/SpeciesCriticalDB.H"
#include "SpeciesRegistry.H"

using namespace Foam;
using namespace Foam::thermo;
using namespace Foam::phaseEquilibrium;

namespace
{

int checks = 0;
int failures = 0;

void check(bool condition, const std::string& message)
{
    ++checks;
    if (!condition)
    {
        ++failures;
        std::cerr << "FAIL: " << message << "\n";
    }
}

bool close(
    double lhs, double rhs, double relative = 1.0e-8,
    double absolute = 1.0e-10)
{
    return std::isfinite(lhs) && std::isfinite(rhs)
        && std::abs(lhs - rhs) <= absolute + relative
            *std::max(std::abs(lhs), std::abs(rhs));
}

bool same(const ConservedVars& lhs, const ConservedVars& rhs)
{
    return lhs.rho == rhs.rho && lhs.rhoU == rhs.rhoU
        && lhs.rhoE == rhs.rhoE && lhs.rhoY == rhs.rhoY
        && lhs.rhoZ == rhs.rhoZ && lhs.rhoC == rhs.rhoC;
}

bool sameExactDouble(double lhs, double rhs)
{
    return std::memcmp(&lhs, &rhs, sizeof(double)) == 0;
}

bool sameExactVector(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs)
{
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t i = 0; i < lhs.size(); ++i)
        if (!sameExactDouble(lhs[i], rhs[i])) return false;
    return true;
}

bool sameExactPhase(
    const EquilibriumPhaseProperties& lhs,
    const EquilibriumPhaseProperties& rhs)
{
    return sameExactVector(lhs.massFractions, rhs.massFractions)
        && sameExactVector(lhs.moleFractions, rhs.moleFractions)
        && sameExactDouble(lhs.molecularWeight, rhs.molecularWeight)
        && sameExactDouble(lhs.rho, rhs.rho)
        && sameExactDouble(lhs.h, rhs.h)
        && sameExactDouble(lhs.e, rhs.e)
        && sameExactDouble(lhs.s, rhs.s)
        && sameExactDouble(lhs.compressibility, rhs.compressibility)
        && lhs.entropyAvailable == rhs.entropyAvailable;
}

bool sameExactOptionalPhase(
    const std::optional<EquilibriumPhaseProperties>& lhs,
    const std::optional<EquilibriumPhaseProperties>& rhs)
{
    return lhs.has_value() == rhs.has_value()
        && (!lhs.has_value() || sameExactPhase(*lhs, *rhs));
}

bool sameExactThermodynamics(
    const std::optional<EquilibriumThermodynamicPayload>& lhs,
    const std::optional<EquilibriumThermodynamicPayload>& rhs)
{
    if (lhs.has_value() != rhs.has_value()) return false;
    if (!lhs.has_value()) return true;
    return lhs->version == rhs->version
        && lhs->closure == rhs->closure
        && sameExactDouble(lhs->Cp, rhs->Cp)
        && sameExactDouble(lhs->Cv, rhs->Cv)
        && sameExactDouble(lhs->dpdrho_T, rhs->dpdrho_T)
        && sameExactDouble(lhs->dpdT_rho, rhs->dpdT_rho)
        && sameExactDouble(lhs->dedrho_T, rhs->dedrho_T)
        && sameExactDouble(lhs->dedT_rho, rhs->dedT_rho);
}

bool sameExactEquilibriumState(
    const EquilibriumState& lhs,
    const EquilibriumState& rhs)
{
    return lhs.region == rhs.region
        && sameExactVector(
            lhs.overallMassFractions, rhs.overallMassFractions)
        && sameExactDouble(lhs.p, rhs.p)
        && sameExactDouble(lhs.T, rhs.T)
        && sameExactDouble(lhs.rho, rhs.rho)
        && sameExactDouble(lhs.h, rhs.h)
        && sameExactDouble(lhs.e, rhs.e)
        && sameExactDouble(lhs.s, rhs.s)
        && sameExactDouble(
            lhs.vaporMoleFraction, rhs.vaporMoleFraction)
        && sameExactDouble(
            lhs.vaporMassFraction, rhs.vaporMassFraction)
        && sameExactDouble(
            lhs.vaporVolumeFraction, rhs.vaporVolumeFraction)
        && lhs.entropyAvailable == rhs.entropyAvailable
        && lhs.acoustic.closure == rhs.acoustic.closure
        && lhs.acoustic.evidence == rhs.acoustic.evidence
        && sameExactDouble(
            lhs.acoustic.pressureDensityDerivative,
            rhs.acoustic.pressureDensityDerivative)
        && sameExactThermodynamics(
            lhs.thermodynamics, rhs.thermodynamics)
        && sameExactOptionalPhase(lhs.liquid, rhs.liquid)
        && sameExactOptionalPhase(lhs.vapor, rhs.vapor)
        && sameExactOptionalPhase(lhs.singlePhase, rhs.singlePhase);
}

BinaryInteractionDB completeTestInteractions(const SpeciesRegistry& registry)
{
    BinaryInteractionDB result = BinaryInteractionDB::createDefault();
    const std::vector<std::string> names = registry.names();
    for (std::size_t i = 0; i < names.size(); ++i)
        for (std::size_t j = i + 1; j < names.size(); ++j)
            if (!result.has(names[i], names[j])) result.add(names[i], names[j], 0.0);
    return result;
}

RhoEFlashOptions twoPhaseOptions()
{
    RhoEFlashOptions options;
    options.minimumPressure = 1.0e6;
    options.maximumPressure = 4.0e6;
    options.minimumTemperature = 250.0;
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

class FailingReflashModel final : public EquilibriumThermoModel
{
    const EquilibriumThermoModel& delegate_;
    int failOnRecoverCall_;
    mutable int recoverCalls_ = 0;

public:
    explicit FailingReflashModel(
        const EquilibriumThermoModel& delegate,
        int failOnRecoverCall = 1)
        : delegate_(delegate), failOnRecoverCall_(failOnRecoverCall)
    {
    }

    EquilibriumResult recoverRhoE(
        double density, double internalEnergy,
        const std::vector<double>& overallMassFractions,
        const EquilibriumSeed* seed = nullptr) const override
    {
        ++recoverCalls_;
        if (recoverCalls_ != failOnRecoverCall_)
            return delegate_.recoverRhoE(
                density, internalEnergy, overallMassFractions, seed);

        EquilibriumResult result;
        result.status = EquilibriumStatus::NOT_CONVERGED;
        result.diagnostics.modelId = delegate_.modelId();
        result.diagnostics.backendPhaseStatus = "TEST_FORCED_REFLASH_FAILURE";
        return result;
    }

    EquilibriumResult initializePTY(
        double pressure, double temperature,
        const std::vector<double>& overallMassFractions) const override
    {
        return delegate_.initializePTY(
            pressure, temperature, overallMassFractions);
    }

    const std::string& modelId() const override { return delegate_.modelId(); }
    const std::vector<std::string>& speciesNames() const override
    {
        return delegate_.speciesNames();
    }
    EnergyConvention energyConvention() const override
    {
        return delegate_.energyConvention();
    }
    bool stateValidForModel(const EquilibriumState& state) const override
    {
        return delegate_.stateValidForModel(state);
    }
    int recoverCalls() const { return recoverCalls_; }
};

enum class BoxCStepStatus
{
    SUCCESS,
    INVALID_INPUT,
    PREDICTOR_FAILED,
    CORRECTOR_FAILED,
    ENERGY_LEDGER_FAILED
};

struct BoxCState
{
    ConservedVars conserved;
    EquilibriumState equilibrium;
    double time = 0.0;
    double reactedMolesPerMixtureVolume = 0.0;
    double energyLedger = 0.0;
    std::size_t acceptedSteps = 0;
    std::size_t rejectedSteps = 0;
};

class NonproductionBoxCVaporReaction
{
    const SpeciesRegistry& registry_;
    HEMRhoEStateOwnershipTransaction transaction_;
    int oxygen_;
    int fuel_;
    int carbonMonoxide_;
    int water_;
    double rateConstant_;
    double heatPerMole_;
    BoxCState committed_;

    HEMVaporConservativeIncrement sourceForExtent(
        double extentPerVaporVolume) const
    {
        HEMVaporConservativeIncrement source;
        source.deltaRhoYPerVaporVolume.assign(
            static_cast<std::size_t>(registry_.size()),0.0);
        source.deltaRhoYPerVaporVolume[fuel_] =
            -extentPerVaporVolume*registry_.molecularWeight(fuel_);
        source.deltaRhoYPerVaporVolume[oxygen_] =
            -3.5*extentPerVaporVolume*registry_.molecularWeight(oxygen_);
        source.deltaRhoYPerVaporVolume[carbonMonoxide_] =
            3.0*extentPerVaporVolume
                *registry_.molecularWeight(carbonMonoxide_);
        source.deltaRhoYPerVaporVolume[water_] =
            4.0*extentPerVaporVolume*registry_.molecularWeight(water_);
        source.deltaRhoEPerVaporVolume =
            heatPerMole_*extentPerVaporVolume;
        return source;
    }

    double vaporExtentRate(const EquilibriumState& equilibrium) const
    {
        if (equilibrium.region!=thermo::EquilibriumRegion::TWO_PHASE
            || !equilibrium.vapor.has_value())
            return std::numeric_limits<double>::quiet_NaN();
        const EquilibriumPhaseProperties& vapor=*equilibrium.vapor;
        if (vapor.massFractions.size()
                !=static_cast<std::size_t>(registry_.size())
            || !std::isfinite(vapor.rho) || vapor.rho<=0.0)
            return std::numeric_limits<double>::quiet_NaN();
        const double fuelMoles=vapor.rho*vapor.massFractions[fuel_]
            /registry_.molecularWeight(fuel_);
        const double oxygenMoles=vapor.rho*vapor.massFractions[oxygen_]
            /registry_.molecularWeight(oxygen_);
        const double limitingMoles=std::min(fuelMoles,oxygenMoles/3.5);
        return std::isfinite(limitingMoles) && limitingMoles>=0.0
            ? rateConstant_*limitingMoles
            : std::numeric_limits<double>::quiet_NaN();
    }

    BoxCStepStatus reject(BoxCStepStatus status)
    {
        ++committed_.rejectedSteps;
        return status;
    }

public:
    NonproductionBoxCVaporReaction(
        const EquilibriumThermoModel& model,
        const SpeciesRegistry& registry,
        const ConservedVars& initialConserved,
        const EquilibriumState& initialEquilibrium,
        double rateConstant,
        double heatPerMole)
    :
        registry_(registry),
        transaction_(model),
        oxygen_(registry.index("O2")),
        fuel_(registry.index("C3H8")),
        carbonMonoxide_(registry.index("CO")),
        water_(registry.index("H2O")),
        rateConstant_(rateConstant),
        heatPerMole_(heatPerMole)
    {
        committed_.conserved=initialConserved;
        committed_.equilibrium=initialEquilibrium;
    }

    BoxCStepStatus advanceRK2(double dt)
    {
        if (!std::isfinite(dt) || dt<=0.0
            || !std::isfinite(rateConstant_) || rateConstant_<0.0
            || !std::isfinite(heatPerMole_)
            || !committed_.equilibrium.valid()
            || committed_.equilibrium.region
                !=thermo::EquilibriumRegion::TWO_PHASE
            || !std::isfinite(committed_.equilibrium.vaporVolumeFraction)
            || committed_.equilibrium.vaporVolumeFraction<=0.0)
            return reject(BoxCStepStatus::INVALID_INPUT);

        const double alphaInitial=committed_.equilibrium.vaporVolumeFraction;
        const double rateInitial=vaporExtentRate(committed_.equilibrium);
        if (!std::isfinite(rateInitial))
            return reject(BoxCStepStatus::INVALID_INPUT);

        const HEMStateOwnershipResult predictor=
            transaction_.stageVaporConservativeUpdate(
                committed_.conserved,committed_.equilibrium,
                sourceForExtent(dt*rateInitial));
        if (!predictor.success())
            return reject(BoxCStepStatus::PREDICTOR_FAILED);

        const double alphaPredictor=predictor.reflashed->vaporVolumeFraction;
        const double ratePredictor=vaporExtentRate(*predictor.reflashed);
        if (!std::isfinite(alphaPredictor) || alphaPredictor<=0.0
            || !std::isfinite(ratePredictor))
            return reject(BoxCStepStatus::PREDICTOR_FAILED);

        const double mixtureExtent=0.5*dt
            *(alphaInitial*rateInitial+alphaPredictor*ratePredictor);
        const double correctedExtentPerInitialVapor=
            mixtureExtent/alphaInitial;
        const HEMStateOwnershipResult corrected=
            transaction_.stageVaporConservativeUpdate(
                committed_.conserved,committed_.equilibrium,
                sourceForExtent(correctedExtentPerInitialVapor));
        if (!corrected.success())
            return reject(BoxCStepStatus::CORRECTOR_FAILED);

        const double energyIncrement=heatPerMole_*mixtureExtent;
        if (!close(
                corrected.candidate->rhoE-committed_.conserved.rhoE,
                energyIncrement,1.0e-10,1.0e-7))
            return reject(BoxCStepStatus::ENERGY_LEDGER_FAILED);

        committed_.conserved=*corrected.candidate;
        committed_.equilibrium=*corrected.reflashed;
        committed_.time+=dt;
        committed_.reactedMolesPerMixtureVolume+=mixtureExtent;
        committed_.energyLedger+=energyIncrement;
        ++committed_.acceptedSteps;
        return BoxCStepStatus::SUCCESS;
    }

    const BoxCState& state() const { return committed_; }
    int fuelIndex() const { return fuel_; }
};

struct BoxCTrajectory
{
    bool success = false;
    BoxCState finalState;
};

BoxCTrajectory runBoxCTrajectory(
    const EquilibriumThermoModel& model,
    const SpeciesRegistry& registry,
    const ConservedVars& initialConserved,
    const EquilibriumState& initialEquilibrium,
    std::size_t stepCount,
    double totalTime,
    double rateConstant,
    double heatPerMole)
{
    BoxCTrajectory result;
    NonproductionBoxCVaporReaction box(
        model,registry,initialConserved,initialEquilibrium,
        rateConstant,heatPerMole);
    const double dt=totalTime/static_cast<double>(stepCount);
    for (std::size_t step=0; step<stepCount; ++step)
        if (box.advanceRK2(dt)!=BoxCStepStatus::SUCCESS)
        {
            result.finalState=box.state();
            return result;
        }
    result.success=true;
    result.finalState=box.state();
    return result;
}

double convergenceRatio(double coarse, double medium, double fine)
{
    const double denominator=std::abs(medium-fine);
    return denominator>0.0
        ? std::abs(coarse-medium)/denominator
        : std::numeric_limits<double>::infinity();
}

bool samePublishedBoxCState(const BoxCState& lhs, const BoxCState& rhs)
{
    return same(lhs.conserved,rhs.conserved)
        && sameExactEquilibriumState(lhs.equilibrium, rhs.equilibrium)
        && sameExactDouble(lhs.time, rhs.time)
        && sameExactDouble(
            lhs.reactedMolesPerMixtureVolume,
            rhs.reactedMolesPerMixtureVolume)
        && sameExactDouble(lhs.energyLedger, rhs.energyLedger)
        && lhs.acceptedSteps == rhs.acceptedSteps;
}

} // namespace

int main()
{
    const SpeciesRegistry registry;
    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB testKij = completeTestInteractions(registry);

    PRTPFlash feedFlash(critical, testKij, {"O2", "C3H8"});
    const TPFlashResult feed =
        feedFlash.flashMoleFractions(2.0e6, 300.0, {0.2, 0.8});
    check(feed.valid()
            && feed.region == Foam::phaseEquilibrium::EquilibriumRegion::TWO_PHASE,
        "unvalidated PR fixture reaches a two-phase O2-C3H8 feed");
    if (failures != 0) return 1;

    std::vector<double> overallY(static_cast<std::size_t>(registry.size()), 0.0);
    constexpr double productTrace = 1.0e-4;
    const double reactantScale = 1.0 - 2.0*productTrace;
    overallY[registry.index("O2")] = reactantScale*feed.feedMassFractions[0];
    overallY[registry.index("C3H8")] = reactantScale*feed.feedMassFractions[1];
    overallY[registry.index("CO")] = productTrace;
    overallY[registry.index("H2O")] = productTrace;

    const auto model = std::make_shared<const PREquilibriumThermoModel>(
        critical, testKij, registry.names(), twoPhaseOptions());
    const EquilibriumResult initial =
        model->initializePTY(2.0e6, 300.0, overallY);
    check(initial.success() && initial.state.has_value()
            && initial.state->region == thermo::EquilibriumRegion::TWO_PHASE,
        "fixture constructs a complete committed two-phase state");
    if (failures != 0) return 1;

    const EquilibriumState& equilibrium = *initial.state;
    ConservedVars committed(registry.size());
    committed.rho = equilibrium.rho;
    committed.rhoU = {2.0, -1.0, 0.5};
    const double kinetic = 0.5*(
        std::pow(committed.rhoU[0]/committed.rho, 2)
        + std::pow(committed.rhoU[1]/committed.rho, 2)
        + std::pow(committed.rhoU[2]/committed.rho, 2));
    committed.rhoE = committed.rho*(equilibrium.e + kinetic);
    for (int species = 0; species < registry.size(); ++species)
        committed.rhoY[species] = committed.rho*overallY[species];
    committed.rhoZ = 0.125*committed.rho;
    committed.rhoC = 0.375*committed.rho;
    const ConservedVars original = committed;

    // This is an integrated vapor-volume source for
    // C3H8 + 3.5 O2 -> 3 CO + 4 H2O.  It exercises only the ownership and
    // transaction contract, not R2 kinetics or a validated kij fit.
    HEMVaporConservativeIncrement vaporSource;
    vaporSource.deltaRhoYPerVaporVolume.assign(
        static_cast<std::size_t>(registry.size()), 0.0);
    const double extent = 1.0e-4;
    vaporSource.deltaRhoYPerVaporVolume[registry.index("C3H8")] =
        -extent*registry.molecularWeight(registry.index("C3H8"));
    vaporSource.deltaRhoYPerVaporVolume[registry.index("O2")] =
        -3.5*extent*registry.molecularWeight(registry.index("O2"));
    vaporSource.deltaRhoYPerVaporVolume[registry.index("CO")] =
        3.0*extent*registry.molecularWeight(registry.index("CO"));
    vaporSource.deltaRhoYPerVaporVolume[registry.index("H2O")] =
        4.0*extent*registry.molecularWeight(registry.index("H2O"));
    vaporSource.deltaRhoEPerVaporVolume = 1.19502e6*extent;

    HEMRhoEStateOwnershipTransaction transaction(*model);
    const HEMStateOwnershipResult staged =
        transaction.stageVaporConservativeUpdate(
            committed, equilibrium, vaporSource);
    check(staged.success(), "vapor source has a conservative rho-e-Y reflash");
    if (staged.success())
    {
        const ConservedVars& candidate = *staged.candidate;
        const double alphaV = equilibrium.vaporVolumeFraction;
        check(candidate.rho == committed.rho && candidate.rhoU == committed.rhoU
                && candidate.rhoZ == committed.rhoZ && candidate.rhoC == committed.rhoC,
            "the transaction preserves non-chemical conserved ownership");
        for (int species = 0; species < registry.size(); ++species)
        {
            check(close(
                    candidate.rhoY[species],
                    committed.rhoY[species] + alphaV
                        *vaporSource.deltaRhoYPerVaporVolume[species]),
                "vapor-volume species source writes to total rhoY");
        }
        check(close(
                candidate.rhoE, committed.rhoE + alphaV
                    *vaporSource.deltaRhoEPerVaporVolume),
            "vapor-volume energy source writes to total rhoE");
        check(staged.reflashed->valid()
                && staged.reflashed->region == thermo::EquilibriumRegion::TWO_PHASE,
            "quality and phase compositions are derived by the accepted reflash");
        for (int species = 0; species < registry.size(); ++species)
        {
            const double candidateY = candidate.rhoY[species]/candidate.rho;
            check(close(
                    staged.reflashed->overallMassFractions[species], candidateY),
                "reflash returns the candidate overall mass basis");
        }
    }
    check(same(committed, original),
        "staging a candidate never mutates the committed rho-e-Y state");

    HEMVaporConservativeIncrement nonconservative = vaporSource;
    nonconservative.deltaRhoYPerVaporVolume[registry.index("O2")] += 1.0e-4;
    const HEMStateOwnershipResult rejectedSource =
        transaction.stageVaporConservativeUpdate(
            committed, equilibrium, nonconservative);
    check(rejectedSource.status
            == HEMStateOwnershipStatus::SOURCE_MASS_NOT_CONSERVATIVE
            && !rejectedSource.candidate.has_value()
            && !rejectedSource.reflashed.has_value(),
        "nonconservative vapor source is rejected before publication");
    check(same(committed, original),
        "source rejection retains the complete committed state");

    HEMVaporConservativeIncrement vaporOverdraw;
    vaporOverdraw.deltaRhoYPerVaporVolume.assign(
        static_cast<std::size_t>(registry.size()), 0.0);
    const int fuel = registry.index("C3H8");
    const int carbonMonoxide = registry.index("CO");
    const double availableVaporFuel = equilibrium.vapor->rho
        *equilibrium.vapor->massFractions[fuel];
    const double overdraw = 1.01*availableVaporFuel + 1.0e-8;
    vaporOverdraw.deltaRhoYPerVaporVolume[fuel] = -overdraw;
    vaporOverdraw.deltaRhoYPerVaporVolume[carbonMonoxide] = overdraw;
    vaporOverdraw.deltaRhoEPerVaporVolume = 0.0;
    const HEMStateOwnershipResult rejectedVaporOverdraw =
        transaction.stageVaporConservativeUpdate(
            committed, equilibrium, vaporOverdraw);
    check(rejectedVaporOverdraw.status
            == HEMStateOwnershipStatus::VAPOR_SPECIES_DEPLETED
            && !rejectedVaporOverdraw.candidate.has_value()
            && !rejectedVaporOverdraw.reflashed.has_value(),
        "vapor-side chemistry cannot consume liquid-owned species inventory");
    check(same(committed, original),
        "vapor inventory rejection retains the complete committed state");

    HEMStateOwnershipTolerances invalidLooseTolerances;
    invalidLooseTolerances.concentrationRelative = 1.0;
    check(!invalidLooseTolerances.valid(),
        "a unit relative tolerance is not a valid ownership configuration");
    const FailingReflashModel invalidToleranceModel(*model);
    const HEMRhoEStateOwnershipTransaction invalidToleranceTransaction(
        invalidToleranceModel, invalidLooseTolerances);
    const HEMStateOwnershipResult invalidToleranceResult =
        invalidToleranceTransaction.stageVaporConservativeUpdate(
            committed, equilibrium, vaporOverdraw);
    check(invalidToleranceResult.status
            == HEMStateOwnershipStatus::INVALID_CONFIGURATION
            && invalidToleranceModel.recoverCalls() == 0,
        "an invalid loose tolerance is rejected before reflash");

    HEMStateOwnershipTolerances maximumTolerances;
    maximumTolerances.concentrationAbsolute =
        HEMStateOwnershipTolerances::maximumConcentrationAbsolute;
    maximumTolerances.concentrationRelative =
        HEMStateOwnershipTolerances::maximumRelative;
    maximumTolerances.energyAbsolute =
        HEMStateOwnershipTolerances::maximumEnergyAbsolute;
    maximumTolerances.energyRelative =
        HEMStateOwnershipTolerances::maximumRelative;
    maximumTolerances.sourceMassAbsolute =
        HEMStateOwnershipTolerances::maximumSourceMassAbsolute;
    maximumTolerances.sourceMassRelative =
        HEMStateOwnershipTolerances::maximumRelative;
    check(maximumTolerances.valid(),
        "the declared maximum ownership tolerances are accepted");
    const FailingReflashModel maximumToleranceModel(*model);
    const HEMRhoEStateOwnershipTransaction maximumToleranceTransaction(
        maximumToleranceModel, maximumTolerances);
    const HEMStateOwnershipResult maximumToleranceOverdraw =
        maximumToleranceTransaction.stageVaporConservativeUpdate(
            committed, equilibrium, vaporOverdraw);
    check(maximumToleranceOverdraw.status
            == HEMStateOwnershipStatus::VAPOR_SPECIES_DEPLETED
            && maximumToleranceModel.recoverCalls() == 0,
        "the loosest valid tolerances cannot bypass vapor inventory");

    const FailingReflashModel failingModel(*model);
    const HEMRhoEStateOwnershipTransaction failingTransaction(failingModel);
    const HEMStateOwnershipResult rejectedReflash =
        failingTransaction.stageVaporConservativeUpdate(
            committed, equilibrium, vaporSource);
    check(rejectedReflash.status == HEMStateOwnershipStatus::REFLASH_FAILED
            && !rejectedReflash.candidate.has_value()
            && !rejectedReflash.reflashed.has_value()
            && failingModel.recoverCalls() == 1,
        "failed reflash rejects the whole candidate before publication");
    check(same(committed, original),
        "reflash failure retains rho-e-Y and all auxiliary conserved fields");

    constexpr double boxCRateConstant=20.0;
    constexpr double boxCHeatPerMole=1.19502e6;
    constexpr double boxCTotalTime=2.0e-3;
    const BoxCTrajectory coarse=runBoxCTrajectory(
        *model,registry,committed,equilibrium,
        4,boxCTotalTime,boxCRateConstant,boxCHeatPerMole);
    const BoxCTrajectory medium=runBoxCTrajectory(
        *model,registry,committed,equilibrium,
        8,boxCTotalTime,boxCRateConstant,boxCHeatPerMole);
    const BoxCTrajectory fine=runBoxCTrajectory(
        *model,registry,committed,equilibrium,
        16,boxCTotalTime,boxCRateConstant,boxCHeatPerMole);
    const BoxCTrajectory reference=runBoxCTrajectory(
        *model,registry,committed,equilibrium,
        32,boxCTotalTime,boxCRateConstant,boxCHeatPerMole);
    check(coarse.success && medium.success && fine.success && reference.success,
        "Box C RK2 refinements complete without source or reflash fallback");
    if (coarse.success && medium.success && fine.success && reference.success)
    {
        const int fuel=registry.index("C3H8");
        const double fuelRatioOne=convergenceRatio(
            coarse.finalState.conserved.rhoY[fuel],
            medium.finalState.conserved.rhoY[fuel],
            fine.finalState.conserved.rhoY[fuel]);
        const double fuelRatioTwo=convergenceRatio(
            medium.finalState.conserved.rhoY[fuel],
            fine.finalState.conserved.rhoY[fuel],
            reference.finalState.conserved.rhoY[fuel]);
        const double temperatureRatioOne=convergenceRatio(
            coarse.finalState.equilibrium.T,
            medium.finalState.equilibrium.T,
            fine.finalState.equilibrium.T);
        const double temperatureRatioTwo=convergenceRatio(
            medium.finalState.equilibrium.T,
            fine.finalState.equilibrium.T,
            reference.finalState.equilibrium.T);
        std::cout << "    Box C RK2 ratios fuel=" << fuelRatioOne
                  << ',' << fuelRatioTwo
                  << " T=" << temperatureRatioOne
                  << ',' << temperatureRatioTwo << '\n';
        check(fuelRatioOne>2.5 && fuelRatioTwo>2.5
                && temperatureRatioOne>2.5 && temperatureRatioTwo>2.5,
            "Box C vapor-rate/write-back/reflash sequence converges at RK2 scale");
        check(close(
                reference.finalState.conserved.rhoE-committed.rhoE,
                reference.finalState.energyLedger,1.0e-10,1.0e-6),
            "Box C total rhoE closes against the accepted reaction ledger");
        double finalSpeciesDensity=0.0;
        for (double rhoY : reference.finalState.conserved.rhoY)
            finalSpeciesDensity+=rhoY;
        check(close(finalSpeciesDensity,reference.finalState.conserved.rho)
                && reference.finalState.conserved.rho==committed.rho
                && reference.finalState.conserved.rhoU==committed.rhoU
                && reference.finalState.conserved.rhoZ==committed.rhoZ
                && reference.finalState.conserved.rhoC==committed.rhoC,
            "Box C conserves total mass and all non-chemical ledgers");
        check(reference.finalState.equilibrium.valid()
                && reference.finalState.equilibrium.overallMassFractions.size()
                    ==reference.finalState.conserved.rhoY.size()
                && reference.finalState.acceptedSteps==32
                && reference.finalState.rejectedSteps==0
                && close(reference.finalState.time,boxCTotalTime),
            "every accepted Box C step publishes only a complete reflash");
    }

    const FailingReflashModel failPredictorModel(*model,1);
    NonproductionBoxCVaporReaction predictorRollbackBox(
        failPredictorModel,registry,committed,equilibrium,
        boxCRateConstant,boxCHeatPerMole);
    const BoxCState beforePredictorFailure=predictorRollbackBox.state();
    const BoxCStepStatus predictorFailure=
        predictorRollbackBox.advanceRK2(boxCTotalTime/4.0);
    const BoxCState afterPredictorFailure=predictorRollbackBox.state();
    check(predictorFailure==BoxCStepStatus::PREDICTOR_FAILED
            && failPredictorModel.recoverCalls()==1
            && samePublishedBoxCState(
                beforePredictorFailure,afterPredictorFailure)
            && afterPredictorFailure.rejectedSteps
                ==beforePredictorFailure.rejectedSteps+1,
        "injected predictor reflash failure changes only rejectedSteps");

    const FailingReflashModel failCorrectorModel(*model,2);
    NonproductionBoxCVaporReaction correctorRollbackBox(
        failCorrectorModel,registry,committed,equilibrium,
        boxCRateConstant,boxCHeatPerMole);
    const BoxCState beforeCorrectorFailure=correctorRollbackBox.state();
    const BoxCStepStatus correctorFailure=
        correctorRollbackBox.advanceRK2(boxCTotalTime/4.0);
    const BoxCState afterCorrectorFailure=correctorRollbackBox.state();
    check(correctorFailure==BoxCStepStatus::CORRECTOR_FAILED
            && failCorrectorModel.recoverCalls()==2
            && samePublishedBoxCState(
                beforeCorrectorFailure,afterCorrectorFailure)
            && afterCorrectorFailure.rejectedSteps
                ==beforeCorrectorFailure.rejectedSteps+1,
        "injected second reflash failure changes only rejectedSteps");

    check(!beforeCorrectorFailure.equilibrium.thermodynamics.has_value()
            && beforeCorrectorFailure.equilibrium.liquid.has_value()
            && beforeCorrectorFailure.equilibrium.vapor.has_value()
            && !beforeCorrectorFailure.equilibrium.singlePhase.has_value(),
        "rollback fixture exposes liquid/vapor and absent optional thermo/single-phase payloads");

    BoxCState mutatedRegion=beforeCorrectorFailure;
    mutatedRegion.equilibrium.region=thermo::EquilibriumRegion::UNKNOWN;
    check(!samePublishedBoxCState(beforeCorrectorFailure,mutatedRegion),
        "whole-state comparator detects equilibrium-region mutation");

    BoxCState mutatedEnthalpy=beforeCorrectorFailure;
    mutatedEnthalpy.equilibrium.h=std::nextafter(
        mutatedEnthalpy.equilibrium.h,
        std::numeric_limits<double>::infinity());
    check(!samePublishedBoxCState(beforeCorrectorFailure,mutatedEnthalpy),
        "whole-state comparator detects mixture-enthalpy mutation");

    BoxCState mutatedEntropy=beforeCorrectorFailure;
    mutatedEntropy.equilibrium.s=std::nextafter(
        mutatedEntropy.equilibrium.s,
        std::numeric_limits<double>::infinity());
    check(!samePublishedBoxCState(beforeCorrectorFailure,mutatedEntropy),
        "whole-state comparator detects mixture-entropy mutation");

    BoxCState mutatedEntropyAvailability=beforeCorrectorFailure;
    mutatedEntropyAvailability.equilibrium.entropyAvailable=
        !mutatedEntropyAvailability.equilibrium.entropyAvailable;
    check(!samePublishedBoxCState(
            beforeCorrectorFailure,mutatedEntropyAvailability),
        "whole-state comparator detects entropy-availability mutation");

    BoxCState mutatedVaporMoleFraction=beforeCorrectorFailure;
    mutatedVaporMoleFraction.equilibrium.vaporMoleFraction=std::nextafter(
        mutatedVaporMoleFraction.equilibrium.vaporMoleFraction,
        std::numeric_limits<double>::infinity());
    check(!samePublishedBoxCState(
            beforeCorrectorFailure,mutatedVaporMoleFraction),
        "whole-state comparator detects vapor-mole-fraction mutation");

    BoxCState mutatedVaporMassFraction=beforeCorrectorFailure;
    mutatedVaporMassFraction.equilibrium.vaporMassFraction=std::nextafter(
        mutatedVaporMassFraction.equilibrium.vaporMassFraction,
        std::numeric_limits<double>::infinity());
    check(!samePublishedBoxCState(
            beforeCorrectorFailure,mutatedVaporMassFraction),
        "whole-state comparator detects vapor-mass-fraction mutation");

    BoxCState mutatedVaporVolumeFraction=beforeCorrectorFailure;
    mutatedVaporVolumeFraction.equilibrium.vaporVolumeFraction=std::nextafter(
        mutatedVaporVolumeFraction.equilibrium.vaporVolumeFraction,
        std::numeric_limits<double>::infinity());
    check(!samePublishedBoxCState(
            beforeCorrectorFailure,mutatedVaporVolumeFraction),
        "whole-state comparator detects vapor-volume-fraction mutation");

    BoxCState mutatedAcousticClosure=beforeCorrectorFailure;
    mutatedAcousticClosure.equilibrium.acoustic.closure=
        thermo::AcousticClosure::UNKNOWN;
    check(!samePublishedBoxCState(
            beforeCorrectorFailure,mutatedAcousticClosure),
        "whole-state comparator detects acoustic-closure mutation");

    BoxCState mutatedAcousticDerivative=beforeCorrectorFailure;
    mutatedAcousticDerivative.equilibrium.acoustic.pressureDensityDerivative=
        std::nextafter(
            mutatedAcousticDerivative.equilibrium.acoustic
                .pressureDensityDerivative,
            std::numeric_limits<double>::infinity());
    check(!samePublishedBoxCState(
            beforeCorrectorFailure,mutatedAcousticDerivative),
        "whole-state comparator detects acoustic-derivative mutation");

    BoxCState mutatedThermodynamicsPresence=beforeCorrectorFailure;
    mutatedThermodynamicsPresence.equilibrium.thermodynamics=
        EquilibriumThermodynamicPayload{};
    check(!samePublishedBoxCState(
            beforeCorrectorFailure,mutatedThermodynamicsPresence),
        "whole-state comparator detects thermodynamic-payload presence mutation");

    BoxCState thermodynamicsBaseline=beforeCorrectorFailure;
    thermodynamicsBaseline.equilibrium.thermodynamics=
        EquilibriumThermodynamicPayload{};
    BoxCState mutatedThermodynamicsField=thermodynamicsBaseline;
    mutatedThermodynamicsField.equilibrium.thermodynamics->Cp=1.0;
    check(!samePublishedBoxCState(
            thermodynamicsBaseline,mutatedThermodynamicsField),
        "whole-state comparator detects thermodynamic-payload field mutation");

    BoxCState mutatedLiquid=beforeCorrectorFailure;
    mutatedLiquid.equilibrium.liquid->s=std::nextafter(
        mutatedLiquid.equilibrium.liquid->s,
        std::numeric_limits<double>::infinity());
    check(!samePublishedBoxCState(beforeCorrectorFailure,mutatedLiquid),
        "whole-state comparator detects liquid-payload mutation");

    BoxCState mutatedVapor=beforeCorrectorFailure;
    mutatedVapor.equilibrium.vapor->molecularWeight=std::nextafter(
        mutatedVapor.equilibrium.vapor->molecularWeight,
        std::numeric_limits<double>::infinity());
    check(!samePublishedBoxCState(beforeCorrectorFailure,mutatedVapor),
        "whole-state comparator detects vapor-payload mutation");

    BoxCState mutatedSinglePhase=beforeCorrectorFailure;
    mutatedSinglePhase.equilibrium.singlePhase=
        *mutatedSinglePhase.equilibrium.vapor;
    check(!samePublishedBoxCState(beforeCorrectorFailure,mutatedSinglePhase),
        "whole-state comparator detects single-phase-payload mutation");

    if (failures != 0)
    {
        std::cerr << "PR_VAPOR_CHEMISTRY_REFLASH_CONTRACT_FAIL checks="
                  << checks << " failures=" << failures << "\n";
        return 1;
    }
    std::cout << "PR_VAPOR_CHEMISTRY_REFLASH_CONTRACT_PASS checks="
              << checks << " scope=nonproduction-ownership\n";
    return 0;
}
