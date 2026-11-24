#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ChemistrySourceOperator.H"
#include "MixingRules/BinaryInteractionDB.H"
#include "MixingRules/SpeciesCriticalDB.H"
#include "PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "SectorChemistryOperator.H"
#include "SpeciesRegistry.H"

using namespace Foam;

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
        std::cerr << "FAIL: " << message << '\n';
    }
}

bool exactScalar(double lhs, double rhs)
{
    return lhs == rhs || (std::isnan(lhs) && std::isnan(rhs));
}

bool exactVector(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs)
{
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t index = 0; index < lhs.size(); ++index)
        if (!exactScalar(lhs[index], rhs[index])) return false;
    return true;
}

bool exactSame(
    const thermo::EquilibriumPhaseProperties& lhs,
    const thermo::EquilibriumPhaseProperties& rhs)
{
    return exactVector(lhs.massFractions, rhs.massFractions)
        && exactVector(lhs.moleFractions, rhs.moleFractions)
        && exactScalar(lhs.molecularWeight, rhs.molecularWeight)
        && exactScalar(lhs.rho, rhs.rho)
        && exactScalar(lhs.h, rhs.h)
        && exactScalar(lhs.e, rhs.e)
        && exactScalar(lhs.s, rhs.s)
        && exactScalar(lhs.compressibility, rhs.compressibility)
        && lhs.entropyAvailable == rhs.entropyAvailable;
}

bool exactSame(
    const thermo::EquilibriumThermodynamicPayload& lhs,
    const thermo::EquilibriumThermodynamicPayload& rhs)
{
    return lhs.version == rhs.version
        && lhs.closure == rhs.closure
        && exactScalar(lhs.Cp, rhs.Cp)
        && exactScalar(lhs.Cv, rhs.Cv)
        && exactScalar(lhs.dpdrho_T, rhs.dpdrho_T)
        && exactScalar(lhs.dpdT_rho, rhs.dpdT_rho)
        && exactScalar(lhs.dedrho_T, rhs.dedrho_T)
        && exactScalar(lhs.dedT_rho, rhs.dedT_rho);
}

bool exactSame(
    const thermo::EquilibriumState& lhs,
    const thermo::EquilibriumState& rhs)
{
    if (lhs.liquid.has_value() != rhs.liquid.has_value()
        || lhs.vapor.has_value() != rhs.vapor.has_value()
        || lhs.singlePhase.has_value() != rhs.singlePhase.has_value()
        || lhs.thermodynamics.has_value() != rhs.thermodynamics.has_value())
        return false;
    return lhs.region == rhs.region
        && exactVector(lhs.overallMassFractions, rhs.overallMassFractions)
        && exactScalar(lhs.p, rhs.p)
        && exactScalar(lhs.T, rhs.T)
        && exactScalar(lhs.rho, rhs.rho)
        && exactScalar(lhs.h, rhs.h)
        && exactScalar(lhs.e, rhs.e)
        && exactScalar(lhs.s, rhs.s)
        && exactScalar(lhs.vaporMoleFraction, rhs.vaporMoleFraction)
        && exactScalar(lhs.vaporMassFraction, rhs.vaporMassFraction)
        && exactScalar(lhs.vaporVolumeFraction, rhs.vaporVolumeFraction)
        && lhs.entropyAvailable == rhs.entropyAvailable
        && lhs.acoustic.closure == rhs.acoustic.closure
        && lhs.acoustic.evidence == rhs.acoustic.evidence
        && exactScalar(
            lhs.acoustic.pressureDensityDerivative,
            rhs.acoustic.pressureDensityDerivative)
        && (!lhs.thermodynamics
            || exactSame(*lhs.thermodynamics, *rhs.thermodynamics))
        && (!lhs.liquid || exactSame(*lhs.liquid, *rhs.liquid))
        && (!lhs.vapor || exactSame(*lhs.vapor, *rhs.vapor))
        && (!lhs.singlePhase
            || exactSame(*lhs.singlePhase, *rhs.singlePhase));
}

bool exactSame(const ConservedVars& lhs, const ConservedVars& rhs)
{
    return exactScalar(lhs.rho, rhs.rho)
        && exactVector(lhs.rhoU, rhs.rhoU)
        && exactScalar(lhs.rhoE, rhs.rhoE)
        && exactVector(lhs.rhoY, rhs.rhoY)
        && exactScalar(lhs.rhoZ, rhs.rhoZ)
        && exactScalar(lhs.rhoC, rhs.rhoC);
}

bool exactSame(const PrimitiveVars& lhs, const PrimitiveVars& rhs)
{
    return exactScalar(lhs.p, rhs.p)
        && exactScalar(lhs.T, rhs.T)
        && exactVector(lhs.U, rhs.U)
        && exactVector(lhs.Y, rhs.Y)
        && exactScalar(lhs.Z, rhs.Z)
        && exactScalar(lhs.C, rhs.C);
}

struct ExactDomainState
{
    std::vector<ConservedVars> current;
    std::vector<ConservedVars> candidate;
    std::vector<PrimitiveVars> primitive;
    std::vector<thermo::EquilibriumState> equilibrium;
};

ExactDomainState captureExactState(const Domain& domain)
{
    ExactDomainState result;
    result.current.reserve(domain.N());
    result.candidate.reserve(domain.N());
    result.primitive.reserve(domain.N());
    result.equilibrium.reserve(domain.N());
    for (int cell = 0; cell < domain.N(); ++cell)
    {
        result.current.push_back(domain.U(cell));
        result.candidate.push_back(domain.Unew(cell));
        result.primitive.push_back(domain.W(cell));
        result.equilibrium.push_back(domain.Q(cell));
    }
    return result;
}

bool exactStateMatches(
    const Domain& domain,
    const ExactDomainState& expected)
{
    if (domain.N() != static_cast<int>(expected.current.size())) return false;
    for (int cell = 0; cell < domain.N(); ++cell)
        if (!exactSame(domain.U(cell), expected.current[cell])
            || !exactSame(domain.Unew(cell), expected.candidate[cell])
            || !exactSame(domain.W(cell), expected.primitive[cell])
            || !exactSame(domain.Q(cell), expected.equilibrium[cell]))
            return false;
    return true;
}

class CountingVaporChemistry final : public ChemistryModel
{
    static constexpr double fuelMolecularWeight_ = 0.04409562;
    static constexpr double oxygenMolecularWeight_ = 0.0319988;
    static constexpr double coMolecularWeight_ = 0.0280101;
    static constexpr double waterMolecularWeight_ = 0.01801528;
    static constexpr double reactionHeat_ = 1.19502e6;

    double rateConstant_;
    mutable std::atomic<int> calls_{0};
    mutable double observedTemperature_ = 0.0;
    mutable double observedPressure_ = 0.0;
    mutable double observedDensity_ = 0.0;
    mutable std::vector<double> observedY_;

public:
    explicit CountingVaporChemistry(
        const SpeciesRegistry& registry,
        double rateConstant = 50.0)
        : rateConstant_(rateConstant)
    {
        speciesNames_ = {"C3H8", "O2", "CO", "CO2", "H2O", "N2"};
        nSpecies_ = static_cast<int>(speciesNames_.size());
        nReactions_ = 1;
        molecularWeights_.reserve(speciesNames_.size());
        for (const std::string& name : speciesNames_)
            molecularWeights_.push_back(
                registry.molecularWeight(registry.index(name)));
    }

    void computeSourceTerms(
        double temperature,
        double pressure,
        double density,
        const std::vector<double>& massFractions,
        std::vector<double>& omega,
        double& heatRelease) const override
    {
        if (calls_.fetch_add(1, std::memory_order_relaxed) == 0)
        {
            observedTemperature_ = temperature;
            observedPressure_ = pressure;
            observedDensity_ = density;
            observedY_ = massFractions;
        }

        omega.assign(static_cast<std::size_t>(nSpecies_), 0.0);
        const int fuel = speciesIndex("C3H8");
        const int oxygen = speciesIndex("O2");
        const int carbonMonoxide = speciesIndex("CO");
        const int water = speciesIndex("H2O");
        const double fuelMoles = density
            * std::max(massFractions[fuel], 0.0) / fuelMolecularWeight_;
        const double oxygenMoles = density
            * std::max(massFractions[oxygen], 0.0)
            / oxygenMolecularWeight_;
        const double extentRate = rateConstant_
            * std::min(fuelMoles, oxygenMoles / 3.5);
        omega[fuel] = -extentRate * fuelMolecularWeight_;
        omega[oxygen] = -3.5 * extentRate * oxygenMolecularWeight_;
        omega[carbonMonoxide] = 3.0 * extentRate * coMolecularWeight_;
        omega[water] = 4.0 * extentRate * waterMolecularWeight_;
        heatRelease = extentRate * reactionHeat_;
    }

    void resetObservation() const
    {
        calls_.store(0, std::memory_order_relaxed);
        observedTemperature_ = 0.0;
        observedPressure_ = 0.0;
        observedDensity_ = 0.0;
        observedY_.clear();
    }

    int calls() const { return calls_.load(std::memory_order_relaxed); }
    double observedTemperature() const { return observedTemperature_; }
    double observedPressure() const { return observedPressure_; }
    double observedDensity() const { return observedDensity_; }
    const std::vector<double>& observedY() const { return observedY_; }
};

class MutableIdentityThermoModel final
    : public thermo::EquilibriumThermoModel
{
    std::shared_ptr<const thermo::EquilibriumThermoModel> delegate_;
    std::string modelId_;
    std::vector<std::string> speciesNames_;
    EnergyConvention energyConvention_;

    thermo::EquilibriumResult retag(
        thermo::EquilibriumResult result) const
    {
        result.diagnostics.modelId = modelId_;
        return result;
    }

public:
    explicit MutableIdentityThermoModel(
        std::shared_ptr<const thermo::EquilibriumThermoModel> delegate)
        : delegate_(std::move(delegate)),
          modelId_(delegate_->modelId() + ":mutable-identity-test"),
          speciesNames_(delegate_->speciesNames()),
          energyConvention_(delegate_->energyConvention())
    {}

    thermo::EquilibriumResult recoverRhoE(
        double density,
        double internalEnergy,
        const std::vector<double>& overallMassFractions,
        const thermo::EquilibriumSeed* seed = nullptr) const override
    {
        return retag(delegate_->recoverRhoE(
            density, internalEnergy, overallMassFractions, seed));
    }

    thermo::EquilibriumResult initializePTY(
        double pressure,
        double temperature,
        const std::vector<double>& overallMassFractions) const override
    {
        return retag(delegate_->initializePTY(
            pressure, temperature, overallMassFractions));
    }

    const std::string& modelId() const override { return modelId_; }
    const std::vector<std::string>& speciesNames() const override
    { return speciesNames_; }
    EnergyConvention energyConvention() const override
    { return energyConvention_; }
    bool stateValidForModel(
        const thermo::EquilibriumState& state) const override
    { return delegate_->stateValidForModel(state); }
    bool supportsPrimitivePTYInitialization() const override
    { return delegate_->supportsPrimitivePTYInitialization(); }

    void driftModelId() { modelId_ += ":drifted"; }
    void driftSpecies()
    {
        if (speciesNames_.size() >= 2)
            std::swap(speciesNames_[0], speciesNames_[1]);
    }
    void driftEnergyConvention()
    {
        energyConvention_ = static_cast<EnergyConvention>(999);
    }
};

class RegionOverrideThermoModel final
    : public thermo::EquilibriumThermoModel
{
    std::shared_ptr<const thermo::EquilibriumThermoModel> delegate_;
    thermo::EquilibriumRegion region_;
    std::string modelId_;

    thermo::EquilibriumResult retag(
        thermo::EquilibriumResult result) const
    {
        result.diagnostics.modelId = modelId_;
        if (result.success()) result.state->region = region_;
        return result;
    }

public:
    RegionOverrideThermoModel(
        std::shared_ptr<const thermo::EquilibriumThermoModel> delegate,
        thermo::EquilibriumRegion region)
        : delegate_(std::move(delegate)),
          region_(region),
          modelId_(delegate_->modelId() + ":forced-region-test")
    {}

    thermo::EquilibriumResult recoverRhoE(
        double density,
        double internalEnergy,
        const std::vector<double>& overallMassFractions,
        const thermo::EquilibriumSeed* seed = nullptr) const override
    {
        return retag(delegate_->recoverRhoE(
            density, internalEnergy, overallMassFractions, seed));
    }

    thermo::EquilibriumResult initializePTY(
        double pressure,
        double temperature,
        const std::vector<double>& overallMassFractions) const override
    {
        return retag(delegate_->initializePTY(
            pressure, temperature, overallMassFractions));
    }

    const std::string& modelId() const override { return modelId_; }
    const std::vector<std::string>& speciesNames() const override
    { return delegate_->speciesNames(); }
    EnergyConvention energyConvention() const override
    { return delegate_->energyConvention(); }
    bool stateValidForModel(
        const thermo::EquilibriumState& state) const override
    { return delegate_->stateValidForModel(state); }
    bool supportsPrimitivePTYInitialization() const override
    { return delegate_->supportsPrimitivePTYInitialization(); }
};

BinaryInteractionDB completeTestInteractions(
    const SpeciesRegistry& registry)
{
    BinaryInteractionDB result = BinaryInteractionDB::createDefault();
    const std::vector<std::string> names = registry.names();
    for (std::size_t first = 0; first < names.size(); ++first)
        for (std::size_t second = first + 1;
             second < names.size(); ++second)
            if (!result.has(names[first], names[second]))
                result.add(names[first], names[second], 0.0);
    return result;
}

phaseEquilibrium::RhoEFlashOptions baseOptions()
{
    phaseEquilibrium::RhoEFlashOptions options;
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

phaseEquilibrium::RhoEFlashOptions vaporOptions()
{
    phaseEquilibrium::RhoEFlashOptions options = baseOptions();
    options.minimumPressure = 1.0e5;
    options.maximumPressure = 1.0e7;
    options.minimumTemperature = 250.0;
    options.maximumTemperature = 2000.0;
    options.initialPressure = 4.0e6;
    options.initialTemperature = 800.0;
    return options;
}

std::vector<double> vaporComposition(const SpeciesRegistry& registry)
{
    std::vector<double> result(registry.size(), 0.0);
    result[registry.index("O2")] = 0.8;
    result[registry.index("C3H8")] = 0.2;
    return result;
}

std::vector<double> twoPhaseComposition(
    const SpeciesRegistry& registry,
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& interactions)
{
    phaseEquilibrium::PRTPFlash flash(
        critical, interactions, {"O2", "C3H8"});
    const phaseEquilibrium::TPFlashResult split =
        flash.flashMoleFractions(2.0e6, 300.0, {0.2, 0.8});
    if (!split.valid()
        || split.region != phaseEquilibrium::EquilibriumRegion::TWO_PHASE)
        throw std::runtime_error("two-phase fixture flash failed");

    constexpr double productTrace = 1.0e-4;
    const double reactantScale = 1.0 - 2.0 * productTrace;
    std::vector<double> result(registry.size(), 0.0);
    result[registry.index("O2")] =
        reactantScale * split.feedMassFractions[0];
    result[registry.index("C3H8")] =
        reactantScale * split.feedMassFractions[1];
    result[registry.index("CO")] = productTrace;
    result[registry.index("H2O")] = productTrace;
    return result;
}

PrimitiveVars primitiveFrom(
    const thermo::EquilibriumState& equilibrium,
    int speciesCount)
{
    PrimitiveVars primitive(speciesCount);
    primitive.p = equilibrium.p;
    primitive.T = equilibrium.T;
    primitive.U = {0.0, 0.0, 0.0};
    primitive.Y = equilibrium.overallMassFractions;
    return primitive;
}

std::unique_ptr<Domain> makeOneDimensionalDomain(
    const std::shared_ptr<const thermo::EquilibriumThermoModel>& owner,
    const thermo::EquilibriumState& equilibrium,
    int speciesCount)
{
    auto domain = std::make_unique<Domain>(
        1, 0.01, speciesCount, owner);
    const PrimitiveVars primitive = primitiveFrom(equilibrium, speciesCount);
    if (!domain->replaceCompleteField({primitive}, {equilibrium}))
        throw std::runtime_error("1-D Domain initialization failed");
    return domain;
}

std::unique_ptr<Domain> makeSectorDomain(
    const std::shared_ptr<const thermo::EquilibriumThermoModel>& owner,
    const thermo::EquilibriumState& equilibrium,
    int speciesCount)
{
    auto domain = std::make_unique<Domain>(
        3, 3, 3, 0.03, 0.03, 0.03, speciesCount, owner);
    const PrimitiveVars primitive = primitiveFrom(equilibrium, speciesCount);
    if (!domain->replaceCompleteField(
            std::vector<PrimitiveVars>(domain->N(), primitive),
            std::vector<thermo::EquilibriumState>(
                domain->N(), equilibrium)))
        throw std::runtime_error("Sector Domain initialization failed");
    return domain;
}

template<class Callable>
bool rejects(Callable&& callable)
{
    try
    {
        callable();
    }
    catch (const std::runtime_error&)
    {
        return true;
    }
    return false;
}

struct SectorLedger
{
    long active;
    long failures;
    long retried;
    long limited;
    long twoPhaseActive;
    long liquidSkipped;
    long phaseReflashFailures;
    long cvodeFailures;
    long stateFailures;
    long temperatureFailures;
    long massFailures;
    long elementFailures;
    double failureElementError;
    int failureI;
    int failureJ;
    int failureK;
    double failureDt;
    double failureRho;
    double failurePressure;
    double failureTemperature;
    std::vector<double> failureChemY;
    double heatMinimum;
    double heatAverage;
    double heatMaximum;
};

SectorLedger captureLedger(const SectorChemistryOperator& op)
{
    return {
        op.lastActiveCells(),
        op.lastFailures(),
        op.lastRetriedCells(),
        op.lastPositivityLimitedCells(),
        op.lastTwoPhaseActiveCells(),
        op.lastLiquidSkippedCells(),
        op.lastPhaseReflashFailures(),
        op.lastCvodeFailures(),
        op.lastStateFailures(),
        op.lastTemperatureFailures(),
        op.lastMassFailures(),
        op.lastElementFailures(),
        op.lastFailureElementError(),
        op.lastFailureI(),
        op.lastFailureJ(),
        op.lastFailureK(),
        op.lastFailureDt(),
        op.lastFailureRho(),
        op.lastFailurePressure(),
        op.lastFailureTemperature(),
        op.lastFailureChemY(),
        op.lastHeatReleaseMinimum(),
        op.lastHeatReleaseAverage(),
        op.lastHeatReleaseMaximum()};
}

bool exactSame(const SectorLedger& lhs, const SectorLedger& rhs)
{
    return lhs.active == rhs.active
        && lhs.failures == rhs.failures
        && lhs.retried == rhs.retried
        && lhs.limited == rhs.limited
        && lhs.twoPhaseActive == rhs.twoPhaseActive
        && lhs.liquidSkipped == rhs.liquidSkipped
        && lhs.phaseReflashFailures == rhs.phaseReflashFailures
        && lhs.cvodeFailures == rhs.cvodeFailures
        && lhs.stateFailures == rhs.stateFailures
        && lhs.temperatureFailures == rhs.temperatureFailures
        && lhs.massFailures == rhs.massFailures
        && lhs.elementFailures == rhs.elementFailures
        && exactScalar(lhs.failureElementError, rhs.failureElementError)
        && lhs.failureI == rhs.failureI
        && lhs.failureJ == rhs.failureJ
        && lhs.failureK == rhs.failureK
        && exactScalar(lhs.failureDt, rhs.failureDt)
        && exactScalar(lhs.failureRho, rhs.failureRho)
        && exactScalar(lhs.failurePressure, rhs.failurePressure)
        && exactScalar(lhs.failureTemperature, rhs.failureTemperature)
        && exactVector(lhs.failureChemY, rhs.failureChemY)
        && exactScalar(lhs.heatMinimum, rhs.heatMinimum)
        && exactScalar(lhs.heatAverage, rhs.heatAverage)
        && exactScalar(lhs.heatMaximum, rhs.heatMaximum);
}

void exerciseChemistryMismatch(
    ChemistrySourceOperator& op,
    Domain& domain,
    const CountingVaporChemistry& chemistry,
    const std::string& label)
{
    domain.beginStep();
    const ExactDomainState before = captureExactState(domain);
    const int callsBefore = chemistry.calls();
    check(
        rejects([&]() { static_cast<void>(op.maxStableDt(domain, 0.1)); }),
        label + " maxStableDt rejects");
    check(
        chemistry.calls() == callsBefore,
        label + " maxStableDt rejects before RHS");
    check(
        exactStateMatches(domain, before),
        label + " maxStableDt leaves Domain exact");
    check(
        rejects([&]() { op.apply(domain, 1.0e-8); }),
        label + " apply rejects");
    check(
        chemistry.calls() == callsBefore,
        label + " apply rejects before RHS");
    check(
        exactStateMatches(domain, before),
        label + " apply leaves U/Unew/W/Q exact");
}

void exerciseSectorMismatch(
    SectorChemistryOperator& op,
    Domain& domain,
    const CountingVaporChemistry& chemistry,
    const SectorLedger& expectedLedger,
    const std::string& label)
{
    domain.beginStep();
    const ExactDomainState before = captureExactState(domain);
    const int callsBefore = chemistry.calls();
    check(
        rejects([&]() { static_cast<void>(op.maxStableDt(domain, 0.1)); }),
        label + " maxStableDt rejects");
    check(
        chemistry.calls() == callsBefore,
        label + " maxStableDt rejects before RHS");
    check(
        exactStateMatches(domain, before),
        label + " maxStableDt leaves Domain exact");
    check(
        exactSame(captureLedger(op), expectedLedger),
        label + " maxStableDt preserves Sector ledger");
    check(
        rejects([&]() { op.apply(domain, 1.0e-8); }),
        label + " apply rejects");
    check(
        chemistry.calls() == callsBefore,
        label + " apply rejects before RHS");
    check(
        exactStateMatches(domain, before),
        label + " apply leaves U/Unew/W/Q exact");
    check(
        exactSame(captureLedger(op), expectedLedger),
        label + " apply preserves Sector counters/failure ledger");
}

void testExactOwnerAndPointerAdmission(
    const SpeciesRegistry& registry,
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& interactions)
{
    const auto owner =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical, interactions, registry.names(), vaporOptions());
    const auto ownerClone =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical, interactions, registry.names(), vaporOptions());
    BinaryInteractionDB changedInteractions = interactions;
    changedInteractions.add("O2", "C3H8", 0.125);
    const auto changedOwner =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical, changedInteractions, registry.names(), vaporOptions());
    const std::vector<double> composition = vaporComposition(registry);
    const thermo::EquilibriumResult ownerState =
        owner->initializePTY(4.0e6, 800.0, composition);
    const thermo::EquilibriumResult cloneState =
        ownerClone->initializePTY(4.0e6, 800.0, composition);
    const thermo::EquilibriumResult changedState =
        changedOwner->initializePTY(4.0e6, 800.0, composition);

    check(ownerState.success(), "exact owner single-vapor fixture succeeds");
    check(cloneState.success(), "clone owner single-vapor fixture succeeds");
    check(changedState.success(), "changed-kij owner fixture succeeds");
    if (!ownerState.success() || !cloneState.success()
        || !changedState.success())
        return;
    check(
        ownerState.state->region == thermo::EquilibriumRegion::SINGLE_VAPOR,
        "exact owner fixture is single vapor");
    check(
        ownerState.diagnostics.modelId == owner->modelId(),
        "exact owner fixture diagnostics carry owner ID");
    check(
        !ownerState.state->thermodynamics.has_value(),
        "exact owner fixture exercises committed-energy dt fallback");
    check(
        owner->stateValidForModel(*ownerState.state),
        "exact owner validates its fixture state");
    check(
        owner->modelId() == ownerClone->modelId()
            && owner.get() != ownerClone.get(),
        "equivalent clone has same model ID but different owner pointer");
    check(
        owner->modelId() != changedOwner->modelId(),
        "changed kij changes configured PR model ID");

    CountingVaporChemistry chemistry(registry);
    ChemistrySourceOperator chemistryOp(&chemistry, registry, owner);
    auto oneDimensional = makeOneDimensionalDomain(
        owner, *ownerState.state, registry.size());
    const double chemistryDt = chemistryOp.maxStableDt(
        *oneDimensional, 0.1);
    if (!(std::isfinite(chemistryDt) && chemistryDt > 0.0
          && chemistryDt < 1.0e29))
        std::cerr << "same-owner chemistry dt=" << chemistryDt
                  << " region="
                  << thermo::equilibriumRegionName(ownerState.state->region)
                  << " rho=" << ownerState.state->rho
                  << " p=" << ownerState.state->p
                  << " T=" << ownerState.state->T
                  << " diagnostic=" << ownerState.diagnostics.modelId
                  << '/' << owner->modelId() << '\n';
    check(
        std::isfinite(chemistryDt) && chemistryDt > 0.0
            && chemistryDt < 1.0e29,
        "ChemistrySourceOperator same-owner dt succeeds");
    oneDimensional->beginStep();
    const ConservedVars chemistryInitial = oneDimensional->U(0);
    chemistryOp.apply(*oneDimensional, 1.0e-3 * chemistryDt);
    check(
        oneDimensional->Unew(0).rhoE > chemistryInitial.rhoE,
        "ChemistrySourceOperator same-owner apply couples energy");
    check(
        chemistryOp.hasBoundThermoOwner()
            && chemistryOp.boundThermoModelId() == owner->modelId(),
        "ChemistrySourceOperator exposes frozen owner identity");
    check(
        oneDimensional->commit(),
        "ChemistrySourceOperator same-owner candidate reflashes on commit");

    CountingVaporChemistry sectorChemistry(registry);
    SectorChemistryOperator sectorOp(
        &sectorChemistry, registry, owner, 0.0);
    auto sector = makeSectorDomain(
        owner, *ownerState.state, registry.size());
    const double sectorDt = sectorOp.maxStableDt(*sector, 0.1);
    check(
        std::isfinite(sectorDt) && sectorDt > 0.0 && sectorDt < 1.0e29,
        "SectorChemistryOperator same-owner dt succeeds");
    sector->beginStep();
    sectorOp.apply(*sector, 1.0e-3 * sectorDt);
    check(
        sectorOp.lastActiveCells() == 1 && sectorOp.lastFailures() == 0,
        "SectorChemistryOperator same-owner apply succeeds");
    check(
        sectorOp.hasBoundThermoOwner()
            && sectorOp.boundThermoModelId() == owner->modelId(),
        "SectorChemistryOperator exposes frozen owner identity");
    check(
        sector->commit(),
        "SectorChemistryOperator same-owner candidate reflashes on commit");
    const SectorLedger successfulLedger = captureLedger(sectorOp);
    check(
        successfulLedger.active == 1,
        "Sector ledger is nonzero before mismatch checks");

    auto changedOneDimensional = makeOneDimensionalDomain(
        changedOwner, *changedState.state, registry.size());
    exerciseChemistryMismatch(
        chemistryOp,
        *changedOneDimensional,
        chemistry,
        "changed-kij/model-ID owner");
    auto changedSector = makeSectorDomain(
        changedOwner, *changedState.state, registry.size());
    exerciseSectorMismatch(
        sectorOp,
        *changedSector,
        sectorChemistry,
        successfulLedger,
        "changed-kij/model-ID owner");

    auto cloneOneDimensional = makeOneDimensionalDomain(
        ownerClone, *cloneState.state, registry.size());
    exerciseChemistryMismatch(
        chemistryOp,
        *cloneOneDimensional,
        chemistry,
        "same-ID different-pointer owner");
    auto cloneSector = makeSectorDomain(
        ownerClone, *cloneState.state, registry.size());
    exerciseSectorMismatch(
        sectorOp,
        *cloneSector,
        sectorChemistry,
        successfulLedger,
        "same-ID different-pointer owner");
}

void testSupercriticalOwnerWiring(
    const SpeciesRegistry& registry,
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& interactions)
{
    const auto delegate =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical,
            interactions,
            registry.names(),
            vaporOptions());
    const auto owner =
        std::make_shared<const RegionOverrideThermoModel>(
            delegate, thermo::EquilibriumRegion::SUPERCRITICAL);
    const thermo::EquilibriumResult equilibrium = owner->initializePTY(
        4.0e6, 800.0, vaporComposition(registry));
    check(
        equilibrium.success()
            && equilibrium.state->region
                == thermo::EquilibriumRegion::SUPERCRITICAL,
        "controlled supercritical shared-owner fixture succeeds");
    if (!equilibrium.success()
        || equilibrium.state->region
            != thermo::EquilibriumRegion::SUPERCRITICAL)
        return;

    CountingVaporChemistry chemistry(registry, 0.0);
    ChemistrySourceOperator chemistryOp(&chemistry, registry, owner);
    auto oneDimensional = makeOneDimensionalDomain(
        owner, *equilibrium.state, registry.size());
    const double chemistryDt = chemistryOp.maxStableDt(
        *oneDimensional, 0.1);
    check(
        std::isfinite(chemistryDt) && chemistryDt >= 1.0e29,
        "supercritical inert ChemistrySourceOperator dt succeeds");
    oneDimensional->beginStep();
    chemistryOp.apply(*oneDimensional, 1.0e-8);
    check(
        chemistry.calls() > 0 && oneDimensional->commit()
            && oneDimensional->Q(0).region
                == thermo::EquilibriumRegion::SUPERCRITICAL,
        "supercritical ChemistrySourceOperator apply/commit succeeds");

    CountingVaporChemistry sectorChemistry(registry, 0.0);
    SectorChemistryOperator sectorOp(
        &sectorChemistry, registry, owner, 0.0);
    auto sector = makeSectorDomain(
        owner, *equilibrium.state, registry.size());
    const double sectorDt = sectorOp.maxStableDt(*sector, 0.1);
    check(
        std::isfinite(sectorDt) && sectorDt >= 1.0e29,
        "supercritical inert SectorChemistryOperator dt succeeds");
    sector->beginStep();
    sectorOp.apply(*sector, 1.0e-8);
    check(
        sectorChemistry.calls() > 0
            && sectorOp.lastActiveCells() == 1
            && sectorOp.lastFailures() == 0
            && sector->commit()
            && sector->Q(1, 1, 1).region
                == thermo::EquilibriumRegion::SUPERCRITICAL,
        "supercritical SectorChemistryOperator apply/commit succeeds");
}

enum class DriftKind
{
    MODEL_ID,
    SPECIES,
    ENERGY
};

const char* driftName(DriftKind kind)
{
    switch (kind)
    {
        case DriftKind::MODEL_ID: return "same-pointer model-ID drift";
        case DriftKind::SPECIES: return "same-pointer species drift";
        case DriftKind::ENERGY: return "same-pointer energy drift";
    }
    return "unknown drift";
}

void exerciseIdentityDrift(
    DriftKind kind,
    const SpeciesRegistry& registry,
    const std::shared_ptr<const thermo::EquilibriumThermoModel>& delegate)
{
    auto mutableOwner =
        std::make_shared<MutableIdentityThermoModel>(delegate);
    const std::vector<double> composition = vaporComposition(registry);
    const thermo::EquilibriumResult equilibrium =
        mutableOwner->initializePTY(4.0e6, 800.0, composition);
    check(equilibrium.success(), std::string(driftName(kind)) + " fixture");
    if (!equilibrium.success()) return;

    CountingVaporChemistry chemistry(registry);
    ChemistrySourceOperator chemistryOp(
        &chemistry, registry, mutableOwner);
    auto oneDimensional = makeOneDimensionalDomain(
        mutableOwner, *equilibrium.state, registry.size());
    const double chemistryDt = chemistryOp.maxStableDt(
        *oneDimensional, 0.1);
    oneDimensional->beginStep();
    chemistryOp.apply(*oneDimensional, 1.0e-3 * chemistryDt);
    check(
        oneDimensional->commit(),
        std::string(driftName(kind)) + " pre-drift 1-D success");

    CountingVaporChemistry sectorChemistry(registry);
    SectorChemistryOperator sectorOp(
        &sectorChemistry, registry, mutableOwner, 0.0);
    auto sector = makeSectorDomain(
        mutableOwner, *equilibrium.state, registry.size());
    const double sectorDt = sectorOp.maxStableDt(*sector, 0.1);
    sector->beginStep();
    sectorOp.apply(*sector, 1.0e-3 * sectorDt);
    check(
        sector->commit(),
        std::string(driftName(kind)) + " pre-drift Sector success");
    const SectorLedger successfulLedger = captureLedger(sectorOp);

    switch (kind)
    {
        case DriftKind::MODEL_ID: mutableOwner->driftModelId(); break;
        case DriftKind::SPECIES: mutableOwner->driftSpecies(); break;
        case DriftKind::ENERGY: mutableOwner->driftEnergyConvention(); break;
    }

    exerciseChemistryMismatch(
        chemistryOp,
        *oneDimensional,
        chemistry,
        driftName(kind));
    exerciseSectorMismatch(
        sectorOp,
        *sector,
        sectorChemistry,
        successfulLedger,
        driftName(kind));
}

void testTwoPhaseFailClose(
    const SpeciesRegistry& registry,
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& interactions)
{
    const auto owner =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical, interactions, registry.names(), baseOptions());
    const std::vector<double> overallY =
        twoPhaseComposition(registry, critical, interactions);
    const thermo::EquilibriumResult initial =
        owner->initializePTY(2.0e6, 300.0, overallY);
    check(initial.success(), "two-phase fail-close fixture succeeds");
    if (!initial.success()) return;
    check(
        initial.state->region == thermo::EquilibriumRegion::TWO_PHASE
            && initial.state->liquid.has_value()
            && initial.state->vapor.has_value(),
        "two-phase fail-close fixture carries both phase payloads");
    if (initial.state->region != thermo::EquilibriumRegion::TWO_PHASE
        || !initial.state->vapor.has_value())
        return;

    const thermo::EquilibriumResult directVaporRecovery =
        owner->recoverRhoE(
            initial.state->vapor->rho,
            initial.state->vapor->e,
            initial.state->vapor->massFractions);
    check(
        directVaporRecovery.status
                == thermo::EquilibriumStatus::ACOUSTIC_CLOSURE_FAILED
            && !directVaporRecovery.state.has_value()
            && directVaporRecovery.diagnostics.backendStatus == "SUCCESS"
            && directVaporRecovery.diagnostics.backendPhaseStatus
                == "SUCCESS"
            && directVaporRecovery.diagnostics.backendAcousticStatus
                == "PHASE_TOPOLOGY_CHANGE",
        "same owner rejects saturated vapor on acoustic topology change");

    auto domain = makeSectorDomain(
        owner, *initial.state, registry.size());
    CountingVaporChemistry chemistry(registry);
    SectorChemistryOperator op(&chemistry, registry, owner, 0.0);
    domain->beginStep();
    const ExactDomainState before = captureExactState(*domain);
    const SectorLedger ledgerBefore = captureLedger(op);
    chemistry.resetObservation();
    check(
        rejects([&]() { static_cast<void>(op.maxStableDt(*domain, 0.1)); }),
        "two-phase limiter rejects unadmitted vapor state");
    check(
        chemistry.calls() == 0,
        "two-phase limiter rejects before chemistry RHS");
    check(
        exactStateMatches(*domain, before),
        "two-phase limiter leaves U/Unew/W/Q exact");
    check(
        exactSame(captureLedger(op), ledgerBefore),
        "two-phase limiter preserves Sector ledger");

    check(
        rejects([&]() { op.apply(*domain, 1.0e-8); }),
        "two-phase apply fails closed on the same owner boundary");
    check(
        chemistry.calls() == 0,
        "two-phase apply rejects before chemistry RHS");
    check(
        exactStateMatches(*domain, before),
        "two-phase apply leaves U/Unew/W/Q exact");
    check(
        op.lastActiveCells() == 1
            && op.lastFailures() == 1
            && op.lastStateFailures() == 1
            && op.lastTwoPhaseActiveCells() == 0
            && op.lastPhaseReflashFailures() == 0,
        "two-phase apply records a state rejection without source commit");
}

} // namespace

int main()
{
    try
    {
        const SpeciesRegistry registry;
        const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
        const BinaryInteractionDB interactions =
            completeTestInteractions(registry);

        testExactOwnerAndPointerAdmission(registry, critical, interactions);
        testSupercriticalOwnerWiring(registry, critical, interactions);

        const auto driftDelegate =
            std::make_shared<const thermo::PREquilibriumThermoModel>(
                critical,
                interactions,
                registry.names(),
                vaporOptions());
        exerciseIdentityDrift(
            DriftKind::MODEL_ID, registry, driftDelegate);
        exerciseIdentityDrift(
            DriftKind::SPECIES, registry, driftDelegate);
        exerciseIdentityDrift(
            DriftKind::ENERGY, registry, driftDelegate);

        testTwoPhaseFailClose(registry, critical, interactions);
    }
    catch (const std::exception& error)
    {
        ++failures;
        std::cerr << "FAIL: unexpected exception: " << error.what() << '\n';
    }

    if (failures != 0)
    {
        std::cerr << "CHEMISTRY_OPERATOR_SHARED_THERMO_WIRING_FAIL checks="
                  << checks << " failures=" << failures << '\n';
        return 1;
    }
    std::cout << "CHEMISTRY_OPERATOR_SHARED_THERMO_WIRING_PASS checks="
              << checks << '\n';
    return 0;
}
