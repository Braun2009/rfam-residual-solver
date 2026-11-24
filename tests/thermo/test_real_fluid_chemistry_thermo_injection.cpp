#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "EquilibriumThermoModel.H"
#include "FiniteRate/RealFluidCVODEWrapper.H"
#include "FiniteRateReactor.H"
#include "MixingRules/BinaryInteractionDB.H"
#include "MixingRules/SpeciesCriticalDB.H"
#include "PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "RealFluidState.H"
#include "SpeciesRegistry.H"

using namespace Foam;
using namespace Foam::thermo;

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

class RecordingInertChemistry final : public ChemistryModel
{
public:
    mutable int calls = 0;
    mutable double lastTemperature =
        std::numeric_limits<double>::quiet_NaN();
    mutable double lastPressure =
        std::numeric_limits<double>::quiet_NaN();

    explicit RecordingInertChemistry(const SpeciesRegistry& registry)
    {
        nSpecies_ = 1;
        nReactions_ = 0;
        speciesNames_ = {"O2"};
        molecularWeights_ = {
            registry.molecularWeight(registry.index("O2"))};
    }

    void computeSourceTerms(
        double temperature,
        double pressure,
        double,
        const std::vector<double>&,
        std::vector<double>& omega,
        double& heatRelease) const override
    {
        ++calls;
        lastTemperature = temperature;
        lastPressure = pressure;
        omega.assign(1, 0.0);
        heatRelease = 0.0;
    }
};

class TrackingThermoModel final : public EquilibriumThermoModel
{
    LegacySinglePhaseThermoModel delegate_;
    std::vector<std::string> speciesNames_;
    std::string modelId_;
    EnergyConvention energyConvention_;
    EquilibriumRegion region_;
    std::string diagnosticsModelIdOverride_;
    int twoPhaseFromRecoverCall_ = 0;
    int diagnosticsMismatchFromRecoverCall_ = 0;
    mutable int recoverCalls_ = 0;
    mutable int initializeCalls_ = 0;

    static void convertToClosedTwoPhase(EquilibriumState& state)
    {
        const EquilibriumPhaseProperties phase = *state.singlePhase;
        state.region = EquilibriumRegion::TWO_PHASE;
        state.singlePhase.reset();
        state.liquid = phase;
        state.vapor = phase;
        state.vaporMoleFraction = 0.5;
        state.vaporMassFraction = 0.5;
        state.vaporVolumeFraction = 0.5;
        state.acoustic.closure = AcousticClosure::HOMOGENEOUS_EQUILIBRIUM;
        state.thermodynamics.reset();
    }

    EquilibriumResult route(
        EquilibriumResult result,
        bool forceTwoPhase = false,
        bool forceDiagnosticsMismatch = false) const
    {
        result.diagnostics.modelId = forceDiagnosticsMismatch
            ? "test-wrong-runtime-result-owner"
            : (diagnosticsModelIdOverride_.empty()
                ? modelId_ : diagnosticsModelIdOverride_);
        if (result.state.has_value())
        {
            if (forceTwoPhase)
                convertToClosedTwoPhase(*result.state);
            else
                result.state->region = region_;
        }
        return result;
    }

public:
    TrackingThermoModel(
        std::vector<std::string> speciesNames,
        EquilibriumRegion region,
        std::string modelId = "test-shared-thermo-v1",
        EnergyConvention energyConvention = EnergyConvention::
            SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT)
        : delegate_(speciesNames),
          speciesNames_(std::move(speciesNames)),
          modelId_(std::move(modelId)),
          energyConvention_(energyConvention),
          region_(region)
    {}

    EquilibriumResult recoverRhoE(
        double density,
        double internalEnergy,
        const std::vector<double>& overallMassFractions,
        const EquilibriumSeed* seed = nullptr) const override
    {
        ++recoverCalls_;
        return route(
            delegate_.recoverRhoE(
                density, internalEnergy, overallMassFractions, seed),
            twoPhaseFromRecoverCall_ > 0
                && recoverCalls_ >= twoPhaseFromRecoverCall_,
            diagnosticsMismatchFromRecoverCall_ > 0
                && recoverCalls_ >= diagnosticsMismatchFromRecoverCall_);
    }

    EquilibriumResult initializePTY(
        double pressure,
        double temperature,
        const std::vector<double>& overallMassFractions) const override
    {
        ++initializeCalls_;
        return route(delegate_.initializePTY(
            pressure, temperature, overallMassFractions));
    }

    const std::string& modelId() const override { return modelId_; }
    const std::vector<std::string>& speciesNames() const override
    {
        return speciesNames_;
    }
    EnergyConvention energyConvention() const override
    {
        return energyConvention_;
    }
    bool stateValidForModel(const EquilibriumState& state) const override
    {
        if (state.region == EquilibriumRegion::TWO_PHASE)
            return EquilibriumThermoModel::stateValidForModel(state);
        return delegate_.stateValidForModel(state);
    }
    bool supportsPrimitivePTYInitialization() const override { return true; }

    int recoverCalls() const { return recoverCalls_; }
    int initializeCalls() const { return initializeCalls_; }
    void setModelId(std::string modelId) { modelId_ = std::move(modelId); }
    void setDiagnosticsModelIdOverride(std::string value)
    {
        diagnosticsModelIdOverride_ = std::move(value);
    }
    void setTwoPhaseFromRecoverCall(int call)
    {
        twoPhaseFromRecoverCall_ = call;
    }
    void setDiagnosticsMismatchFromRecoverCall(int call)
    {
        diagnosticsMismatchFromRecoverCall_ = call;
    }
};

template<class Callable>
bool throwsInvalidArgument(Callable&& callable)
{
    try
    {
        callable();
    }
    catch (const std::invalid_argument&)
    {
        return true;
    }
    catch (...)
    {
        return false;
    }
    return false;
}

phaseEquilibrium::RhoEFlashOptions prInjectionTestOptions()
{
    phaseEquilibrium::RhoEFlashOptions options;
    options.minimumPressure = 1.0e5;
    options.maximumPressure = 1.0e7;
    options.minimumTemperature = 250.0;
    options.maximumTemperature = 2000.0;
    options.initialPressure = 4.0e6;
    options.initialTemperature = 800.0;
    options.pressureScanIntervals = 64;
    options.maximumIterations = 120;
    return options;
}

} // namespace

int main()
{
    SpeciesRegistry registry;
    RecordingInertChemistry chemistry(registry);
    SpeciesMapper mapper(registry, chemistry.speciesNames());
    std::vector<double> composition(registry.size(), 0.0);
    composition[registry.index("O2")] = 1.0;

    RealFluidState reference;
    reference.setSpecies(registry.names());
    check(
        reference.updatePTY(1.0e6, 500.0, composition)
            == RealFluidState::Status::SUCCESS,
        "reference rho-e state is available");
    if (failures != 0) return 1;

    RealFluidCVODEWrapper legacy(&chemistry, registry, mapper);
    check(
        !legacy.hasInjectedThermoModel()
            && !legacy.injectedThermoIdentityValid()
            && legacy.boundThermoModelId().empty(),
        "legacy constructor remains an unbound diagnostic path");

    const auto vaporModel = std::make_shared<TrackingThermoModel>(
        registry.names(), EquilibriumRegion::SINGLE_VAPOR);
    RealFluidCVODEWrapper vaporSolver(
        &chemistry, registry, mapper, vaporModel);
    check(
        vaporSolver.hasInjectedThermoModel()
            && vaporSolver.injectedThermoIdentityValid(),
        "injected thermo identity is bound");
    check(
        vaporSolver.boundThermoModelId() == vaporModel->modelId()
            && vaporSolver.boundThermoSpeciesNames() == registry.names()
            && vaporSolver.boundEnergyConvention() == EnergyConvention::
                SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT,
        "model id, ordered species, and energy convention are exact");

    std::vector<double> advancedComposition = composition;
    double advancedEnergy = reference.e();
    double advancedTemperature = reference.T();
    double advancedPressure = reference.p();
    vaporSolver.init(
        0.0, reference.rho(), composition, reference.e(), reference.T());
    const int recoveriesBeforeAdvance = vaporModel->recoverCalls();
    const bool vaporAdvance = vaporSolver.advanceTo(
        1.0e-8, advancedComposition, advancedEnergy,
        advancedTemperature, advancedPressure);
    check(vaporAdvance, "single-vapor injected CVODE path advances");
    check(
        vaporModel->recoverCalls() > recoveriesBeforeAdvance
            && chemistry.calls > 0,
        "RHS and final recovery use the injected thermo model");
    check(
        std::isfinite(chemistry.lastTemperature)
            && std::isfinite(chemistry.lastPressure)
            && chemistry.lastTemperature > 0.0
            && chemistry.lastPressure > 0.0,
        "chemistry receives injected-model temperature and pressure");

    SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    BinaryInteractionDB configuredKij = BinaryInteractionDB::createDefault();
    configuredKij.add("O2", "C3H8", 0.20);
    const auto configuredPrModel =
        std::make_shared<const PREquilibriumThermoModel>(
            critical, configuredKij, registry.names(),
            prInjectionTestOptions());
    std::vector<double> prComposition(registry.size(), 0.0);
    prComposition[registry.index("O2")] = 0.80;
    prComposition[registry.index("C3H8")] = 0.20;
    const EquilibriumResult prInitial = configuredPrModel->initializePTY(
        4.0e6, 800.0, prComposition);
    check(
        prInitial.success()
            && (prInitial.state->region == EquilibriumRegion::SINGLE_VAPOR
                || prInitial.state->region
                    == EquilibriumRegion::SUPERCRITICAL),
        "configured PR model provides an admitted test state");
    if (prInitial.success())
    {
        RealFluidCVODEWrapper configuredPrSolver(
            &chemistry, registry, mapper, configuredPrModel);
        configuredPrSolver.init(
            0.0, prInitial.state->rho, prComposition,
            prInitial.state->e, prInitial.state->T);
        std::vector<double> prAdvancedComposition = prComposition;
        double prAdvancedEnergy = prInitial.state->e;
        double prAdvancedTemperature = prInitial.state->T;
        double prAdvancedPressure = prInitial.state->p;
        const bool prAdvance = configuredPrSolver.advanceTo(
            1.0e-8, prAdvancedComposition, prAdvancedEnergy,
            prAdvancedTemperature, prAdvancedPressure);
        check(
            prAdvance
                && configuredPrSolver.boundThermoModelId()
                    == configuredPrModel->modelId()
                && std::abs(prAdvancedTemperature-prInitial.state->T)
                    <= 1.0e-8*prInitial.state->T
                && std::abs(prAdvancedPressure-prInitial.state->p)
                    <= 1.0e-8*prInitial.state->p,
            "CVODE uses the exact injected configured-PR identity");
    }

    const auto supercriticalModel = std::make_shared<TrackingThermoModel>(
        registry.names(), EquilibriumRegion::SUPERCRITICAL,
        "test-shared-thermo-supercritical-v1");
    RealFluidCVODEWrapper supercriticalSolver(
        &chemistry, registry, mapper, supercriticalModel);
    check(
        !throwsInvalidArgument([&]
        {
            supercriticalSolver.init(
                0.0, reference.rho(), composition,
                reference.e(), reference.T());
        }),
        "supercritical injected state is admitted");

    for (const EquilibriumRegion rejectedRegion : {
             EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE,
             EquilibriumRegion::SINGLE_LIQUID})
    {
        const auto rejectedModel = std::make_shared<TrackingThermoModel>(
            registry.names(), rejectedRegion,
            "test-rejected-region-" +
                std::to_string(static_cast<int>(rejectedRegion)));
        RealFluidCVODEWrapper rejectedSolver(
            &chemistry, registry, mapper, rejectedModel);
        check(
            throwsInvalidArgument([&]
            {
                rejectedSolver.init(
                    0.0, reference.rho(), composition,
                    reference.e(), reference.T());
            }),
            "non-admitted injected region fails before CVODE RHS");
    }

    check(
        throwsInvalidArgument([&]
        {
            RealFluidCVODEWrapper invalid(
                &chemistry, registry, mapper,
                std::shared_ptr<const EquilibriumThermoModel>());
        }),
        "null injected thermo model is rejected");
    check(
        throwsInvalidArgument([&]
        {
            auto emptyIdentity = std::make_shared<TrackingThermoModel>(
                registry.names(), EquilibriumRegion::SINGLE_VAPOR, "");
            RealFluidCVODEWrapper invalid(
                &chemistry, registry, mapper, emptyIdentity);
        }),
        "empty model identity is rejected");
    check(
        throwsInvalidArgument([&]
        {
            auto blankIdentity = std::make_shared<TrackingThermoModel>(
                registry.names(), EquilibriumRegion::SINGLE_VAPOR,
                " \t\r\n");
            RealFluidCVODEWrapper invalid(
                &chemistry, registry, mapper, blankIdentity);
        }),
        "blank model identity is rejected");
    check(
        throwsInvalidArgument([&]
        {
            std::vector<std::string> wrongOrder = registry.names();
            std::swap(wrongOrder[0], wrongOrder[1]);
            auto wrongSpecies = std::make_shared<TrackingThermoModel>(
                wrongOrder, EquilibriumRegion::SINGLE_VAPOR,
                "test-wrong-species-v1");
            RealFluidCVODEWrapper invalid(
                &chemistry, registry, mapper, wrongSpecies);
        }),
        "ordered species identity mismatch is rejected");
    check(
        throwsInvalidArgument([&]
        {
            auto wrongEnergy = std::make_shared<TrackingThermoModel>(
                registry.names(), EquilibriumRegion::SINGLE_VAPOR,
                "test-wrong-energy-v1",
                static_cast<EnergyConvention>(91));
            RealFluidCVODEWrapper invalid(
                &chemistry, registry, mapper, wrongEnergy);
        }),
        "energy convention mismatch is rejected");

    const auto diagnosticsMismatch = std::make_shared<TrackingThermoModel>(
        registry.names(), EquilibriumRegion::SINGLE_VAPOR,
        "test-diagnostics-binding-v1");
    diagnosticsMismatch->setDiagnosticsModelIdOverride("wrong-result-owner");
    RealFluidCVODEWrapper diagnosticsMismatchSolver(
        &chemistry, registry, mapper, diagnosticsMismatch);
    check(
        throwsInvalidArgument([&]
        {
            diagnosticsMismatchSolver.init(
                0.0, reference.rho(), composition,
                reference.e(), reference.T());
        }),
        "recovery diagnostics model id mismatch is rejected");

    const auto driftingModel = std::make_shared<TrackingThermoModel>(
        registry.names(), EquilibriumRegion::SINGLE_VAPOR,
        "test-drift-v1");
    RealFluidCVODEWrapper driftingSolver(
        &chemistry, registry, mapper, driftingModel);
    driftingSolver.init(
        0.0, reference.rho(), composition, reference.e(), reference.T());
    const int driftRecoveries = driftingModel->recoverCalls();
    driftingModel->setModelId("test-drift-v2");
    std::vector<double> driftComposition = composition;
    double driftEnergy = reference.e();
    double driftTemperature = reference.T();
    double driftPressure = reference.p();
    const bool driftAdvance = driftingSolver.advanceTo(
        1.0e-8, driftComposition, driftEnergy,
        driftTemperature, driftPressure);
    check(
        !driftAdvance && !driftingSolver.injectedThermoIdentityValid(),
        "runtime thermo identity drift fails closed before chemistry RHS");
    check(
        driftingModel->recoverCalls() == driftRecoveries
            && driftComposition == composition
            && driftEnergy == reference.e()
            && driftTemperature == reference.T()
            && driftPressure == reference.p(),
        "identity-drift failure performs no recovery or output publication");

    const auto reactorModel = std::make_shared<TrackingThermoModel>(
        registry.names(), EquilibriumRegion::SINGLE_VAPOR,
        "test-reactor-shared-thermo-v1");
    FiniteRateReactor reactor(&chemistry, registry, reactorModel);
    const double estimatedDt = reactor.estimateConservativeMaxDt(
        reference.rho(), reference.p(), reference.T(), composition, 0.1);
    check(
        reactor.hasInjectedThermoModel()
            && reactor.boundThermoModelId() == reactorModel->modelId()
            && reactorModel->initializeCalls() > 0
            && estimatedDt >= 1.0e29,
        "reactor timestep route uses the same injected thermo identity");

    std::vector<double> reactorComposition = composition;
    double reactorEnergyIncrement = 123.0;
    double reactorTemperature = -1.0;
    double reactorPressure = -1.0;
    const int reactorRecoveriesBefore = reactorModel->recoverCalls();
    const bool zeroStep = reactor.reactConservative(
        reference.rho(), reference.e(), reactorComposition, 0.0,
        reactorEnergyIncrement, reactorTemperature, reactorPressure);
    check(
        zeroStep && reactorEnergyIncrement == 0.0
            && reactorModel->recoverCalls() > reactorRecoveriesBefore
            && reactorTemperature > 0.0 && reactorPressure > 0.0,
        "reactor initial recovery uses shared model at zero step");

    const auto liquidModel = std::make_shared<TrackingThermoModel>(
        registry.names(), EquilibriumRegion::SINGLE_LIQUID,
        "test-reactor-liquid-rejection-v1");
    FiniteRateReactor liquidReactor(&chemistry, registry, liquidModel);
    const int chemistryCallsBeforeLiquidEstimate = chemistry.calls;
    check(
        liquidReactor.estimateConservativeMaxDt(
            reference.rho(), reference.p(), reference.T(),
            composition, 0.1) == 0.0
            && chemistry.calls == chemistryCallsBeforeLiquidEstimate,
        "timestep route rejects liquid state before chemistry evaluation");
    std::vector<double> rejectedComposition = composition;
    double rejectedIncrement = 456.0;
    double rejectedTemperature = 777.0;
    double rejectedPressure = 888.0;
    const bool liquidAccepted = liquidReactor.reactConservative(
        reference.rho(), reference.e(), rejectedComposition, 0.0,
        rejectedIncrement, rejectedTemperature, rejectedPressure);
    check(
        !liquidAccepted && liquidReactor.lastFailureCode() == 2
            && rejectedComposition == composition
            && rejectedIncrement == 0.0
            && rejectedTemperature == 777.0
            && rejectedPressure == 888.0,
        "reactor rejects non-vapor model state without publishing outputs");

    const auto closedTwoPhaseModel = std::make_shared<TrackingThermoModel>(
        registry.names(), EquilibriumRegion::SINGLE_VAPOR,
        "test-closed-two-phase-result-v1");
    closedTwoPhaseModel->setTwoPhaseFromRecoverCall(1);
    const EquilibriumResult closedTwoPhase =
        closedTwoPhaseModel->recoverRhoE(
            reference.rho(), reference.e(), composition);
    check(
        closedTwoPhase.success()
            && closedTwoPhase.state->region == EquilibriumRegion::TWO_PHASE,
        "forced mid-RHS two-phase result is structurally valid");

    const auto checkPersistentMidRhsFailure = [&]
    (
        const std::shared_ptr<TrackingThermoModel>& model,
        const std::string& label
    )
    {
        FiniteRateReactor failingReactor(&chemistry, registry, model);
        std::vector<double> candidateComposition = composition;
        const std::vector<double> beforeComposition = candidateComposition;
        double candidateIncrement = 987.0;
        double candidateTemperature = reference.T();
        double candidatePressure = reference.p();
        const double beforeTemperature = candidateTemperature;
        const double beforePressure = candidatePressure;
        const int chemistryCallsBefore = chemistry.calls;
        const bool accepted = failingReactor.reactConservative(
            reference.rho(), reference.e(), candidateComposition, 1.0e-8,
            candidateIncrement, candidateTemperature, candidatePressure);
        check(
            !accepted
                && failingReactor.lastFailureCode() == 2
                && failingReactor.lastSubcycles() == -1
                && model->recoverCalls() > 3,
            label + " remains rejected through all timestep refinements");
        check(
            candidateComposition == beforeComposition
                && candidateIncrement == 0.0
                && candidateTemperature == beforeTemperature
                && candidatePressure == beforePressure
                && chemistry.calls == chemistryCallsBefore,
            label + " publishes no species, energy, T, p, or RHS source");
    };

    const auto midRhsTwoPhaseModel =
        std::make_shared<TrackingThermoModel>(
            registry.names(), EquilibriumRegion::SINGLE_VAPOR,
            "test-mid-rhs-two-phase-v1");
    midRhsTwoPhaseModel->setTwoPhaseFromRecoverCall(3);
    checkPersistentMidRhsFailure(
        midRhsTwoPhaseModel, "mid-RHS TWO_PHASE recovery");

    const auto midRhsDiagnosticsDriftModel =
        std::make_shared<TrackingThermoModel>(
            registry.names(), EquilibriumRegion::SINGLE_VAPOR,
            "test-mid-rhs-diagnostics-v1");
    midRhsDiagnosticsDriftModel->setDiagnosticsMismatchFromRecoverCall(3);
    checkPersistentMidRhsFailure(
        midRhsDiagnosticsDriftModel,
        "mid-RHS diagnostics modelId drift");

    if (failures == 0)
    {
        std::cout << "REAL_FLUID_CHEMISTRY_THERMO_INJECTION_PASS checks="
                  << checks << '\n';
        return 0;
    }
    std::cout << "REAL_FLUID_CHEMISTRY_THERMO_INJECTION_FAIL failures="
              << failures << " checks=" << checks << '\n';
    return 1;
}
