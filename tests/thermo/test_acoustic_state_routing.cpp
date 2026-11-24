#include "Domain.H"
#include "FluxScheme.H"
#include "NonReflectingOutlet.H"
#include "Operator.H"
#include "PointImplicit.H"
#include "PseudoTransient.H"
#include "SectorConvectionOperator.H"
#include "PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "PhaseEquilibrium/PRTPFlash.H"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Foam;
using namespace Foam::phaseEquilibrium;

namespace
{

int failures = 0;

void check(bool condition, const std::string& label)
{
    std::cout << "  [" << (condition ? "PASS" : "FAIL") << "] "
              << label << '\n';
    if (!condition) ++failures;
}

bool close(double lhs, double rhs, double relativeTolerance = 1.0e-11)
{
    return std::isfinite(lhs) && std::isfinite(rhs)
        && std::abs(lhs - rhs) <= relativeTolerance
            *std::max({std::abs(lhs),std::abs(rhs),1.0});
}

bool same(const ConservedVars& lhs, const ConservedVars& rhs)
{
    return lhs.rho == rhs.rho
        && lhs.rhoU == rhs.rhoU
        && lhs.rhoE == rhs.rhoE
        && lhs.rhoY == rhs.rhoY
        && lhs.rhoZ == rhs.rhoZ
        && lhs.rhoC == rhs.rhoC;
}

struct DomainSnapshot
{
    std::vector<ConservedVars> committed;
    std::vector<ConservedVars> staged;
    std::uint64_t acceptedSteps = 0;
    std::uint64_t rejectedSteps = 0;
};

DomainSnapshot snapshot(const Domain& domain)
{
    DomainSnapshot result;
    for (int cell=0; cell<domain.N(); ++cell)
    {
        result.committed.push_back(domain.U(cell));
        result.staged.push_back(domain.Unew(cell));
    }
    result.acceptedSteps=domain.acceptedSteps();
    result.rejectedSteps=domain.rejectedSteps();
    return result;
}

bool same(const DomainSnapshot& before, const Domain& after)
{
    if (before.committed.size()!=static_cast<std::size_t>(after.N())
        || before.staged.size()!=static_cast<std::size_t>(after.N())
        || before.acceptedSteps!=after.acceptedSteps()
        || before.rejectedSteps!=after.rejectedSteps())
        return false;
    for (int cell=0; cell<after.N(); ++cell)
        if (!same(before.committed[cell],after.U(cell))
            || !same(before.staged[cell],after.Unew(cell)))
            return false;
    return true;
}

class PrescribedAcousticModel final : public thermo::EquilibriumThermoModel
{
    std::vector<std::string> species_ = {"O2"};
    thermo::LegacySinglePhaseThermoModel delegate_{species_};
    double derivative_;
    bool primitiveInitialization_;
    thermo::AcousticClosure closure_;
    std::string id_;

    thermo::EquilibriumResult prescribe(
        thermo::EquilibriumResult result) const
    {
        if (!result.success()) return result;
        result.state->acoustic.closure = closure_;
        result.state->acoustic.pressureDensityDerivative = derivative_;
        result.state->thermodynamics.reset();
        result.diagnostics.modelId = id_;
        if (!result.state->valid())
        {
            result.status = thermo::EquilibriumStatus::ACOUSTIC_CLOSURE_FAILED;
            result.state.reset();
        }
        return result;
    }

public:
    PrescribedAcousticModel(
        double soundSpeed,
        bool primitiveInitialization,
        thermo::AcousticClosure closure=thermo::AcousticClosure::UNKNOWN)
    :
        derivative_(soundSpeed*soundSpeed),
        primitiveInitialization_(primitiveInitialization),
        closure_(closure==thermo::AcousticClosure::UNKNOWN
            ? (primitiveInitialization
                ? thermo::AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE
                : thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM)
            : closure),
        id_("prescribed-acoustic-test:" + std::to_string(soundSpeed)
            + (primitiveInitialization ? ":pty" : ":rhoe-only"))
    {}

    thermo::EquilibriumResult recoverRhoE(
        double density,
        double internalEnergy,
        const std::vector<double>& massFractions,
        const thermo::EquilibriumSeed* seed) const override
    {
        return prescribe(delegate_.recoverRhoE(
            density,internalEnergy,massFractions,seed));
    }

    thermo::EquilibriumResult initializePTY(
        double pressure,
        double temperature,
        const std::vector<double>& massFractions) const override
    {
        if (!primitiveInitialization_)
            return thermo::EquilibriumThermoModel::initializePTY(
                pressure,temperature,massFractions);
        return prescribe(delegate_.initializePTY(
            pressure,temperature,massFractions));
    }

    const std::string& modelId() const override { return id_; }
    const std::vector<std::string>& speciesNames() const override
    { return species_; }
    EnergyConvention energyConvention() const override
    {
        return EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT;
    }
    bool supportsPrimitivePTYInitialization() const override
    { return primitiveInitialization_; }
};

class RecordingFlux final : public FluxScheme
{
    mutable std::mutex mutex_;

public:
    mutable int primitiveCalls = 0;
    mutable int cachedCalls = 0;
    mutable std::vector<double> soundSpeeds;

    void computeFlux(
        const PrimitiveVars& left,
        const PrimitiveVars&,
        const std::vector<double>&,
        ConservedVars& flux) const override
    {
        const std::lock_guard<std::mutex> lock(mutex_);
        ++primitiveCalls;
        flux = ConservedVars(static_cast<int>(left.Y.size()));
    }

    void computeFluxCached(
        const PrimitiveVars& left,
        const PrimitiveVars&,
        const std::vector<double>&,
        ConservedVars& flux,
        double,
        double soundLeft,
        double,
        double,
        double soundRight,
        double) const override
    {
        const std::lock_guard<std::mutex> lock(mutex_);
        ++cachedCalls;
        soundSpeeds.push_back(soundLeft);
        soundSpeeds.push_back(soundRight);
        flux = ConservedVars(static_cast<int>(left.Y.size()));
    }

    std::string name() const override { return "recording equilibrium flux"; }
};

class MomentumKick final : public Operator
{
    double source_;

public:
    explicit MomentumKick(double source) : source_(source) {}

    void apply(Domain& domain, double dt) const override
    {
        for (int cell = 0; cell < domain.N(); ++cell)
            domain.Unew(cell).rhoU[0] += dt*source_;
    }

    std::string name() const override { return "momentum kick"; }
};

class DtRecorder final : public Operator
{
public:
    mutable double appliedDt = std::numeric_limits<double>::quiet_NaN();

    void apply(Domain&, double dt) const override { appliedDt = dt; }
    std::string name() const override { return "dt recorder"; }
};

class TransactionProbe final : public Operator
{
public:
    mutable int transactionBegins = 0;
    mutable int applies = 0;
    mutable int stableDtQueries = 0;

    void beginStepTransaction(const Domain&) const override
    { ++transactionBegins; }
    void apply(Domain&, double) const override { ++applies; }
    double maxStableDt(const Domain&, double) const override
    {
        ++stableDtQueries;
        return 1.0;
    }
    std::string name() const override { return "transaction probe"; }
};

void routeSyntheticHemForDiagnostic(
    const FluxScheme& flux,
    const PrimitiveVars& left,
    const thermo::EquilibriumState& leftThermo,
    const PrimitiveVars& right,
    const thermo::EquilibriumState& rightThermo,
    ConservedVars& result)
{
    if (leftThermo.acoustic.closure
            !=thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM
        || rightThermo.acoustic.closure
            !=thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM
        || !leftThermo.acoustic.valid() || !rightThermo.acoustic.valid())
        throw std::runtime_error(
            "diagnostic HEM route requires valid synthetic HEM acoustics");
    flux.computeFluxCached(
        left,right,{1.0,0.0,0.0},result,
        leftThermo.rho,leftThermo.acoustic.soundSpeed(),leftThermo.h,
        rightThermo.rho,rightThermo.acoustic.soundSpeed(),rightThermo.h);
}

PrimitiveVars primitive(double axialVelocity)
{
    PrimitiveVars result(1);
    result.p = 4.0e6;
    result.T = 300.0;
    result.U = {axialVelocity,0.0,0.0};
    result.Y = {1.0};
    result.Z = 0.0;
    result.C = 0.0;
    return result;
}

thermo::EquilibriumState syntheticTwoPhaseHemState(double soundSpeed)
{
    constexpr double pressure=1.0e6;
    constexpr double vaporFraction=0.25;
    constexpr double temperature=200.0;
    constexpr double molecularWeight=0.0319988;
    constexpr double gasConstant=8.314462618;

    thermo::EquilibriumPhaseProperties liquid;
    liquid.massFractions={1.0};
    liquid.moleFractions={1.0};
    liquid.molecularWeight=molecularWeight;
    liquid.rho=1000.0;
    liquid.h=100000.0;
    liquid.e=liquid.h-pressure/liquid.rho;
    liquid.compressibility=
        pressure*molecularWeight/(liquid.rho*gasConstant*temperature);

    thermo::EquilibriumPhaseProperties vapor;
    vapor.massFractions={1.0};
    vapor.moleFractions={1.0};
    vapor.molecularWeight=molecularWeight;
    vapor.rho=10.0;
    vapor.h=500000.0;
    vapor.e=vapor.h-pressure/vapor.rho;
    vapor.compressibility=
        pressure*molecularWeight/(vapor.rho*gasConstant*temperature);

    const double specificVolume=
        (1.0-vaporFraction)/liquid.rho+vaporFraction/vapor.rho;

    thermo::EquilibriumState state;
    state.region=thermo::EquilibriumRegion::TWO_PHASE;
    state.overallMassFractions={1.0};
    state.p=pressure;
    state.T=temperature;
    state.rho=1.0/specificVolume;
    state.h=(1.0-vaporFraction)*liquid.h+vaporFraction*vapor.h;
    state.e=(1.0-vaporFraction)*liquid.e+vaporFraction*vapor.e;
    state.vaporMoleFraction=vaporFraction;
    state.vaporMassFraction=vaporFraction;
    state.vaporVolumeFraction=
        (vaporFraction/vapor.rho)/specificVolume;
    state.acoustic.closure=
        thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM;
    state.acoustic.pressureDensityDerivative=soundSpeed*soundSpeed;
    state.liquid=liquid;
    state.vapor=vapor;
    return state;
}

bool initializeWithEquilibrium(
    Domain& domain,
    const thermo::EquilibriumState& equilibrium,
    double axialVelocity=0.0)
{
    PrimitiveVars state(static_cast<int>(equilibrium.overallMassFractions.size()));
    state.p=equilibrium.p;
    state.T=equilibrium.T;
    state.U={axialVelocity,0.0,0.0};
    state.Y=equilibrium.overallMassFractions;
    state.Z=0.0;
    state.C=0.0;
    return domain.replaceCompleteField(
        std::vector<PrimitiveVars>(domain.N(),state),
        std::vector<thermo::EquilibriumState>(domain.N(),equilibrium));
}

RhoEFlashOptions realPrRoutingOptions()
{
    RhoEFlashOptions options;
    options.minimumPressure=1.0e6;
    options.maximumPressure=6.0e6;
    options.minimumTemperature=280.0;
    options.maximumTemperature=600.0;
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

bool initialize(Domain& domain, double axialVelocity)
{
    const PrimitiveVars state = primitive(axialVelocity);
    for (int cell = 0; cell < domain.N(); ++cell)
        domain.stagedW(cell) = state;
    return domain.primitiveToConserved();
}

} // namespace

int main()
{
    std::cout << "=== Equilibrium acoustic routing contract ===\n";
    constexpr double prescribedSoundSpeed = 123.0;
    const auto model = std::make_shared<const PrescribedAcousticModel>(
        prescribedSoundSpeed,true);

    Domain convectionDomain(3,0.3,1,model);
    check(initialize(convectionDomain,10.0),
        "prescribed-acoustic Domain initializes");
    check(close(
            convectionDomain.Q(0).acoustic.soundSpeed(),
            prescribedSoundSpeed),
        "Domain Q retains the model acoustic payload");

    std::unique_ptr<RecordingFlux> recordingFlux(new RecordingFlux());
    RecordingFlux* recording = recordingFlux.get();
    ConvectionOperator convection(std::move(recordingFlux));
    convectionDomain.beginStep();
    convection.apply(convectionDomain,1.0e-7);
    bool routedSound = recording->cachedCalls > 0
        && recording->primitiveCalls == 0;
    for (double sound : recording->soundSpeeds)
        routedSound = routedSound && close(sound,prescribedSoundSpeed);
    check(routedSound,
        "first-order AUSM consumes only the committed Q acoustics");

    const double cfl = 0.4;
    const double expectedDt = cfl*convectionDomain.dx()
        /(10.0 + prescribedSoundSpeed);
    check(close(convection.maxStableDt(convectionDomain,cfl),expectedDt),
        "convection CFL consumes Q sound speed");

    AUSMFluxAdapter ausm;
    thermo::EquilibriumState invalidThermo = convectionDomain.Q(0);
    invalidThermo.acoustic.pressureDensityDerivative =
        std::numeric_limits<double>::quiet_NaN();
    ConservedVars untouched(1);
    untouched.rho = 42.0;
    bool invalidRejected = false;
    try
    {
        ausm.computeFluxFromEquilibrium(
            convectionDomain.W(0),invalidThermo,
            convectionDomain.W(1),convectionDomain.Q(1),
            {1.0,0.0,0.0},untouched);
    }
    catch (const std::runtime_error&)
    {
        invalidRejected = true;
    }
    check(invalidRejected && untouched.rho == 42.0,
        "invalid cached acoustics fail before flux mutation");

    Domain implicitDomain(3,0.3,1,model);
    check(initialize(implicitDomain,0.0),
        "point-implicit Domain initializes");
    constexpr double source = 1000.0;
    constexpr double implicitDt = 1.0e-6;
    std::vector<std::unique_ptr<Operator>> implicitOperators;
    implicitOperators.emplace_back(new MomentumKick(source));
    PointImplicitIntegrator pointImplicit;
    const double damping = 1.0/(1.0
        + implicitDt*prescribedSoundSpeed/implicitDomain.dx());
    check(pointImplicit.step(
            implicitDomain,implicitOperators,implicitDt),
        "point-implicit step accepts a Q-acoustic state");
    check(close(
            implicitDomain.U(1).rhoU[0],
            implicitDt*source*damping,2.0e-10),
        "point-implicit spectral radius consumes Q sound speed");

    Domain pseudoDomain(3,0.3,1,model);
    check(initialize(pseudoDomain,10.0),
        "pseudo-transient Domain initializes");
    std::unique_ptr<DtRecorder> dtRecorder(new DtRecorder());
    DtRecorder* recorder = dtRecorder.get();
    std::vector<std::unique_ptr<Operator>> pseudoOperators;
    pseudoOperators.push_back(std::move(dtRecorder));
    PseudoTransientIntegrator pseudo(cfl);
    check(pseudo.step(pseudoDomain,pseudoOperators,1.0),
        "pseudo-transient fallback step succeeds");
    check(close(recorder->appliedDt,expectedDt),
        "pseudo-transient timestep consumes Q sound speed");

    Domain outletDomain(3,0.3,1,model);
    outletDomain.stagedW(0) = primitive(10.0);
    outletDomain.stagedW(1) = primitive(10.0);
    outletDomain.stagedW(2) = primitive(20.0);
    check(outletDomain.primitiveToConserved(),
        "nonreflecting outlet Domain initializes");
    constexpr double outletDt = 1.0e-6;
    const double outletRho = outletDomain.Q(2).rho;
    const double outletVelocityGradient = 10.0/outletDomain.dx();
    const double outgoingWave = (20.0 + prescribedSoundSpeed)
        *outletRho*prescribedSoundSpeed*outletVelocityGradient;
    const double expectedOutletPressure = 4.0e6
        -0.5*outletDt*outgoingWave;
    NonReflectingOutlet outlet(4.0e6,0.3,{"O2"},0.0);
    outlet.apply(outletDomain,outletDt);
    check(close(outletDomain.W(2).p,expectedOutletPressure,2.0e-10),
        "nonreflecting characteristic update consumes committed Q acoustics");

    std::vector<ConservedVars> field;
    std::vector<double> seeds;
    for (int cell = 0; cell < convectionDomain.N(); ++cell)
    {
        field.push_back(convectionDomain.U(cell));
        seeds.push_back(convectionDomain.W(cell).T);
    }
    const auto rhoEOnlyModel =
        std::make_shared<const PrescribedAcousticModel>(
            prescribedSoundSpeed,false);
    Domain phaseAwareDomain(3,0.3,1,rhoEOnlyModel);
    check(phaseAwareDomain.replaceConservedField(field,seeds),
        "rho-e-only model recovers a complete Q field");
    check(!solverConsumerAdmitsAcousticState(phaseAwareDomain.Q(0)),
        "production solver predicate rejects synthetic HEM acoustics");

    RecordingFlux diagnosticFlux;
    ConservedVars diagnosticResult(1);
    routeSyntheticHemForDiagnostic(
        diagnosticFlux,
        phaseAwareDomain.W(0),phaseAwareDomain.Q(0),
        phaseAwareDomain.W(1),phaseAwareDomain.Q(1),
        diagnosticResult);
    check(diagnosticFlux.cachedCalls==1
            && diagnosticFlux.primitiveCalls==0,
        "test-only diagnostic harness can inspect synthetic HEM acoustics");

    RecordingFlux directHemFlux;
    ConservedVars directHemUntouched(1);
    directHemUntouched.rho=43.0;
    bool directHemRejected=false;
    try
    {
        directHemFlux.computeFluxFromEquilibrium(
            phaseAwareDomain.W(0),phaseAwareDomain.Q(0),
            phaseAwareDomain.W(1),phaseAwareDomain.Q(1),
            {1.0,0.0,0.0},directHemUntouched);
    }
    catch(const std::runtime_error&)
    {
        directHemRejected=true;
    }
    check(directHemRejected && directHemFlux.cachedCalls==0
            && directHemFlux.primitiveCalls==0
            && directHemUntouched.rho==43.0,
        "direct production flux rejects HEM before flux mutation");

    std::unique_ptr<RecordingFlux> hemConvectionFlux(new RecordingFlux());
    RecordingFlux* hemConvectionRecording=hemConvectionFlux.get();
    ConvectionOperator hemConvection(std::move(hemConvectionFlux));
    const DomainSnapshot hemConvectionBefore=snapshot(phaseAwareDomain);
    bool hemConvectionApplyRejected=false;
    try
    {
        hemConvection.apply(phaseAwareDomain,1.0e-7);
    }
    catch(const std::runtime_error&)
    {
        hemConvectionApplyRejected=true;
    }
    bool hemConvectionCflRejected=false;
    try
    {
        (void)hemConvection.maxStableDt(phaseAwareDomain,cfl);
    }
    catch(const std::runtime_error&)
    {
        hemConvectionCflRejected=true;
    }
    check(hemConvectionApplyRejected && hemConvectionCflRejected
            && hemConvectionRecording->cachedCalls==0
            && hemConvectionRecording->primitiveCalls==0
            && same(hemConvectionBefore,phaseAwareDomain),
        "generic convection apply and CFL reject HEM before mutation");

    const thermo::EquilibriumState twoPhaseHem=
        syntheticTwoPhaseHemState(prescribedSoundSpeed);
    check(twoPhaseHem.valid()
            && !solverConsumerAdmitsAcousticState(twoPhaseHem),
        "production solver predicate rejects a valid two-phase HEM state");
    Domain twoPhaseConvectionDomain(3,0.3,1,rhoEOnlyModel);
    check(initializeWithEquilibrium(
            twoPhaseConvectionDomain,twoPhaseHem,10.0),
        "two-phase HEM convection Domain initializes");
    RecordingFlux directTwoPhaseFlux;
    ConservedVars directTwoPhaseUntouched(1);
    directTwoPhaseUntouched.rho=44.0;
    bool directTwoPhaseRejected=false;
    try
    {
        directTwoPhaseFlux.computeFluxFromEquilibrium(
            twoPhaseConvectionDomain.W(0),twoPhaseConvectionDomain.Q(0),
            twoPhaseConvectionDomain.W(1),twoPhaseConvectionDomain.Q(1),
            {1.0,0.0,0.0},directTwoPhaseUntouched);
    }
    catch(const std::runtime_error&)
    {
        directTwoPhaseRejected=true;
    }
    std::unique_ptr<RecordingFlux> twoPhaseConvectionFlux(new RecordingFlux());
    RecordingFlux* twoPhaseConvectionRecording=twoPhaseConvectionFlux.get();
    ConvectionOperator twoPhaseConvection(std::move(twoPhaseConvectionFlux));
    const DomainSnapshot twoPhaseConvectionBefore=
        snapshot(twoPhaseConvectionDomain);
    bool twoPhaseConvectionApplyRejected=false;
    try
    {
        twoPhaseConvection.apply(twoPhaseConvectionDomain,1.0e-7);
    }
    catch(const std::runtime_error&)
    {
        twoPhaseConvectionApplyRejected=true;
    }
    bool twoPhaseConvectionCflRejected=false;
    try
    {
        (void)twoPhaseConvection.maxStableDt(twoPhaseConvectionDomain,cfl);
    }
    catch(const std::runtime_error&)
    {
        twoPhaseConvectionCflRejected=true;
    }
    check(directTwoPhaseRejected
            && directTwoPhaseFlux.cachedCalls==0
            && directTwoPhaseFlux.primitiveCalls==0
            && directTwoPhaseUntouched.rho==44.0
            && twoPhaseConvectionApplyRejected
            && twoPhaseConvectionCflRejected
            && twoPhaseConvectionRecording->cachedCalls==0
            && twoPhaseConvectionRecording->primitiveCalls==0
            && same(twoPhaseConvectionBefore,twoPhaseConvectionDomain),
        "direct and generic consumers reject two-phase HEM before mutation");

    Domain hemImplicitDomain(3,0.3,1,rhoEOnlyModel);
    check(hemImplicitDomain.replaceConservedField(field,seeds),
        "HEM point-implicit Domain initializes");
    std::unique_ptr<TransactionProbe> implicitProbe(new TransactionProbe());
    TransactionProbe* implicitProbeView=implicitProbe.get();
    std::vector<std::unique_ptr<Operator>> hemImplicitOperators;
    hemImplicitOperators.push_back(std::move(implicitProbe));
    const DomainSnapshot hemImplicitBefore=snapshot(hemImplicitDomain);
    check(!pointImplicit.step(
            hemImplicitDomain,hemImplicitOperators,implicitDt)
            && implicitProbeView->transactionBegins==0
            && implicitProbeView->applies==0
            && implicitProbeView->stableDtQueries==0
            && same(hemImplicitBefore,hemImplicitDomain),
        "point-implicit rejects HEM before StepTransaction or Unew mutation");

    Domain hemPseudoDomain(3,0.3,1,rhoEOnlyModel);
    check(hemPseudoDomain.replaceConservedField(field,seeds),
        "HEM pseudo-transient Domain initializes");
    std::unique_ptr<TransactionProbe> pseudoProbe(new TransactionProbe());
    TransactionProbe* pseudoProbeView=pseudoProbe.get();
    std::vector<std::unique_ptr<Operator>> hemPseudoOperators;
    hemPseudoOperators.push_back(std::move(pseudoProbe));
    const DomainSnapshot hemPseudoBefore=snapshot(hemPseudoDomain);
    check(!pseudo.step(hemPseudoDomain,hemPseudoOperators,1.0)
            && pseudoProbeView->transactionBegins==0
            && pseudoProbeView->applies==0
            && pseudoProbeView->stableDtQueries==0
            && same(hemPseudoBefore,hemPseudoDomain),
        "pseudo-transient rejects HEM before transaction or dt cache use");

    Domain twoPhaseImplicitDomain(3,0.3,1,rhoEOnlyModel);
    check(initializeWithEquilibrium(twoPhaseImplicitDomain,twoPhaseHem),
        "two-phase HEM point-implicit Domain initializes");
    std::unique_ptr<TransactionProbe> twoPhaseImplicitProbe(
        new TransactionProbe());
    TransactionProbe* twoPhaseImplicitProbeView=twoPhaseImplicitProbe.get();
    std::vector<std::unique_ptr<Operator>> twoPhaseImplicitOperators;
    twoPhaseImplicitOperators.push_back(std::move(twoPhaseImplicitProbe));
    const DomainSnapshot twoPhaseImplicitBefore=snapshot(twoPhaseImplicitDomain);
    check(!pointImplicit.step(
            twoPhaseImplicitDomain,twoPhaseImplicitOperators,implicitDt)
            && twoPhaseImplicitProbeView->transactionBegins==0
            && twoPhaseImplicitProbeView->applies==0
            && twoPhaseImplicitProbeView->stableDtQueries==0
            && same(twoPhaseImplicitBefore,twoPhaseImplicitDomain),
        "point-implicit rejects two-phase HEM before StepTransaction mutation");

    Domain twoPhasePseudoDomain(3,0.3,1,rhoEOnlyModel);
    check(initializeWithEquilibrium(twoPhasePseudoDomain,twoPhaseHem),
        "two-phase HEM pseudo-transient Domain initializes");
    std::unique_ptr<TransactionProbe> twoPhasePseudoProbe(
        new TransactionProbe());
    TransactionProbe* twoPhasePseudoProbeView=twoPhasePseudoProbe.get();
    std::vector<std::unique_ptr<Operator>> twoPhasePseudoOperators;
    twoPhasePseudoOperators.push_back(std::move(twoPhasePseudoProbe));
    const DomainSnapshot twoPhasePseudoBefore=snapshot(twoPhasePseudoDomain);
    check(!pseudo.step(twoPhasePseudoDomain,twoPhasePseudoOperators,1.0)
            && twoPhasePseudoProbeView->transactionBegins==0
            && twoPhasePseudoProbeView->applies==0
            && twoPhasePseudoProbeView->stableDtQueries==0
            && same(twoPhasePseudoBefore,twoPhasePseudoDomain),
        "pseudo-transient rejects two-phase HEM before transaction or dt cache use");

    constexpr int sectorCells=4;
    SectorMeshMetrics mesh(
        vortexProfilePath(),sectorCells,sectorCells,sectorCells,120.0);
    Domain sectorSinglePhase(
        sectorCells,sectorCells,sectorCells,
        mesh.length(),0.065,mesh.sectorAngle(),1,model);
    check(initialize(sectorSinglePhase,10.0),
        "single-phase sector acoustic Domain initializes");
    std::unique_ptr<RecordingFlux> sectorSinglePhaseFlux(new RecordingFlux());
    RecordingFlux* sectorSinglePhaseRecording=sectorSinglePhaseFlux.get();
    SectorConvectionOperator sectorSinglePhaseConvection(
        std::move(sectorSinglePhaseFlux),mesh);
    const double sectorSinglePhaseDt=
        sectorSinglePhaseConvection.maxStableDt(sectorSinglePhase,cfl);
    sectorSinglePhase.beginStep();
    sectorSinglePhaseConvection.apply(sectorSinglePhase,1.0e-7);
    bool sectorSinglePhaseRouted=
        std::isfinite(sectorSinglePhaseDt) && sectorSinglePhaseDt>0.0
        && sectorSinglePhaseRecording->cachedCalls>0
        && sectorSinglePhaseRecording->primitiveCalls==0;
    for(double sound:sectorSinglePhaseRecording->soundSpeeds)
        sectorSinglePhaseRouted=
            sectorSinglePhaseRouted && close(sound,prescribedSoundSpeed);
    check(sectorSinglePhaseRouted,
        "sector convection apply and CFL retain legacy single-phase routing");
    std::vector<ConservedVars> sectorField;
    std::vector<double> sectorSeeds;
    for (int cell=0; cell<sectorSinglePhase.N(); ++cell)
    {
        sectorField.push_back(sectorSinglePhase.U(cell));
        sectorSeeds.push_back(sectorSinglePhase.W(cell).T);
    }
    Domain sectorHem(
        sectorCells,sectorCells,sectorCells,
        mesh.length(),0.065,mesh.sectorAngle(),1,rhoEOnlyModel);
    check(sectorHem.replaceConservedField(sectorField,sectorSeeds),
        "synthetic HEM sector Domain initializes");
    std::unique_ptr<RecordingFlux> sectorFlux(new RecordingFlux());
    RecordingFlux* sectorRecording=sectorFlux.get();
    SectorConvectionOperator sectorConvection(std::move(sectorFlux),mesh);
    const DomainSnapshot sectorBefore=snapshot(sectorHem);
    bool sectorApplyRejected=false;
    try
    {
        sectorConvection.apply(sectorHem,1.0e-7);
    }
    catch(const std::runtime_error&)
    {
        sectorApplyRejected=true;
    }
    bool sectorCflRejected=false;
    try
    {
        (void)sectorConvection.maxStableDt(sectorHem,cfl);
    }
    catch(const std::runtime_error&)
    {
        sectorCflRejected=true;
    }
    check(sectorApplyRejected && sectorCflRejected
            && sectorRecording->cachedCalls==0
            && sectorRecording->primitiveCalls==0
            && same(sectorBefore,sectorHem),
        "sector convection apply and CFL reject HEM before cache or flux use");

    Domain sectorTwoPhase(
        sectorCells,sectorCells,sectorCells,
        mesh.length(),0.065,mesh.sectorAngle(),1,rhoEOnlyModel);
    check(initializeWithEquilibrium(sectorTwoPhase,twoPhaseHem,10.0),
        "two-phase HEM sector Domain initializes");
    std::unique_ptr<RecordingFlux> sectorTwoPhaseFlux(new RecordingFlux());
    RecordingFlux* sectorTwoPhaseRecording=sectorTwoPhaseFlux.get();
    SectorConvectionOperator sectorTwoPhaseConvection(
        std::move(sectorTwoPhaseFlux),mesh);
    const DomainSnapshot sectorTwoPhaseBefore=snapshot(sectorTwoPhase);
    bool sectorTwoPhaseApplyRejected=false;
    try
    {
        sectorTwoPhaseConvection.apply(sectorTwoPhase,1.0e-7);
    }
    catch(const std::runtime_error&)
    {
        sectorTwoPhaseApplyRejected=true;
    }
    bool sectorTwoPhaseCflRejected=false;
    try
    {
        (void)sectorTwoPhaseConvection.maxStableDt(sectorTwoPhase,cfl);
    }
    catch(const std::runtime_error&)
    {
        sectorTwoPhaseCflRejected=true;
    }
    check(sectorTwoPhaseApplyRejected && sectorTwoPhaseCflRejected
            && sectorTwoPhaseRecording->cachedCalls==0
            && sectorTwoPhaseRecording->primitiveCalls==0
            && same(sectorTwoPhaseBefore,sectorTwoPhase),
        "sector convection rejects two-phase HEM before cache or flux use");

    const SpeciesCriticalDB realCritical=SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB realKij=BinaryInteractionDB::createDefault();
    const std::vector<std::string> realSpecies={"O2","C3H8"};
    const RhoEFlashOptions realOptions=realPrRoutingOptions();
    PRTPFlash realTp(realCritical,realKij,realSpecies);
    const TPFlashResult realSplit=
        realTp.flashMoleFractions(2.0e6,300.0,{0.2,0.8});
    const auto realModel=
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            realCritical,realKij,realSpecies,realOptions);
    const thermo::EquilibriumResult realTwoPhase=
        realSplit.phaseSplitReady()
        ? realModel->initializePTY(
            2.0e6,300.0,realSplit.feedMassFractions)
        : thermo::EquilibriumResult();
    const thermo::EquilibriumResult realSinglePhase=
        realModel->initializePTY(4.0e6,450.0,{0.0,1.0});
    check(realTwoPhase.success()
            && realTwoPhase.state->region==thermo::EquilibriumRegion::TWO_PHASE
            && realTwoPhase.state->acoustic.evidence
                ==thermo::AcousticEvidence::CONVERGED_RHOE_ISENTROPE
            && (ModelCapabilities::phaseAwareSoundSpeed
                ? solverConsumerAdmitsAcousticState(*realTwoPhase.state)
                : !solverConsumerAdmitsAcousticState(*realTwoPhase.state)),
        ModelCapabilities::phaseAwareSoundSpeed
            ? "real two-phase PR rho-e-Y acoustics are solver-admitted"
            : "real two-phase PR rho-e-Y acoustics remain fail-closed");
    check(realSinglePhase.success()
            && realSinglePhase.state->region!=thermo::EquilibriumRegion::TWO_PHASE
            && realSinglePhase.state->acoustic.evidence
                ==thermo::AcousticEvidence::CONVERGED_RHOE_ISENTROPE
            && (ModelCapabilities::phaseAwareSoundSpeed
                ? solverConsumerAdmitsAcousticState(*realSinglePhase.state)
                : !solverConsumerAdmitsAcousticState(*realSinglePhase.state)),
        ModelCapabilities::phaseAwareSoundSpeed
            ? "real single-phase PR rho-e-Y acoustics are solver-admitted"
            : "real single-phase PR rho-e-Y acoustics remain fail-closed");

    if(ModelCapabilities::phaseAwareSoundSpeed
        && realTwoPhase.success() && realSinglePhase.success())
    {
        Domain realSingleDomain(3,0.3,2,realModel);
        check(initializeWithEquilibrium(
                realSingleDomain,*realSinglePhase.state,10.0),
            "real single-phase PR routing Domain initializes");
        RecordingFlux realSingleFlux;
        ConservedVars realSingleResult(2);
        realSingleFlux.computeFluxFromEquilibrium(
            realSingleDomain.W(0),realSingleDomain.Q(0),
            realSingleDomain.W(1),realSingleDomain.Q(1),
            {1.0,0.0,0.0},realSingleResult);
        check(realSingleFlux.cachedCalls==1
                && realSingleFlux.primitiveCalls==0,
            "direct flux accepts verified single-phase PR HEM acoustics");

        Domain realGenericDomain(3,0.3,2,realModel);
        check(initializeWithEquilibrium(
                realGenericDomain,*realTwoPhase.state,10.0),
            "real two-phase PR generic Domain initializes");
        std::unique_ptr<RecordingFlux> realGenericFlux(new RecordingFlux());
        RecordingFlux* realGenericRecording=realGenericFlux.get();
        ConvectionOperator realGenericConvection(std::move(realGenericFlux));
        const DomainSnapshot realGenericBefore=snapshot(realGenericDomain);
        const double realGenericDt=
            realGenericConvection.maxStableDt(realGenericDomain,cfl);
        realGenericDomain.beginStep();
        realGenericConvection.apply(realGenericDomain,1.0e-7);
        check(std::isfinite(realGenericDt) && realGenericDt>0.0
                && realGenericRecording->cachedCalls>0
                && realGenericRecording->primitiveCalls==0
                && same(realGenericBefore,realGenericDomain),
            "generic apply and CFL route verified two-phase PR acoustics");

        Domain realImplicitDomain(3,0.3,2,realModel);
        check(initializeWithEquilibrium(
                realImplicitDomain,*realTwoPhase.state),
            "real two-phase PR implicit Domain initializes");
        std::unique_ptr<TransactionProbe> realImplicitProbe(
            new TransactionProbe());
        TransactionProbe* realImplicitProbeView=realImplicitProbe.get();
        std::vector<std::unique_ptr<Operator>> realImplicitOperators;
        realImplicitOperators.push_back(std::move(realImplicitProbe));
        check(pointImplicit.step(
                realImplicitDomain,realImplicitOperators,implicitDt)
                && realImplicitProbeView->transactionBegins==1
                && realImplicitProbeView->applies==1
                && realImplicitProbeView->stableDtQueries==0
                && realImplicitDomain.acceptedSteps()==1
                && realImplicitDomain.Q(0).acoustic.evidence
                    ==thermo::AcousticEvidence::CONVERGED_RHOE_ISENTROPE,
            "point-implicit routes verified two-phase PR acoustics transactionally");

        Domain realPseudoDomain(3,0.3,2,realModel);
        check(initializeWithEquilibrium(
                realPseudoDomain,*realTwoPhase.state),
            "real two-phase PR pseudo-transient Domain initializes");
        std::unique_ptr<TransactionProbe> realPseudoProbe(
            new TransactionProbe());
        TransactionProbe* realPseudoProbeView=realPseudoProbe.get();
        std::vector<std::unique_ptr<Operator>> realPseudoOperators;
        realPseudoOperators.push_back(std::move(realPseudoProbe));
        check(pseudo.step(realPseudoDomain,realPseudoOperators,1.0)
                && realPseudoProbeView->transactionBegins==1
                && realPseudoProbeView->applies==1
                && realPseudoProbeView->stableDtQueries==1
                && realPseudoDomain.acceptedSteps()==1
                && realPseudoDomain.Q(0).acoustic.evidence
                    ==thermo::AcousticEvidence::CONVERGED_RHOE_ISENTROPE,
            "pseudo-transient routes verified two-phase PR acoustics transactionally");

        Domain realSectorDomain(
            sectorCells,sectorCells,sectorCells,
            mesh.length(),0.065,mesh.sectorAngle(),2,realModel);
        check(initializeWithEquilibrium(
                realSectorDomain,*realTwoPhase.state,10.0),
            "real two-phase PR sector Domain initializes");
        std::unique_ptr<RecordingFlux> realSectorFlux(new RecordingFlux());
        RecordingFlux* realSectorRecording=realSectorFlux.get();
        SectorConvectionOperator realSectorConvection(
            std::move(realSectorFlux),mesh);
        const DomainSnapshot realSectorBefore=snapshot(realSectorDomain);
        const double realSectorDt=
            realSectorConvection.maxStableDt(realSectorDomain,cfl);
        realSectorDomain.beginStep();
        realSectorConvection.apply(realSectorDomain,1.0e-7);
        check(std::isfinite(realSectorDt) && realSectorDt>0.0
                && realSectorRecording->cachedCalls>0
                && realSectorRecording->primitiveCalls==0
                && same(realSectorBefore,realSectorDomain),
            "sector apply and CFL route verified two-phase PR acoustics");
    }

    NonReflectingOutlet phaseAwareOutlet(4.0e6,0.3,{"O2"});
    bool phaseAwareOutletRejected = false;
    try
    {
        phaseAwareOutlet.apply(phaseAwareDomain,1.0e-7);
    }
    catch (const std::runtime_error&)
    {
        phaseAwareOutletRejected = true;
    }
    bool phaseAwareOutletUnchanged = !phaseAwareDomain.tripped();
    for (int cell = 0; cell < phaseAwareDomain.N(); ++cell)
        phaseAwareOutletUnchanged = phaseAwareOutletUnchanged
            && same(phaseAwareDomain.U(cell),field[cell])
            && close(
                phaseAwareDomain.Q(cell).acoustic.soundSpeed(),
                prescribedSoundSpeed);
    check(phaseAwareOutletRejected && phaseAwareOutletUnchanged,
        "phase-aware characteristic outlet rejects before Domain mutation");
    const auto rhoEOnlyFrozenModel=
        std::make_shared<const PrescribedAcousticModel>(
            prescribedSoundSpeed,false,
            thermo::AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE);
    Domain rhoEOnlyFrozenDomain(3,0.3,1,rhoEOnlyFrozenModel);
    check(rhoEOnlyFrozenDomain.replaceConservedField(field,seeds),
        "rho-e-only frozen-acoustic Domain initializes");
    std::unique_ptr<RecordingFlux> musclFlux(new RecordingFlux());
    ConvectionOperator phaseAwareMuscl(
        std::move(musclFlux),false,true);
    bool musclRejected = false;
    try
    {
        phaseAwareMuscl.apply(rhoEOnlyFrozenDomain,1.0e-7);
    }
    catch (const std::runtime_error&)
    {
        musclRejected = true;
    }
    check(musclRejected,
        "phase-aware MUSCL rejects underdetermined PTY face states");

    Domain multidimensionalMusclDomain(3,3,0.3,0.3,1,model);
    check(initialize(multidimensionalMusclDomain,10.0),
        "multidimensional MUSCL test Domain initializes");
    std::unique_ptr<RecordingFlux> multidimensionalFlux(new RecordingFlux());
    RecordingFlux* multidimensionalRecording=multidimensionalFlux.get();
    ConvectionOperator multidimensionalMuscl(
        std::move(multidimensionalFlux),false,true);
    bool multidimensionalMusclRejected=false;
    try
    {
        multidimensionalMuscl.apply(multidimensionalMusclDomain,1.0e-7);
    }
    catch(const std::runtime_error&)
    {
        multidimensionalMusclRejected=true;
    }
    check(multidimensionalMusclRejected
            && multidimensionalRecording->primitiveCalls==0
            && multidimensionalRecording->cachedCalls==0,
        "multidimensional MUSCL rejects instead of silently using first order");

    if (failures != 0)
    {
        std::cerr << "ACOUSTIC_STATE_ROUTING_FAIL failures="
                  << failures << '\n';
        return 1;
    }
    std::cout << "ACOUSTIC_STATE_ROUTING_PASS\n";
    return 0;
}
