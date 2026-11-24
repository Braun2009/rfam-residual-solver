#include "ModelCapabilities.H"
#include "Operator.H"
#include "PhaseEquilibrium/ConservedPhaseMixture.H"
#include "PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "PhaseEquilibrium/PRPhaseCaloric.H"
#include "PhaseEquilibrium/PRTPFlash.H"
#include "SolverBase.H"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef RFAM_USE_OPENMP
#include <omp.h>
#endif

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

bool close(double lhs, double rhs, double relativeTolerance = 1.0e-8)
{
    return std::isfinite(lhs) && std::isfinite(rhs)
        && std::abs(lhs - rhs) <= relativeTolerance
            *std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

bool sameConserved(const ConservedVars& lhs, const ConservedVars& rhs)
{
    if (!close(lhs.rho, rhs.rho) || !close(lhs.rhoE, rhs.rhoE)
        || lhs.rhoU.size() != rhs.rhoU.size()
        || lhs.rhoY.size() != rhs.rhoY.size())
        return false;
    for (std::size_t direction = 0; direction < lhs.rhoU.size(); ++direction)
        if (!close(lhs.rhoU[direction], rhs.rhoU[direction])) return false;
    for (std::size_t species = 0; species < lhs.rhoY.size(); ++species)
        if (!close(lhs.rhoY[species], rhs.rhoY[species])) return false;
    return close(lhs.rhoZ, rhs.rhoZ) && close(lhs.rhoC, rhs.rhoC);
}

struct ReferenceState
{
    bool valid = false;
    double pressure = 2.0e6;
    double temperature = 300.0;
    double density = std::numeric_limits<double>::quiet_NaN();
    double energy = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> massFractions;
};

ConservedPhaseProperties conservedPhase(const PRPhaseCaloricProperties& phase)
{
    ConservedPhaseProperties result;
    result.molecularWeight = phase.molecularWeight;
    result.rho = phase.rho;
    result.h = phase.h;
    result.e = phase.e;
    result.s = phase.s;
    result.massFractions = phase.massFractions;
    return result;
}

ReferenceState makeReference(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij)
{
    ReferenceState reference;
    const std::vector<std::string> species = {"O2", "C3H8"};
    PRTPFlash tp(critical, kij, species);
    const TPFlashResult split = tp.flashMoleFractions(
        reference.pressure, reference.temperature, {0.2, 0.8});
    if (!split.valid() || split.region != EquilibriumRegion::TWO_PHASE
        || !split.phaseSplitReady())
        return reference;

    PRPhaseCaloric caloric(critical, kij, species);
    const PRPhaseCaloricResult liquid = caloric.evaluate(
        reference.pressure, reference.temperature,
        split.liquidMoleFractions, PhaseRoot::LIQUID);
    const PRPhaseCaloricResult vapor = caloric.evaluate(
        reference.pressure, reference.temperature,
        split.vaporMoleFractions, PhaseRoot::VAPOR);
    if (!liquid.success() || !vapor.success()) return reference;

    reference.massFractions = split.feedMassFractions;
    const ConservedMixtureResult mixture = assembleConservedPhaseMixture(
        reference.pressure, split.vaporMoleFraction,
        conservedPhase(*liquid.state), conservedPhase(*vapor.state),
        reference.massFractions);
    if (!mixture.success()) return reference;

    reference.density = mixture.state->rho;
    reference.energy = mixture.state->e;
    reference.valid = true;
    return reference;
}

RhoEFlashOptions flashOptions(const ReferenceState& reference)
{
    RhoEFlashOptions options;
    options.minimumPressure = 1.0e6;
    options.maximumPressure = 3.0e6;
    options.minimumTemperature = 280.0;
    options.maximumTemperature = 320.0;
    options.initialPressure = reference.pressure;
    options.initialTemperature = reference.temperature;
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

ConservedVars makeConserved(const ReferenceState& reference)
{
    constexpr double velocity = 8.0;
    ConservedVars result(2);
    result.rho = reference.density;
    result.rhoU[0] = reference.density*velocity;
    result.rhoU[1] = 0.0;
    result.rhoU[2] = 0.0;
    result.rhoE = reference.density*(reference.energy
        + 0.5*velocity*velocity);
    for (std::size_t species = 0;
         species < reference.massFractions.size(); ++species)
        result.rhoY[species] = reference.density*reference.massFractions[species];
    result.rhoZ = 0.0;
    result.rhoC = 0.0;
    return result;
}

struct RunResult
{
    bool initialized = false;
    bool stepped = false;
    bool tripped = false;
    double maximumStableDt = std::numeric_limits<double>::quiet_NaN();
    double time = 0.0;
    std::uint64_t acceptedSteps = 0;
    std::uint64_t rejectedSteps = 0;
    std::vector<ConservedVars> field;
    std::vector<thermo::EquilibriumState> states;
};

RunResult runUniformDiagnostic(
    std::shared_ptr<const thermo::PREquilibriumThermoModel> model,
    const ReferenceState& reference,
    const ConservedVars& initial)
{
    RunResult result;
    Domain domain(6, 0.06, 2, std::move(model));
    const std::vector<ConservedVars> field(domain.N(), initial);
    const std::vector<double> seeds(domain.N(), reference.temperature);
    result.initialized = domain.replaceConservedField(field, seeds);
    if (!result.initialized) return result;

    std::unique_ptr<FluxScheme> flux(
        new AUSMFluxAdapter(domain.speciesNames()));
    std::unique_ptr<Operator> convection(
        new ConvectionOperator(std::move(flux)));
    result.maximumStableDt = convection->maxStableDt(domain, 0.1);
    SolverBase solver(domain);
    solver.addOperator(std::move(convection));
    result.stepped = solver.step(1.0e-8);
    result.tripped = domain.tripped();
    result.time = solver.time();
    result.acceptedSteps = domain.acceptedSteps();
    result.rejectedSteps = domain.rejectedSteps();
    for (int cell = 0; cell < domain.N(); ++cell)
    {
        result.field.push_back(domain.U(cell));
        result.states.push_back(domain.Q(cell));
    }
    return result;
}

bool validUniformState(
    const RunResult& result,
    const ConservedVars& initial,
    const thermo::PREquilibriumThermoModel& model)
{
    if (result.field.size() != result.states.size() || result.field.empty())
        return false;
    for (std::size_t cell = 0; cell < result.field.size(); ++cell)
    {
        const thermo::EquilibriumState& state = result.states[cell];
        if (!sameConserved(result.field[cell], initial)
            || state.region != thermo::EquilibriumRegion::TWO_PHASE
            || state.acoustic.closure
                != thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM
            || !state.acoustic.valid() || !model.stateValidForModel(state))
            return false;
    }
    return true;
}

bool sameRun(const RunResult& lhs, const RunResult& rhs)
{
    if (lhs.field.size() != rhs.field.size()
        || lhs.states.size() != rhs.states.size()
        || !close(lhs.maximumStableDt, rhs.maximumStableDt)
        || !close(lhs.time, rhs.time)
        || lhs.acceptedSteps != rhs.acceptedSteps
        || lhs.rejectedSteps != rhs.rejectedSteps)
        return false;
    for (std::size_t cell = 0; cell < lhs.field.size(); ++cell)
    {
        if (!sameConserved(lhs.field[cell], rhs.field[cell])
            || !close(lhs.states[cell].p, rhs.states[cell].p)
            || !close(lhs.states[cell].T, rhs.states[cell].T)
            || !close(lhs.states[cell].rho, rhs.states[cell].rho)
            || !close(lhs.states[cell].e, rhs.states[cell].e)
            || !close(lhs.states[cell].acoustic.soundSpeed(),
                rhs.states[cell].acoustic.soundSpeed()))
            return false;
    }
    return true;
}

} // namespace

int main()
{
    std::cout << "=== PR HEM first-order short-step diagnostic ===\n";
    check(!ModelCapabilities::twoPhaseFlash
            && !ModelCapabilities::phaseAwareEnergyInversion
            && !ModelCapabilities::latentHeatClosure
            && ModelCapabilities::phaseAwareSoundSpeed
            && !ModelCapabilities::phaseAwareBoundaryStates
            && !ModelCapabilities::phaseAwareCharacteristicBoundary
            && !ModelCapabilities::denseCryogenicTransport
            && !ModelCapabilities::phaseAwareChemistry,
        "only phase-aware sound speed is enabled; complete production remains closed");

    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB kij = BinaryInteractionDB::createDefault();
    const ReferenceState reference = makeReference(critical, kij);
    check(reference.valid, "two-phase O2/C3H8 diagnostic reference is available");
    if (!reference.valid) return 1;

    const auto model = std::make_shared<const thermo::PREquilibriumThermoModel>(
        critical, kij, std::vector<std::string>{"O2", "C3H8"},
        flashOptions(reference));
    check(model->supportsPrimitivePTYInitialization(),
        "diagnostic model exposes complete PTY initialization");

    const ConservedVars initial = makeConserved(reference);
#ifdef RFAM_USE_OPENMP
    const int originalThreadCount = omp_get_max_threads();
    omp_set_dynamic(0);
    omp_set_num_threads(1);
#endif
    const RunResult serial = runUniformDiagnostic(model, reference, initial);
#ifdef RFAM_USE_OPENMP
    omp_set_num_threads(std::min(4, omp_get_num_procs()));
#endif
    const RunResult parallel = runUniformDiagnostic(model, reference, initial);
#ifdef RFAM_USE_OPENMP
    omp_set_num_threads(originalThreadCount);
#endif

    check(serial.initialized && serial.stepped && !serial.tripped,
        "PR model completes one accepted HEM solver-consumer step");
    check(std::isfinite(serial.maximumStableDt)
            && serial.maximumStableDt > 1.0e-8,
        "production CFL evaluates the committed verified HEM acoustic payload");
    check(serial.time == 1.0e-8 && serial.acceptedSteps == 1
            && serial.rejectedSteps == 0,
        "accepted diagnostic step advances only its local transaction");
    check(validUniformState(serial, initial, *model),
        "uniform two-phase U/W/Q field remains model-valid after commit");
    check(parallel.initialized && parallel.stepped && !parallel.tripped
            && sameRun(serial, parallel),
        "serial and OpenMP short diagnostics are deterministic");

    if (failures != 0)
    {
        std::cerr << "PR_EQUILIBRIUM_SHORT_STEP_FAIL failures="
                  << failures << '\n';
        return 1;
    }
    std::cout << "PR_EQUILIBRIUM_SHORT_STEP_PASS\n"
              << "NOTE: first-order uniform phase diagnostic only; it does not "
              << "authorize boundary, transport, chemistry, or sector advance.\n";
    return 0;
}
