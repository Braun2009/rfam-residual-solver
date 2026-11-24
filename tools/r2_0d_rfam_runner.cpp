/*---------------------------------------------------------------------------*\
    RFAM R2 0-D trajectory runner.

    This is an offline implementation diagnostic for the current TwoStepCO
    candidate.  It advances an ideal-gas, constant-pressure reactor using
    species mass fractions plus sensible mixture enthalpy.  It does not select
    a detailed mechanism, fit kinetic parameters, or validate R2 physics.
\*---------------------------------------------------------------------------*/

#include "FiniteRate/TwoStepCO.H"
#include "NASA7Thermo.H"
#include "SpeciesRegistry.H"

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sundials/sundials_types.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr double universalGasConstant = 8.314462618; // J/(mol K)
constexpr double minimumTemperature = 200.0;
constexpr double maximumTemperature = 8000.0;

struct SolverTolerances
{
    const char* profile = "protocol-v1";
    double relative = 1.0e-9;
    double speciesAbsolute = 1.0e-13;
    double enthalpyAbsolute = 1.0e-5;
};

SolverTolerances solverTolerances(bool strictEquivalence)
{
    if (!strictEquivalence) return SolverTolerances{};
    return SolverTolerances{
        "strict-equivalence-v1",
        1.0e-12,
        1.0e-16,
        1.0e-8};
}

double parseFinite(const char* text, const std::string& label)
{
    std::size_t consumed = 0;
    double value = 0.0;
    try
    {
        value = std::stod(text, &consumed);
    }
    catch (const std::exception&)
    {
        throw std::invalid_argument("invalid-number:" + label);
    }
    if (consumed != std::string(text).size() || !std::isfinite(value))
        throw std::invalid_argument("invalid-number:" + label);
    return value;
}

std::size_t parseCount(const char* text)
{
    std::size_t consumed = 0;
    unsigned long value = 0;
    try
    {
        value = std::stoul(text, &consumed);
    }
    catch (const std::exception&)
    {
        throw std::invalid_argument("invalid-sample-count");
    }
    if (consumed != std::string(text).size() || value < 2 || value > 1000)
        throw std::invalid_argument("invalid-sample-count");
    return static_cast<std::size_t>(value);
}

bool validCaseId(const std::string& value)
{
    if (value.empty() || value.size() > 80) return false;
    for (const unsigned char character : value)
    {
        if ((character >= 'a' && character <= 'z')
            || (character >= 'A' && character <= 'Z')
            || (character >= '0' && character <= '9')
            || character == '_' || character == '-')
            continue;
        return false;
    }
    return true;
}

struct Snapshot
{
    double time = 0.0;
    double temperature = 0.0;
    double pressure = 0.0;
    double density = 0.0;
    double sensibleEnthalpy = 0.0;
    double totalEnthalpy = 0.0;
    double heatRelease = 0.0;
    double temperatureRate = 0.0;
    std::vector<double> massFractions;
    std::vector<double> massFractionRates;
};

class ConstantPressureSensibleEnthalpyReactor
{
    static constexpr int speciesCount = 6;
    static constexpr int equationCount = speciesCount + 1;
    static constexpr int rootCount = 2;

    Foam::TwoStepCO chemistry_;
    Foam::SpeciesRegistry registry_;
    const Foam::NASA7Thermo& thermo_;
    std::vector<std::string> speciesNames_;
    double pressure_ = 0.0;
    double maxInternalStep_ = 0.0;
    SolverTolerances tolerances_;
    double ignitionTemperature_ = 0.0;
    bool ignitionDetected_ = false;
    double ignitionDelay_ = -1.0;
    double currentTime_ = 0.0;
    bool fuelDepleted_ = false;
    int depletionEventCount_ = 0;
    int intervalBisectionCount_ = 0;
    int roundoffClampCount_ = 0;

    SUNContext context_ = nullptr;
    void* memory_ = nullptr;
    N_Vector state_ = nullptr;
    N_Vector absoluteTolerance_ = nullptr;
    SUNMatrix matrix_ = nullptr;
    SUNLinearSolver linearSolver_ = nullptr;

    std::vector<double> rhsMassFractions_;
    std::vector<double> rhsSources_;

    static int rhs(
        sunrealtype,
        N_Vector state,
        N_Vector derivative,
        void* userData)
    {
        auto* self = static_cast<ConstantPressureSensibleEnthalpyReactor*>(
            userData);
        try
        {
            const sunrealtype* values = N_VGetArrayPointer(state);
            sunrealtype* rates = N_VGetArrayPointer(derivative);
            self->evaluationComposition(values, self->rhsMassFractions_);
            const double sensibleEnthalpy = values[speciesCount];
            const double temperature = self->recoverTemperature(
                self->rhsMassFractions_, sensibleEnthalpy);
            const double density = self->density(
                self->rhsMassFractions_, temperature);
            double heatRelease = 0.0;
            self->chemistry_.computeSourceTerms(
                temperature,
                self->pressure_,
                density,
                self->rhsMassFractions_,
                self->rhsSources_,
                heatRelease);
            if (!std::isfinite(heatRelease)) return -1;
            for (int species = 0; species < speciesCount; ++species)
            {
                rates[species] = self->rhsSources_[species] / density;
                if (!std::isfinite(rates[species])) return -1;
            }
            rates[speciesCount] = heatRelease / density;
            return std::isfinite(rates[speciesCount]) ? 0 : -1;
        }
        catch (const std::exception&)
        {
            return -1;
        }
    }

    static int root(
        sunrealtype,
        N_Vector state,
        sunrealtype* roots,
        void* userData)
    {
        auto* self = static_cast<ConstantPressureSensibleEnthalpyReactor*>(
            userData);
        try
        {
            const sunrealtype* values = N_VGetArrayPointer(state);
            self->evaluationComposition(values, self->rhsMassFractions_);
            roots[0] = self->recoverTemperature(
                self->rhsMassFractions_, values[speciesCount])
                - self->ignitionTemperature_;
            roots[1] = self->fuelDepleted_ ? 1.0 : values[0];
            return std::isfinite(roots[0]) && std::isfinite(roots[1])
                ? 0 : -1;
        }
        catch (const std::exception&)
        {
            return -1;
        }
    }

    void evaluationComposition(
        const sunrealtype* state,
        std::vector<double>& composition) const
    {
        composition.resize(speciesCount);
        double sum = 0.0;
        for (int species = 0; species < speciesCount; ++species)
        {
            if (!std::isfinite(state[species]))
                throw std::runtime_error("nonfinite-composition-state");
            // CVODE can evaluate a Newton trial just outside the physical
            // composition domain as a reactant is depleted.  The production
            // wrappers use the same RHS projection.  The protocol validator
            // independently rejects material negative drift in accepted
            // samples.
            composition[species] = std::max(0.0, double(state[species]));
            sum += composition[species];
        }
        if (!std::isfinite(sum) || sum <= 0.0)
            throw std::runtime_error("empty-composition-state");
    }

    double sensibleEnthalpy(
        const std::vector<double>& composition,
        double temperature) const
    {
        return thermo_.hsens_mass_mix(
            speciesNames_, composition, temperature);
    }

    double recoverTemperature(
        const std::vector<double>& composition,
        double targetEnthalpy) const
    {
        if (!std::isfinite(targetEnthalpy))
            throw std::runtime_error("nonfinite-sensible-enthalpy");
        double lower = minimumTemperature;
        double upper = maximumTemperature;
        double lowerResidual = sensibleEnthalpy(composition, lower)
            - targetEnthalpy;
        double upperResidual = sensibleEnthalpy(composition, upper)
            - targetEnthalpy;
        if (!std::isfinite(lowerResidual) || !std::isfinite(upperResidual)
            || lowerResidual > 0.0 || upperResidual < 0.0)
            throw std::runtime_error("temperature-recovery-out-of-envelope");
        for (int iteration = 0; iteration < 100; ++iteration)
        {
            const double middle = 0.5 * (lower + upper);
            const double residual = sensibleEnthalpy(composition, middle)
                - targetEnthalpy;
            if (!std::isfinite(residual))
                throw std::runtime_error("temperature-recovery-nonfinite");
            if (residual < 0.0)
                lower = middle;
            else
                upper = middle;
        }
        const double result = 0.5 * (lower + upper);
        if (!std::isfinite(result))
            throw std::runtime_error("temperature-recovery-failed");
        return result;
    }

    double density(
        const std::vector<double>& composition,
        double temperature) const
    {
        double inverseMolecularWeight = 0.0;
        for (int species = 0; species < speciesCount; ++species)
            inverseMolecularWeight += composition[species]
                / chemistry_.molecularWeights()[species];
        const double gasConstant = universalGasConstant
            * inverseMolecularWeight;
        const double result = pressure_ / (gasConstant * temperature);
        if (!std::isfinite(result) || result <= 0.0)
            throw std::runtime_error("density-recovery-failed");
        return result;
    }

    double totalEnthalpy(
        const std::vector<double>& composition,
        double sensible) const
    {
        double formation = 0.0;
        for (int species = 0; species < speciesCount; ++species)
        {
            const int global = registry_.index(speciesNames_[species]);
            formation += composition[species]
                * registry_.formationEnthalpy(global)
                / registry_.molecularWeight(global);
        }
        return sensible + formation;
    }

    void destroy()
    {
        if (memory_) CVodeFree(&memory_);
        if (linearSolver_) SUNLinSolFree(linearSolver_);
        if (matrix_) SUNMatDestroy(matrix_);
        if (absoluteTolerance_) N_VDestroy(absoluteTolerance_);
        if (state_) N_VDestroy(state_);
        if (context_) SUNContext_Free(&context_);
        memory_ = nullptr;
        linearSolver_ = nullptr;
        matrix_ = nullptr;
        absoluteTolerance_ = nullptr;
        state_ = nullptr;
        context_ = nullptr;
    }

    void reinitializeAt(double time)
    {
        if (CVodeReInit(memory_, time, state_)
            || CVodeSVtolerances(
                memory_, tolerances_.relative, absoluteTolerance_)
            || CVodeSetUserData(memory_, this)
            || CVodeSetMaxNumSteps(memory_, 2000000)
            || CVodeSetMaxStep(memory_, maxInternalStep_))
            throw std::runtime_error("cvode-reinitialize-failed");
        int rootDirections[rootCount] = {1, -1};
        if (CVodeSetRootDirection(memory_, rootDirections))
            throw std::runtime_error("cvode-root-direction-failed");
    }

    void rawAdvanceTo(double targetTime)
    {
        while (targetTime > currentTime_)
        {
            sunrealtype reached = currentTime_;
            const int status = CVode(
                memory_, targetTime, state_, &reached, CV_NORMAL);
            currentTime_ = reached;
            if (status == CV_ROOT_RETURN)
            {
                int roots[rootCount] = {0, 0};
                if (CVodeGetRootInfo(memory_, roots))
                    throw std::runtime_error("cvode-root-info-failed");
                if (roots[0] != 0 && !ignitionDetected_)
                {
                    ignitionDetected_ = true;
                    ignitionDelay_ = currentTime_;
                }
                if (roots[1] != 0 && !fuelDepleted_)
                {
                    sunrealtype* values = N_VGetArrayPointer(state_);
                    // The root-interpolated complete state already contains
                    // the corresponding product and enthalpy increments.  Set
                    // only the root variable to exact zero before restarting.
                    if (std::abs(values[0]) > 1.0e-10)
                        throw std::runtime_error(
                            "fuel-depletion-root-residual-too-large");
                    values[0] = 0.0;
                    fuelDepleted_ = true;
                    ++depletionEventCount_;
                    reinitializeAt(currentTime_);
                }
                continue;
            }
            if (status < 0)
                throw std::runtime_error(
                    "cvode-advance-failed:" + std::to_string(status));
        }
    }

    bool acceptedCompositionIsNonnegative()
    {
        sunrealtype* values = N_VGetArrayPointer(state_);
        bool clamped = false;
        for (int species = 0; species < speciesCount; ++species)
        {
            if (!std::isfinite(values[species]))
                throw std::runtime_error("nonfinite-accepted-composition");
            if (values[species] < -1.0e-12) return false;
            if (values[species] < 0.0)
            {
                values[species] = 0.0;
                clamped = true;
            }
        }
        if (clamped)
        {
            ++roundoffClampCount_;
            reinitializeAt(currentTime_);
        }
        return true;
    }

    void advanceInterval(double targetTime, int depth)
    {
        if (targetTime <= currentTime_) return;
        const double startTime = currentTime_;
        const sunrealtype* startValues = N_VGetArrayPointer(state_);
        const std::vector<double> savedState(
            startValues, startValues + equationCount);
        const bool savedIgnitionDetected = ignitionDetected_;
        const double savedIgnitionDelay = ignitionDelay_;
        const bool savedFuelDepleted = fuelDepleted_;
        const int savedDepletionEvents = depletionEventCount_;
        const int savedRoundoffClamps = roundoffClampCount_;

        rawAdvanceTo(targetTime);
        if (acceptedCompositionIsNonnegative()) return;

        if (depth >= 40)
            throw std::runtime_error("positivity-interval-refinement-exhausted");
        const double midpoint = startTime + 0.5 * (targetTime - startTime);
        if (!(midpoint > startTime && midpoint < targetTime))
            throw std::runtime_error("positivity-interval-not-divisible");

        sunrealtype* restored = N_VGetArrayPointer(state_);
        std::copy(savedState.begin(), savedState.end(), restored);
        currentTime_ = startTime;
        ignitionDetected_ = savedIgnitionDetected;
        ignitionDelay_ = savedIgnitionDelay;
        fuelDepleted_ = savedFuelDepleted;
        depletionEventCount_ = savedDepletionEvents;
        roundoffClampCount_ = savedRoundoffClamps;
        ++intervalBisectionCount_;
        reinitializeAt(currentTime_);
        advanceInterval(midpoint, depth + 1);
        advanceInterval(targetTime, depth + 1);
    }

public:
    ConstantPressureSensibleEnthalpyReactor()
    :
        thermo_(Foam::nasa7Thermo()),
        speciesNames_(chemistry_.speciesNames()),
        rhsMassFractions_(speciesCount, 0.0),
        rhsSources_(speciesCount, 0.0)
    {}

    ~ConstantPressureSensibleEnthalpyReactor() { destroy(); }

    ConstantPressureSensibleEnthalpyReactor(
        const ConstantPressureSensibleEnthalpyReactor&) = delete;
    ConstantPressureSensibleEnthalpyReactor& operator=(
        const ConstantPressureSensibleEnthalpyReactor&) = delete;

    void initialize(
        double pressure,
        double temperature,
        const std::vector<double>& composition,
        double maxInternalStep,
        const SolverTolerances& tolerances)
    {
        if (!std::isfinite(pressure) || pressure <= 0.0
            || !std::isfinite(temperature)
            || temperature < minimumTemperature
            || temperature > maximumTemperature
            || composition.size() != speciesCount
            || !std::isfinite(maxInternalStep)
            || maxInternalStep < 0.0)
            throw std::invalid_argument("invalid-initial-state");
        pressure_ = pressure;
        maxInternalStep_ = maxInternalStep;
        tolerances_ = tolerances;
        ignitionTemperature_ = temperature + 400.0;
        ignitionDetected_ = false;
        ignitionDelay_ = -1.0;
        currentTime_ = 0.0;
        fuelDepleted_ = composition[0] <= 0.0;
        depletionEventCount_ = 0;
        intervalBisectionCount_ = 0;
        roundoffClampCount_ = 0;

        if (SUNContext_Create(SUN_COMM_NULL, &context_))
            throw std::runtime_error("sun-context-create-failed");
        state_ = N_VNew_Serial(equationCount, context_);
        absoluteTolerance_ = N_VNew_Serial(equationCount, context_);
        if (!state_ || !absoluteTolerance_)
            throw std::runtime_error("sun-vector-create-failed");

        sunrealtype* values = N_VGetArrayPointer(state_);
        sunrealtype* absoluteTolerances =
            N_VGetArrayPointer(absoluteTolerance_);
        for (int species = 0; species < speciesCount; ++species)
        {
            values[species] = composition[species];
            absoluteTolerances[species] = tolerances_.speciesAbsolute;
        }
        values[speciesCount] = sensibleEnthalpy(composition, temperature);
        absoluteTolerances[speciesCount] = tolerances_.enthalpyAbsolute;

        memory_ = CVodeCreate(CV_BDF, context_);
        if (!memory_ || CVodeInit(memory_, &rhs, 0.0, state_)
            || CVodeSVtolerances(
                memory_, tolerances_.relative, absoluteTolerance_)
            || CVodeSetUserData(memory_, this)
            || CVodeRootInit(memory_, rootCount, &root))
            throw std::runtime_error("cvode-initialize-failed");
        int rootDirections[rootCount] = {1, -1};
        if (CVodeSetRootDirection(memory_, rootDirections))
            throw std::runtime_error("cvode-root-direction-failed");

        matrix_ = SUNDenseMatrix(equationCount, equationCount, context_);
        linearSolver_ = SUNLinSol_Dense(state_, matrix_, context_);
        if (!matrix_ || !linearSolver_
            || CVodeSetLinearSolver(memory_, linearSolver_, matrix_)
            || CVodeSetMaxNumSteps(memory_, 2000000)
            || CVodeSetMaxStep(memory_, maxInternalStep_))
            throw std::runtime_error("cvode-linear-solver-initialize-failed");
    }

    void advanceTo(double targetTime)
    {
        if (!memory_ || !std::isfinite(targetTime)
            || targetTime < currentTime_)
            throw std::invalid_argument("invalid-target-time");
        advanceInterval(targetTime, 0);
        if (currentTime_ != targetTime)
            throw std::runtime_error("accepted-physical-time-mismatch");
    }

    Snapshot snapshot(double requestedTime)
    {
        const sunrealtype* values = N_VGetArrayPointer(state_);
        std::vector<double> composition;
        evaluationComposition(values, composition);
        const double sensible = values[speciesCount];
        const double temperature = recoverTemperature(composition, sensible);
        const double recoveredDensity = density(composition, temperature);
        std::vector<double> sources;
        double heatRelease = 0.0;
        chemistry_.computeSourceTerms(
            temperature,
            pressure_,
            recoveredDensity,
            composition,
            sources,
            heatRelease);
        std::vector<double> massFractionRates(speciesCount, 0.0);
        double temperatureRate = heatRelease / recoveredDensity;
        for (int species = 0; species < speciesCount; ++species)
        {
            massFractionRates[species] =
                sources[species] / recoveredDensity;
            temperatureRate -= thermo_.hsens_mole(
                speciesNames_[species], temperature)
                / thermo_.MW(speciesNames_[species])
                * massFractionRates[species];
        }
        const double heatCapacity = thermo_.cp_mass_mix(
            speciesNames_, composition, temperature);
        if (!std::isfinite(temperatureRate)
            || !std::isfinite(heatCapacity) || heatCapacity <= 0.0)
            throw std::runtime_error("snapshot-temperature-rate-failed");
        temperatureRate /= heatCapacity;
        return Snapshot{
            requestedTime,
            temperature,
            pressure_,
            recoveredDensity,
            sensible,
            totalEnthalpy(composition, sensible),
            heatRelease,
            temperatureRate,
            std::move(composition),
            std::move(massFractionRates)};
    }

    bool ignitionDetected() const { return ignitionDetected_; }
    double ignitionDelay() const { return ignitionDelay_; }
    int depletionEventCount() const { return depletionEventCount_; }
    int intervalBisectionCount() const { return intervalBisectionCount_; }
    int roundoffClampCount() const { return roundoffClampCount_; }
    double acceptedPhysicalTime() const { return currentTime_; }
};

void emitSnapshot(
    std::size_t index,
    const Snapshot& sample,
    bool strictEquivalence)
{
    std::cout << "sample=" << index
        << "|time_s=" << sample.time
        << "|temperature_k=" << sample.temperature
        << "|pressure_pa=" << sample.pressure
        << "|density_kg_m3=" << sample.density
        << "|sensible_enthalpy_j_kg=" << sample.sensibleEnthalpy
        << "|total_enthalpy_j_kg=" << sample.totalEnthalpy
        << "|qdot_w_m3=" << sample.heatRelease;
    if (strictEquivalence)
    {
        std::cout << "|temperature_rate_k_s=" << sample.temperatureRate
            << "|mass_fraction_rates_s_inv=";
        for (std::size_t species = 0;
             species < sample.massFractionRates.size(); ++species)
        {
            if (species) std::cout << ',';
            std::cout << sample.massFractionRates[species];
        }
    }
    std::cout << "|mass_fractions=";
    for (std::size_t species = 0;
         species < sample.massFractions.size(); ++species)
    {
        if (species) std::cout << ',';
        std::cout << sample.massFractions[species];
    }
    std::cout << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc < 14)
            throw std::invalid_argument("arguments-required");
        const std::string caseId = argv[1];
        if (!validCaseId(caseId))
            throw std::invalid_argument("invalid-case-id");
        const double initialTemperature = parseFinite(argv[2], "temperature");
        const double pressure = parseFinite(argv[3], "pressure");
        const std::size_t sampleCount = parseCount(argv[4]);
        const double maxInternalStep = parseFinite(
            argv[5], "max-internal-step");
        if (maxInternalStep < 0.0)
            throw std::invalid_argument("negative-max-internal-step");
        const int expectedArguments = 12 + static_cast<int>(sampleCount);
        const bool strictEquivalence = argc == expectedArguments + 1
            && std::string(argv[expectedArguments])
                == "--solver-tolerance-profile=strict-equivalence-v1";
        if (argc != expectedArguments && !strictEquivalence)
            throw std::invalid_argument("argument-count-mismatch");
        const SolverTolerances tolerances = solverTolerances(
            strictEquivalence);

        std::vector<double> composition(6, 0.0);
        double compositionSum = 0.0;
        for (int species = 0; species < 6; ++species)
        {
            composition[species] = parseFinite(
                argv[6 + species], "mass-fraction");
            if (composition[species] < 0.0)
                throw std::invalid_argument("negative-mass-fraction");
            compositionSum += composition[species];
        }
        if (std::abs(compositionSum - 1.0) > 2.0e-10)
            throw std::invalid_argument("mass-fraction-sum-mismatch");

        std::vector<double> sampleTimes(sampleCount, 0.0);
        for (std::size_t sample = 0; sample < sampleCount; ++sample)
        {
            sampleTimes[sample] = parseFinite(
                argv[12 + static_cast<int>(sample)], "sample-time");
            if (sampleTimes[sample] < 0.0
                || (sample > 0
                    && sampleTimes[sample] <= sampleTimes[sample - 1]))
                throw std::invalid_argument("sample-times-not-increasing");
        }
        if (sampleTimes.front() != 0.0 || sampleTimes.back() > 1.0)
            throw std::invalid_argument("sample-time-envelope-mismatch");

        ConstantPressureSensibleEnthalpyReactor reactor;
        reactor.initialize(
            pressure,
            initialTemperature,
            composition,
            maxInternalStep,
            tolerances);
        std::vector<Snapshot> samples;
        samples.reserve(sampleCount);
        for (const double sampleTime : sampleTimes)
        {
            reactor.advanceTo(sampleTime);
            samples.push_back(reactor.snapshot(sampleTime));
        }

        std::cout << std::setprecision(17);
        std::cout << "schema=rfam.r2-0d-rfam-trajectory/v1\n";
        std::cout << "status=implementation-protocol-only-not-fit-eligible\n";
        std::cout << "case_id=" << caseId << '\n';
        std::cout << "reactor=ideal-gas-constant-pressure-sensible-enthalpy\n";
        std::cout << "chemistry=rfam-twostepco-candidate\n";
        std::cout << "energy_basis=sensible-enthalpy-plus-explicit-formation-heat\n";
        std::cout << "max_internal_step_s=" << maxInternalStep << '\n';
        if (strictEquivalence)
        {
            std::cout << "solver_tolerance_profile="
                << tolerances.profile << '\n';
            std::cout << "solver_relative_tolerance="
                << tolerances.relative << '\n';
            std::cout << "solver_species_absolute_tolerance="
                << tolerances.speciesAbsolute << '\n';
            std::cout << "solver_enthalpy_absolute_tolerance="
                << tolerances.enthalpyAbsolute << '\n';
        }
        std::cout << "species_order=C3H8,O2,CO,CO2,H2O,N2\n";
        std::cout << "sample_count=" << samples.size() << '\n';
        std::cout << "ignition_definition=first-root-T-minus-T0-equals-400K\n";
        std::cout << "ignition_detected="
            << (reactor.ignitionDetected() ? 1 : 0) << '\n';
        std::cout << "ignition_delay_s=" << reactor.ignitionDelay() << '\n';
        std::cout << "positivity_strategy="
            << "fuel-depletion-root-plus-rejected-interval-bisection\n";
        std::cout << "depletion_event_count="
            << reactor.depletionEventCount() << '\n';
        std::cout << "interval_bisection_count="
            << reactor.intervalBisectionCount() << '\n';
        std::cout << "roundoff_clamp_count="
            << reactor.roundoffClampCount() << '\n';
        std::cout << "accepted_physical_time_s="
            << reactor.acceptedPhysicalTime() << '\n';
        std::cout << "implicit_time_extension_s=0\n";
        std::cout << "fit_approved=0\n";
        std::cout << "physical_validation_approved=0\n";
        std::cout << "production_ready=0\n";
        for (std::size_t sample = 0; sample < samples.size(); ++sample)
            emitSnapshot(sample, samples[sample], strictEquivalence);
        std::cout << "R2_0D_RFAM_RUNNER_OK\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "R2_0D_RFAM_RUNNER_FAIL " << error.what() << '\n';
        return 2;
    }
}
