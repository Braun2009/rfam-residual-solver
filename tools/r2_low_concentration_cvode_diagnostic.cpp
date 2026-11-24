/*
    Offline CVODE stiffness and transaction diagnostic for trace-fuel R1.

    The production TwoStepCO closure is compared with an exact fractional-
    power variant. Both run through the unchanged RealFluidCVODEWrapper, so
    the production trace-state smoothing, rho-e-Y thermodynamics, and solver
    setup remain active. No rate cap or kinetic parameter is changed.
*/

#include "FiniteRate/ArrheniusRate.H"
#include "FiniteRate/ChemistryModel.H"
#include "FiniteRate/MassActionThermo.H"
#include "FiniteRate/RealFluidCVODEWrapper.H"
#include "FiniteRate/TwoStepCO.H"
#include "SpeciesMapper.H"
#include "SpeciesRegistry.H"
#include "RealFluidState.H"

#include <cvode/cvode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

template<class Tag, typename Tag::type Member>
struct PrivateAccess
{
    friend typename Tag::type get(Tag) { return Member; }
};

struct CvodeMemoryTag
{
    using type = void* Foam::RealFluidCVODEWrapper::*;
    friend type get(CvodeMemoryTag);
};

template struct PrivateAccess<
    CvodeMemoryTag,
    &Foam::RealFluidCVODEWrapper::cvodeMem_>;

namespace
{

constexpr int speciesCount = 6;
constexpr double transition = 1.0e-6;
constexpr double fuelMolecularWeight = 0.04409562;
constexpr double oxygenMolecularWeight = 0.0319988;

const std::array<const char*, speciesCount> species = {
    "C3H8", "O2", "CO", "CO2", "H2O", "N2"};
const std::array<double, speciesCount> speciesMolecularWeights = {
    0.04409562, 0.0319988, 0.0280101, 0.0440095, 0.01801528, 0.0280134};

struct DiagnosticCase
{
    const char* id;
    double temperature;
    double oxygenConcentration;
};

struct ToleranceProfile
{
    const char* id;
    double relative;
    double speciesAbsolute;
    double energyAbsolute;
};

struct SolverStats
{
    long steps = -1;
    long rhsEvaluations = -1;
    long errorTestFailures = -1;
    long nonlinearIterations = -1;
    long nonlinearFailures = -1;
    long stepSolveFailures = -1;
};

struct RunResult
{
    std::string caseId;
    std::string mode;
    std::string profile;
    std::string horizon;
    double initialTemperature = 0.0;
    double initialPressure = 0.0;
    double density = 0.0;
    double referenceDepletionTime = 0.0;
    double targetTime = 0.0;
    bool forcedFailure = false;
    bool advanceSucceeded = false;
    bool committed = false;
    bool rollbackVerified = false;
    bool acceptedNonnegative = false;
    bool candidateStateWritten = false;
    bool candidateFinite = false;
    bool candidateNonnegative = false;
    bool cvodeAdvancedToOrBeyondTarget = false;
    std::string failureStage;
    double cvodeCurrentTime = 0.0;
    SolverStats stats;
    long chemistryCalls = 0;
    double minimumChemistryMassFraction = 0.0;
    double minimumChemistryFuelMassFraction = 0.0;
    double candidateMassFractionSum = 0.0;
    double minimumCandidateMassFraction = 0.0;
    double candidateFuelConcentration = 0.0;
    double publishedFuelConcentration = 0.0;
    double finalTemperature = 0.0;
    double finalPressure = 0.0;
    double finalSpecificEnergy = 0.0;
    double minimumFinalMassFraction = 0.0;
    double massResidual = 0.0;
    std::array<double, 4> elementResidual = {0.0, 0.0, 0.0, 0.0};
    double formationEnergyResidual = 0.0;
    std::vector<double> finalGlobalY;
};

class ExactFractionalTwoStepCO : public Foam::ChemistryModel
{
    Foam::ArrheniusRate rate1_{1.0e10, 0.0, 125520.0};
    Foam::ArrheniusRate rate2_{2.0e12, 0.0, 167360.0};

public:
    ExactFractionalTwoStepCO()
    {
        nSpecies_ = speciesCount;
        nReactions_ = 2;
        speciesNames_.assign(species.begin(), species.end());
        molecularWeights_.assign(
            speciesMolecularWeights.begin(), speciesMolecularWeights.end());
    }

    void computeSourceTerms(
        double temperature,
        double,
        double density,
        const std::vector<double>& massFractions,
        std::vector<double>& sources,
        double& heatRelease) const override
    {
        if (!std::isfinite(temperature) || temperature <= 0.0
            || !std::isfinite(density) || density <= 0.0
            || massFractions.size() != speciesCount)
            throw std::runtime_error("exact-chemistry-invalid-state");
        std::array<double, speciesCount> concentration{};
        for (int index = 0; index < speciesCount; ++index)
        {
            if (!std::isfinite(massFractions[index])
                || massFractions[index] < 0.0)
                throw std::runtime_error("exact-chemistry-invalid-composition");
            concentration[index] = density * massFractions[index]
                / speciesMolecularWeights[index];
        }

        double r1 = 0.0;
        if (concentration[0] > 1.0e-20 && concentration[1] > 1.0e-20)
            r1 = rate1_.rate(temperature)
                * std::pow(concentration[0], 0.1)
                * std::pow(concentration[1], 1.65);

        const double kf = rate2_.rate(temperature);
        static const std::vector<std::string> r2Species = {"CO", "O2", "CO2"};
        static const std::vector<double> r2Stoich = {-1.0, -0.5, 1.0};
        const double kc = Foam::concentrationEquilibriumConstant(
            temperature, r2Species, r2Stoich);
        double r2Forward = 0.0;
        if (concentration[2] > 0.0 && concentration[1] > 0.0)
            r2Forward = kf * concentration[2] * std::sqrt(concentration[1]);
        const double r2Reverse = (kf / kc) * concentration[3];
        const double r2 = r2Forward - r2Reverse;

        sources.assign(speciesCount, 0.0);
        sources[0] -= r1 * speciesMolecularWeights[0];
        sources[1] -= 3.5 * r1 * speciesMolecularWeights[1];
        sources[2] += 3.0 * r1 * speciesMolecularWeights[2];
        sources[4] += 4.0 * r1 * speciesMolecularWeights[4];
        sources[1] -= 0.5 * r2 * speciesMolecularWeights[1];
        sources[2] -= r2 * speciesMolecularWeights[2];
        sources[3] += r2 * speciesMolecularWeights[3];
        heatRelease = r1 * 1.19502e6 + r2 * 2.8298e5;
        if (!std::isfinite(heatRelease))
            throw std::runtime_error("exact-chemistry-nonfinite-source");
    }
};

class RecordingChemistry : public Foam::ChemistryModel
{
    Foam::ChemistryModel& delegate_;
    mutable long calls_ = 0;
    mutable double minimumMassFraction_ = std::numeric_limits<double>::infinity();
    mutable double minimumFuelMassFraction_ =
        std::numeric_limits<double>::infinity();

public:
    explicit RecordingChemistry(Foam::ChemistryModel& delegate)
    : delegate_(delegate)
    {
        nSpecies_ = delegate.nSpecies();
        nReactions_ = delegate.nReactions();
        speciesNames_ = delegate.speciesNames();
        molecularWeights_ = delegate.molecularWeights();
    }

    void computeSourceTerms(
        double temperature,
        double pressure,
        double density,
        const std::vector<double>& massFractions,
        std::vector<double>& sources,
        double& heatRelease) const override
    {
        ++calls_;
        for (double value : massFractions)
            minimumMassFraction_ = std::min(minimumMassFraction_, value);
        if (!massFractions.empty())
            minimumFuelMassFraction_ = std::min(
                minimumFuelMassFraction_, massFractions.front());
        delegate_.computeSourceTerms(
            temperature, pressure, density, massFractions, sources, heatRelease);
    }

    long calls() const { return calls_; }
    double minimumMassFraction() const
    {
        return std::isfinite(minimumMassFraction_) ? minimumMassFraction_ : 0.0;
    }
    double minimumFuelMassFraction() const
    {
        return std::isfinite(minimumFuelMassFraction_)
            ? minimumFuelMassFraction_ : 0.0;
    }
};

SolverStats solverStats(void* memory)
{
    SolverStats result;
    if (!memory) return result;
    CVodeGetNumSteps(memory, &result.steps);
    CVodeGetNumRhsEvals(memory, &result.rhsEvaluations);
    CVodeGetNumErrTestFails(memory, &result.errorTestFailures);
    CVodeGetNumNonlinSolvIters(memory, &result.nonlinearIterations);
    CVodeGetNumNonlinSolvConvFails(memory, &result.nonlinearFailures);
    CVodeGetNumStepSolveFails(memory, &result.stepSolveFailures);
    return result;
}

bool sameState(
    const std::vector<double>& firstY,
    double firstEnergy,
    double firstTemperature,
    double firstPressure,
    const std::vector<double>& secondY,
    double secondEnergy,
    double secondTemperature,
    double secondPressure)
{
    return firstY == secondY
        && firstEnergy == secondEnergy
        && firstTemperature == secondTemperature
        && firstPressure == secondPressure;
}

void closureResiduals(
    const Foam::SpeciesRegistry& registry,
    double density,
    const std::vector<double>& initialY,
    double initialEnergy,
    const std::vector<double>& finalY,
    double finalEnergy,
    RunResult& result)
{
    long double mass = 0.0L;
    std::array<long double, 4> elements = {0.0L, 0.0L, 0.0L, 0.0L};
    long double formation = 0.0L;
    for (int index = 0; index < registry.size(); ++index)
    {
        const long double massChange = static_cast<long double>(density)
            * static_cast<long double>(finalY[index] - initialY[index]);
        mass += massChange;
        const long double moles = massChange
            / static_cast<long double>(registry.molecularWeight(index));
        elements[0] += moles * registry.nC(index);
        elements[1] += moles * registry.nH(index);
        elements[2] += moles * registry.nO(index);
        elements[3] += moles * registry.nN(index);
        formation += moles * registry.formationEnthalpy(index);
    }
    result.massResidual = static_cast<double>(mass);
    for (int element = 0; element < 4; ++element)
        result.elementResidual[element] = static_cast<double>(elements[element]);
    result.formationEnergyResidual = static_cast<double>(
        static_cast<long double>(density) * (finalEnergy - initialEnergy)
        + formation);
}

RunResult runCase(
    const DiagnosticCase& diagnosticCase,
    const ToleranceProfile& profile,
    const std::string& mode,
    const std::string& horizon,
    double horizonFactor,
    bool forcedFailure)
{
    RunResult result;
    result.caseId = diagnosticCase.id;
    result.mode = mode;
    result.profile = profile.id;
    result.horizon = horizon;
    result.initialTemperature = diagnosticCase.temperature;
    result.forcedFailure = forcedFailure;

    Foam::SpeciesRegistry registry;
    std::unique_ptr<Foam::ChemistryModel> implementation;
    if (mode == "current-composed")
        implementation = std::make_unique<Foam::TwoStepCO>();
    else if (mode == "exact-power-state-smoothing")
        implementation = std::make_unique<ExactFractionalTwoStepCO>();
    else
        throw std::runtime_error("unknown-mode");
    RecordingChemistry chemistry(*implementation);
    Foam::SpeciesMapper mapper(registry, chemistry.speciesNames());
    const int fuelGlobalIndex = mapper.globalIndex(0);

    std::vector<double> localY(speciesCount, 0.0);
    const double fuelMassDensity = transition * fuelMolecularWeight;
    const double oxygenMassDensity = diagnosticCase.oxygenConcentration
        * oxygenMolecularWeight;
    result.density = fuelMassDensity + oxygenMassDensity;
    localY[0] = fuelMassDensity / result.density;
    localY[1] = oxygenMassDensity / result.density;
    std::vector<double> initialY;
    mapper.mapChemToGlobal(localY, initialY);

    Foam::RealFluidState initialState;
    initialState.setSpecies(registry.names());
    if (initialState.updateRhoTY(
            result.density, diagnosticCase.temperature, initialY)
        != Foam::RealFluidState::Status::SUCCESS)
        throw std::runtime_error("initial-thermo-state-failed:" + result.caseId);
    result.initialPressure = initialState.p();
    const double initialEnergy = initialState.e();

    Foam::TwoStepCO referenceChemistry;
    const auto rates = referenceChemistry.reactionRates(
        diagnosticCase.temperature, result.density, localY);
    if (!std::isfinite(rates.r1Forward) || rates.r1Forward <= 0.0)
        throw std::runtime_error("initial-r1-rate-invalid:" + result.caseId);
    result.referenceDepletionTime = transition / rates.r1Forward;
    result.targetTime = horizonFactor * result.referenceDepletionTime;

    Foam::RealFluidCVODEWrapper wrapper(&chemistry, registry, mapper);
    wrapper.setTolerances(
        profile.relative, profile.speciesAbsolute, profile.energyAbsolute);
    wrapper.init(
        0.0, result.density, initialY, initialEnergy,
        diagnosticCase.temperature);
    void* memory = wrapper.*get(CvodeMemoryTag{});
    if (forcedFailure)
    {
        if (CVodeSetMaxNumSteps(memory, 1)
            || CVodeSetMaxStep(memory, result.targetTime * 1.0e-6))
            throw std::runtime_error("forced-failure-configuration-failed");
    }

    std::vector<double> publishedY = initialY;
    double publishedEnergy = initialEnergy;
    double publishedTemperature = diagnosticCase.temperature;
    double publishedPressure = initialState.p();
    std::vector<double> candidateY = publishedY;
    double candidateEnergy = publishedEnergy;
    double candidateTemperature = publishedTemperature;
    double candidatePressure = publishedPressure;

    result.advanceSucceeded = wrapper.advanceTo(
        result.targetTime,
        candidateY,
        candidateEnergy,
        candidateTemperature,
        candidatePressure);
    result.stats = solverStats(memory);
    result.chemistryCalls = chemistry.calls();
    result.minimumChemistryMassFraction = chemistry.minimumMassFraction();
    result.minimumChemistryFuelMassFraction =
        chemistry.minimumFuelMassFraction();
    sunrealtype cvodeCurrentTime = 0.0;
    if (CVodeGetCurrentTime(memory, &cvodeCurrentTime) == CV_SUCCESS)
        result.cvodeCurrentTime = static_cast<double>(cvodeCurrentTime);
    result.cvodeAdvancedToOrBeyondTarget =
        result.cvodeCurrentTime >= result.targetTime;
    result.candidateStateWritten = candidateY != initialY
        || candidateEnergy != initialEnergy;

    result.candidateFinite = candidateY.size() == initialY.size()
        && std::isfinite(candidateEnergy)
        && std::isfinite(candidateTemperature) && candidateTemperature > 0.0
        && std::isfinite(candidatePressure) && candidatePressure > 0.0;
    double sum = 0.0;
    double minimum = std::numeric_limits<double>::infinity();
    for (double value : candidateY)
    {
        result.candidateFinite = result.candidateFinite && std::isfinite(value);
        sum += value;
        minimum = std::min(minimum, value);
    }
    result.candidateMassFractionSum = sum;
    result.minimumCandidateMassFraction = minimum;
    result.candidateNonnegative = result.candidateFinite && minimum >= 0.0;
    result.candidateFuelConcentration = result.density
        * candidateY[fuelGlobalIndex] / fuelMolecularWeight;
    bool candidateValid = result.advanceSucceeded && result.candidateFinite
        && minimum >= 0.0
        && std::abs(sum - 1.0) <= 1.0e-8;
    result.acceptedNonnegative = candidateValid;

    if (result.advanceSucceeded)
        result.failureStage = "none";
    else if (result.candidateStateWritten)
        result.failureStage = "post-integrator-state-rejected";
    else
        result.failureStage = "integrator-did-not-reach-target";

    if (candidateValid && !forcedFailure)
    {
        publishedY = candidateY;
        publishedEnergy = candidateEnergy;
        publishedTemperature = candidateTemperature;
        publishedPressure = candidatePressure;
        result.committed = true;
    }
    result.rollbackVerified = result.committed || sameState(
        publishedY, publishedEnergy, publishedTemperature, publishedPressure,
        initialY, initialEnergy, diagnosticCase.temperature, initialState.p());

    result.finalGlobalY = publishedY;
    result.publishedFuelConcentration = result.density
        * publishedY[fuelGlobalIndex] / fuelMolecularWeight;
    result.finalSpecificEnergy = publishedEnergy;
    result.finalTemperature = publishedTemperature;
    result.finalPressure = publishedPressure;
    result.minimumFinalMassFraction = *std::min_element(
        publishedY.begin(), publishedY.end());
    if (result.committed)
        closureResiduals(
            registry, result.density, initialY, initialEnergy,
            publishedY, publishedEnergy, result);
    return result;
}

void writeNumber(double value)
{
    if (!std::isfinite(value))
        throw std::runtime_error("nonfinite-output");
    std::cout << value;
}

void writeResult(const RunResult& result)
{
    std::cout << "{\"case\":\"" << result.caseId
        << "\",\"mode\":\"" << result.mode
        << "\",\"profile\":\"" << result.profile
        << "\",\"horizon\":\"" << result.horizon << "\"";
    std::cout << ",\"initial_temperature_k\":";
    writeNumber(result.initialTemperature);
    std::cout << ",\"initial_pressure_pa\":";
    writeNumber(result.initialPressure);
    std::cout << ",\"density_kg_m3\":";
    writeNumber(result.density);
    std::cout << ",\"reference_depletion_time_s\":";
    writeNumber(result.referenceDepletionTime);
    std::cout << ",\"target_time_s\":";
    writeNumber(result.targetTime);
    std::cout << ",\"forced_failure\":"
        << (result.forcedFailure ? "true" : "false")
        << ",\"advance_succeeded\":"
        << (result.advanceSucceeded ? "true" : "false")
        << ",\"committed\":" << (result.committed ? "true" : "false")
        << ",\"rollback_verified\":"
        << (result.rollbackVerified ? "true" : "false")
        << ",\"accepted_nonnegative\":"
        << (result.acceptedNonnegative ? "true" : "false")
        << ",\"candidate_state_written\":"
        << (result.candidateStateWritten ? "true" : "false")
        << ",\"candidate_finite\":"
        << (result.candidateFinite ? "true" : "false")
        << ",\"candidate_nonnegative\":"
        << (result.candidateNonnegative ? "true" : "false")
        << ",\"cvode_advanced_to_or_beyond_target\":"
        << (result.cvodeAdvancedToOrBeyondTarget ? "true" : "false")
        << ",\"failure_stage\":\"" << result.failureStage << "\"";
    std::cout << ",\"cvode_current_time_s\":";
    writeNumber(result.cvodeCurrentTime);
    std::cout << ",\"stats\":{";
    std::cout << "\"steps\":" << result.stats.steps
        << ",\"rhs_evaluations\":" << result.stats.rhsEvaluations
        << ",\"error_test_failures\":" << result.stats.errorTestFailures
        << ",\"nonlinear_iterations\":" << result.stats.nonlinearIterations
        << ",\"nonlinear_failures\":" << result.stats.nonlinearFailures
        << ",\"step_solve_failures\":" << result.stats.stepSolveFailures
        << '}';
    std::cout << ",\"chemistry_calls\":" << result.chemistryCalls
        << ",\"minimum_chemistry_mass_fraction\":";
    writeNumber(result.minimumChemistryMassFraction);
    std::cout << ",\"minimum_chemistry_fuel_mass_fraction\":";
    writeNumber(result.minimumChemistryFuelMassFraction);
    std::cout << ",\"candidate_mass_fraction_sum\":";
    writeNumber(result.candidateMassFractionSum);
    std::cout << ",\"minimum_candidate_mass_fraction\":";
    writeNumber(result.minimumCandidateMassFraction);
    std::cout << ",\"candidate_fuel_concentration_mol_m3\":";
    writeNumber(result.candidateFuelConcentration);
    std::cout << ",\"published_fuel_concentration_mol_m3\":";
    writeNumber(result.publishedFuelConcentration);
    std::cout << ",\"final_temperature_k\":";
    writeNumber(result.finalTemperature);
    std::cout << ",\"final_pressure_pa\":";
    writeNumber(result.finalPressure);
    std::cout << ",\"final_specific_energy_j_kg\":";
    writeNumber(result.finalSpecificEnergy);
    std::cout << ",\"minimum_final_mass_fraction\":";
    writeNumber(result.minimumFinalMassFraction);
    std::cout << ",\"mass_residual_kg_m3\":";
    writeNumber(result.massResidual);
    std::cout << ",\"element_residual_mol_m3\":[";
    for (int index = 0; index < 4; ++index)
    {
        if (index != 0) std::cout << ',';
        writeNumber(result.elementResidual[index]);
    }
    std::cout << "],\"formation_energy_residual_j_m3\":";
    writeNumber(result.formationEnergyResidual);
    std::cout << ",\"final_global_mass_fractions\":[";
    for (std::size_t index = 0; index < result.finalGlobalY.size(); ++index)
    {
        if (index != 0) std::cout << ',';
        writeNumber(result.finalGlobalY[index]);
    }
    std::cout << "]}";
}

} // namespace

int main(int argc, char**)
{
    try
    {
        if (argc != 1)
            throw std::runtime_error("unexpected-command-line-arguments");
        const std::array<DiagnosticCase, 5> cases = {{
            {"low_900_o2_1", 900.0, 1.0},
            {"mid_1200_o2_8", 1200.0, 8.0},
            {"hot_1600_o2_8", 1600.0, 8.0},
            {"hot_2000_o2_70", 2000.0, 70.0},
            {"ceiling_3400_o2_422", 3400.0, 422.49812895050843},
        }};
        const std::array<ToleranceProfile, 3> profiles = {{
            {"loose", 1.0e-4, 1.0e-9, 1.0},
            {"default", 1.0e-6, 1.0e-11, 1.0e-2},
            {"tight", 1.0e-9, 1.0e-14, 1.0e-5},
        }};
        const std::array<std::string, 2> modes = {{
            "current-composed", "exact-power-state-smoothing"}};

        std::vector<RunResult> results;
        for (const auto& diagnosticCase : cases)
        {
            for (const auto& mode : modes)
            {
                for (const auto& profile : profiles)
                    results.push_back(runCase(
                        diagnosticCase, profile, mode,
                        "comparison-half-reference", 0.5, false));
                results.push_back(runCase(
                    diagnosticCase, profiles[1], mode,
                    "pre-extinction-reference", 1.0, false));
                results.push_back(runCase(
                    diagnosticCase, profiles[1], mode,
                    "near-extinction-1p1-reference", 1.1, false));
                results.push_back(runCase(
                    diagnosticCase, profiles[1], mode,
                    "post-extinction-1p2-reference", 1.2, false));
                results.push_back(runCase(
                    diagnosticCase, profiles[1], mode,
                    "stress-twice-reference", 2.0, false));
            }
        }
        for (const auto& mode : modes)
            results.push_back(runCase(
                cases.back(), profiles[1], mode,
                "forced-rollback", 0.5, true));

        std::cout << std::setprecision(17);
        std::cout
            << "{\"schema\":\"rfam.r2-low-concentration-cvode-diagnostic/v1\""
            << ",\"status\":\"offline-stiffness-diagnostic-production-blocked\""
            << ",\"transition_mol_m3\":";
        writeNumber(transition);
        std::cout
            << ",\"options\":["
            << "{\"id\":\"current-composed\",\"classification\":"
            << "\"current-production-numerical-closure\"},"
            << "{\"id\":\"exact-power-state-smoothing\",\"classification\":"
            << "\"reference-only-changes-subtransition-kinetic-evaluation\"}]"
            << ",\"runs\":[";
        for (std::size_t index = 0; index < results.size(); ++index)
        {
            if (index != 0) std::cout << ',';
            writeResult(results[index]);
        }
        std::cout
            << "],\"rate_cap_used\":false"
            << ",\"kinetic_parameters_changed\":false"
            << ",\"production_ready\":false}\n"
            << "R2_LOW_CONCENTRATION_CVODE_DIAGNOSTIC_OK\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr
            << "R2_LOW_CONCENTRATION_CVODE_DIAGNOSTIC_FAIL "
            << error.what() << '\n';
        return 2;
    }
}
