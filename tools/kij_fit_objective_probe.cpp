#include "PhaseEquilibrium/ConservedPhaseMixture.H"
#include "PhaseEquilibrium/PRPhaseCaloric.H"
#include "PhaseEquilibrium/PRTPFlash.H"
#include "MixingRules/BinaryInteractionDB.H"
#include "MixingRules/SpeciesCriticalDB.H"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Foam;
using namespace Foam::phaseEquilibrium;

namespace
{

constexpr const char* schema = "rfam.kij-fit-objective-probe/v1";
constexpr const char* watermark = "SYNTHETIC_NONPHYSICAL_TEST_ONLY";
constexpr const char* inputHeader =
    "watermark,case_id,split,pressure_pa,temperature_k,"
    "pressure_standard_uncertainty_pa,liquid_x_o2,"
    "liquid_x_o2_standard_uncertainty,vapor_y_o2,"
    "vapor_y_o2_standard_uncertainty";

struct Options
{
    double referenceTestKij = std::numeric_limits<double>::quiet_NaN();
    double trialTestKij = std::numeric_limits<double>::quiet_NaN();
    double densityRelativeStandardUncertainty =
        std::numeric_limits<double>::quiet_NaN();
    double energyStandardUncertainty =
        std::numeric_limits<double>::quiet_NaN();
};

struct InputRow
{
    std::string id;
    std::string split;
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double pressureUncertainty = std::numeric_limits<double>::quiet_NaN();
    double liquidXOxygen = std::numeric_limits<double>::quiet_NaN();
    double liquidXUncertainty = std::numeric_limits<double>::quiet_NaN();
    double vaporYOxygen = std::numeric_limits<double>::quiet_NaN();
    double vaporYUncertainty = std::numeric_limits<double>::quiet_NaN();
};

struct MixtureState
{
    double rho = std::numeric_limits<double>::quiet_NaN();
    double e = std::numeric_limits<double>::quiet_NaN();
    double materialResidual = std::numeric_limits<double>::infinity();
    double energyResidual = std::numeric_limits<double>::infinity();
};

struct ObjectiveRow
{
    double referenceX = std::numeric_limits<double>::quiet_NaN();
    double referenceY = std::numeric_limits<double>::quiet_NaN();
    double trialX = std::numeric_limits<double>::quiet_NaN();
    double trialY = std::numeric_limits<double>::quiet_NaN();
    double referenceTieLineResidual = std::numeric_limits<double>::infinity();
    double logFugacityOxygen = std::numeric_limits<double>::quiet_NaN();
    double logFugacityPropane = std::numeric_limits<double>::quiet_NaN();
    MixtureState referenceMixture;
    MixtureState trialMixture;
    double equilibriumChi2 = std::numeric_limits<double>::quiet_NaN();
    double liquidXChi2 = std::numeric_limits<double>::quiet_NaN();
    double vaporYChi2 = std::numeric_limits<double>::quiet_NaN();
    double rhoChi2 = std::numeric_limits<double>::quiet_NaN();
    double energyChi2 = std::numeric_limits<double>::quiet_NaN();
    double totalChi2 = std::numeric_limits<double>::quiet_NaN();
};

double parseDouble(const std::string& text, const std::string& label)
{
    std::size_t consumed = 0;
    double value = std::numeric_limits<double>::quiet_NaN();
    try
    {
        value = std::stod(text, &consumed);
    }
    catch (const std::exception&)
    {
        throw std::runtime_error("invalid-number:" + label);
    }
    if (consumed != text.size() || !std::isfinite(value))
        throw std::runtime_error("invalid-number:" + label);
    return value;
}

double parsePositive(const std::string& text, const std::string& label)
{
    const double value = parseDouble(text, label);
    if (value <= 0.0) throw std::runtime_error("nonpositive:" + label);
    return value;
}

std::vector<std::string> splitCsv(const std::string& line)
{
    std::vector<std::string> fields;
    std::istringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(field);
    if (!line.empty() && line.back() == ',') fields.emplace_back();
    return fields;
}

Options parseOptions(int argc, char** argv)
{
    if (argc != 9)
        throw std::runtime_error("expected-four-explicit-test-only-options");
    Options result;
    for (int index = 1; index < argc; index += 2)
    {
        const std::string name = argv[index];
        const std::string value = argv[index + 1];
        if (name == "--reference-test-kij")
            result.referenceTestKij = parseDouble(value, name);
        else if (name == "--trial-test-kij")
            result.trialTestKij = parseDouble(value, name);
        else if (name == "--density-relative-standard-uncertainty")
            result.densityRelativeStandardUncertainty =
                parsePositive(value, name);
        else if (name == "--energy-standard-uncertainty-j-kg")
            result.energyStandardUncertainty = parsePositive(value, name);
        else
            throw std::runtime_error("unknown-option:" + name);
    }
    if (!std::isfinite(result.referenceTestKij)
        || !std::isfinite(result.trialTestKij)
        || !std::isfinite(result.densityRelativeStandardUncertainty)
        || !std::isfinite(result.energyStandardUncertainty))
        throw std::runtime_error("missing-option");
    if (result.referenceTestKij == result.trialTestKij)
        throw std::runtime_error("reference-and-trial-must-differ");
    return result;
}

InputRow parseInput(const std::string& line, std::size_t row)
{
    const std::vector<std::string> fields = splitCsv(line);
    if (fields.size() != 10)
        throw std::runtime_error("input-column-count:row-" + std::to_string(row));
    if (fields[0] != watermark)
        throw std::runtime_error("input-watermark:row-" + std::to_string(row));
    if (fields[1].find("SYNTHETIC_NONPHYSICAL_") != 0)
        throw std::runtime_error("input-case-id:row-" + std::to_string(row));
    if (fields[2] != "training" && fields[2] != "holdout")
        throw std::runtime_error("input-split:row-" + std::to_string(row));

    InputRow result;
    result.id = fields[1];
    result.split = fields[2];
    result.pressure = parsePositive(fields[3], "pressure");
    result.temperature = parsePositive(fields[4], "temperature");
    result.pressureUncertainty = parsePositive(fields[5], "pressure-uncertainty");
    result.liquidXOxygen = parseDouble(fields[6], "liquid-x-o2");
    result.liquidXUncertainty = parsePositive(fields[7], "liquid-x-uncertainty");
    result.vaporYOxygen = parseDouble(fields[8], "vapor-y-o2");
    result.vaporYUncertainty = parsePositive(fields[9], "vapor-y-uncertainty");
    if (result.liquidXOxygen <= 0.0 || result.liquidXOxygen >= 1.0
        || result.vaporYOxygen <= 0.0 || result.vaporYOxygen >= 1.0)
        throw std::runtime_error("input-composition-bound:row-" + std::to_string(row));
    return result;
}

ConservedPhaseProperties conservedPhase(
    const PRPhaseCaloricProperties& phase)
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

MixtureState assembleMixture(
    const PRTPFlash& flash,
    const PRPhaseCaloric& caloric,
    const InputRow& input,
    TPFlashResult& split)
{
    const double overallOxygen =
        0.5*(input.liquidXOxygen + input.vaporYOxygen);
    split = flash.flashMoleFractions(
        input.pressure,
        input.temperature,
        {overallOxygen, 1.0 - overallOxygen});
    if (!split.valid() || split.region != EquilibriumRegion::TWO_PHASE
        || !split.phaseSplitReady())
        throw std::runtime_error("two-phase-flash-required:" + input.id);

    const PRPhaseCaloricResult liquid = caloric.evaluate(
        input.pressure,
        input.temperature,
        split.liquidMoleFractions,
        PhaseRoot::LIQUID);
    const PRPhaseCaloricResult vapor = caloric.evaluate(
        input.pressure,
        input.temperature,
        split.vaporMoleFractions,
        PhaseRoot::VAPOR);
    if (!liquid.success() || !vapor.success())
        throw std::runtime_error("phase-caloric-failure:" + input.id);
    const ConservedMixtureResult mixture = assembleConservedPhaseMixture(
        input.pressure,
        split.vaporMoleFraction,
        conservedPhase(*liquid.state),
        conservedPhase(*vapor.state),
        split.feedMassFractions);
    if (!mixture.success())
        throw std::runtime_error("mixture-assembly-failure:" + input.id);

    MixtureState result;
    result.rho = mixture.state->rho;
    result.e = mixture.state->e;
    result.materialResidual = mixture.materialResidual;
    result.energyResidual = std::max(
        mixture.phaseEnergyIdentityResidual,
        mixture.mixtureEnergyIdentityResidual);
    return result;
}

double square(double value)
{
    return value*value;
}

ObjectiveRow evaluate(const InputRow& input, const Options& options)
{
    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    BinaryInteractionDB referenceKij = BinaryInteractionDB::createDefault();
    BinaryInteractionDB trialKij = BinaryInteractionDB::createDefault();
    referenceKij.add("O2", "C3H8", options.referenceTestKij);
    trialKij.add("O2", "C3H8", options.trialTestKij);
    const std::vector<std::string> species = {"O2", "C3H8"};

    const PRTPFlash referenceFlash(critical, referenceKij, species);
    const PRTPFlash trialFlash(critical, trialKij, species);
    const PRPhaseCaloric referenceCaloric(critical, referenceKij, species);
    const PRPhaseCaloric trialCaloric(critical, trialKij, species);
    if (!referenceCaloric.configured() || !trialCaloric.configured())
        throw std::runtime_error("caloric-not-configured");

    ObjectiveRow result;
    TPFlashResult referenceSplit;
    TPFlashResult trialSplit;
    result.referenceMixture = assembleMixture(
        referenceFlash, referenceCaloric, input, referenceSplit);
    result.trialMixture = assembleMixture(
        trialFlash, trialCaloric, input, trialSplit);
    result.referenceX = referenceSplit.liquidMoleFractions[0];
    result.referenceY = referenceSplit.vaporMoleFractions[0];
    result.trialX = trialSplit.liquidMoleFractions[0];
    result.trialY = trialSplit.vaporMoleFractions[0];
    result.referenceTieLineResidual = std::max(
        std::abs(result.referenceX - input.liquidXOxygen),
        std::abs(result.referenceY - input.vaporYOxygen));

    const PRPhaseState observedLiquid = trialFlash.kernel().phaseState(
        input.pressure,
        input.temperature,
        {input.liquidXOxygen, 1.0 - input.liquidXOxygen},
        PhaseRoot::LIQUID);
    const PRPhaseState observedVapor = trialFlash.kernel().phaseState(
        input.pressure,
        input.temperature,
        {input.vaporYOxygen, 1.0 - input.vaporYOxygen},
        PhaseRoot::VAPOR);
    if (!observedLiquid.valid() || !observedVapor.valid())
        throw std::runtime_error("observed-tie-line-root-failure:" + input.id);
    result.logFugacityOxygen =
        std::log(input.liquidXOxygen)
        + observedLiquid.lnFugacityCoefficients[0]
        - std::log(input.vaporYOxygen)
        - observedVapor.lnFugacityCoefficients[0];
    result.logFugacityPropane =
        std::log(1.0 - input.liquidXOxygen)
        + observedLiquid.lnFugacityCoefficients[1]
        - std::log(1.0 - input.vaporYOxygen)
        - observedVapor.lnFugacityCoefficients[1];

    // For this scaffold only, relative pressure uncertainty is the declared
    // dimensionless scale for the P-x-y log-fugacity equilibrium residual.
    const double equilibriumScale =
        input.pressureUncertainty/input.pressure;
    const double rhoScale = options.densityRelativeStandardUncertainty
        *std::abs(result.referenceMixture.rho);
    result.equilibriumChi2 =
        square(result.logFugacityOxygen/equilibriumScale)
        + square(result.logFugacityPropane/equilibriumScale);
    result.liquidXChi2 = square(
        (result.trialX - input.liquidXOxygen)/input.liquidXUncertainty);
    result.vaporYChi2 = square(
        (result.trialY - input.vaporYOxygen)/input.vaporYUncertainty);
    result.rhoChi2 = square(
        (result.trialMixture.rho - result.referenceMixture.rho)/rhoScale);
    result.energyChi2 = square(
        (result.trialMixture.e - result.referenceMixture.e)
        /options.energyStandardUncertainty);
    result.totalChi2 = result.equilibriumChi2 + result.liquidXChi2
        + result.vaporYChi2 + result.rhoChi2 + result.energyChi2;
    return result;
}

void writeHeader(const Options& options)
{
    std::cout << "# " << schema << '\n';
    std::cout << "# " << watermark
              << " reference_test_kij=" << options.referenceTestKij
              << " trial_test_kij=" << options.trialTestKij
              << " not_fitted=1 runtime_publishable=0\n";
    std::cout
        << "case_id,split,pressure_pa,temperature_k,observed_liquid_x_o2,"
        << "observed_vapor_y_o2,reference_liquid_x_o2,reference_vapor_y_o2,"
        << "trial_liquid_x_o2,trial_vapor_y_o2,reference_tie_line_max_abs_residual,"
        << "trial_log_fugacity_residual_o2,trial_log_fugacity_residual_c3h8,"
        << "reference_rho_kg_m3,trial_rho_kg_m3,reference_e_j_kg,trial_e_j_kg,"
        << "phase_equilibrium_chi2,liquid_x_chi2,vapor_y_chi2,rho_chi2,e_chi2,"
        << "total_chi2,reference_material_residual,trial_material_residual,"
        << "reference_energy_residual,trial_energy_residual\n";
}

void writeRow(const InputRow& input, const ObjectiveRow& value)
{
    std::cout << input.id << ',' << input.split << ','
              << input.pressure << ',' << input.temperature << ','
              << input.liquidXOxygen << ',' << input.vaporYOxygen << ','
              << value.referenceX << ',' << value.referenceY << ','
              << value.trialX << ',' << value.trialY << ','
              << value.referenceTieLineResidual << ','
              << value.logFugacityOxygen << ','
              << value.logFugacityPropane << ','
              << value.referenceMixture.rho << ','
              << value.trialMixture.rho << ','
              << value.referenceMixture.e << ','
              << value.trialMixture.e << ','
              << value.equilibriumChi2 << ','
              << value.liquidXChi2 << ','
              << value.vaporYChi2 << ','
              << value.rhoChi2 << ','
              << value.energyChi2 << ','
              << value.totalChi2 << ','
              << value.referenceMixture.materialResidual << ','
              << value.trialMixture.materialResidual << ','
              << value.referenceMixture.energyResidual << ','
              << value.trialMixture.energyResidual << '\n';
}

} // End anonymous namespace

int main(int argc, char** argv)
{
    try
    {
        std::cout.imbue(std::locale::classic());
        std::cout << std::setprecision(17);
        const Options options = parseOptions(argc, argv);
        std::string line;
        if (!std::getline(std::cin, line) || line != inputHeader)
            throw std::runtime_error("input-header-mismatch");
        writeHeader(options);

        std::size_t row = 1;
        std::size_t trainingRecords = 0;
        std::size_t holdoutRecords = 0;
        double trainingObjective = 0.0;
        double holdoutObjective = 0.0;
        while (std::getline(std::cin, line))
        {
            ++row;
            if (line.empty()) throw std::runtime_error("empty-input-row");
            const InputRow input = parseInput(line, row);
            const ObjectiveRow output = evaluate(input, options);
            if (!std::isfinite(output.totalChi2)
                || output.referenceTieLineResidual > 5.0e-9
                || output.referenceMixture.materialResidual > 1.0e-8
                || output.trialMixture.materialResidual > 1.0e-8
                || output.referenceMixture.energyResidual > 1.0e-8
                || output.trialMixture.energyResidual > 1.0e-8)
                throw std::runtime_error("objective-invariant-failure:" + input.id);
            writeRow(input, output);
            if (input.split == "training")
            {
                ++trainingRecords;
                trainingObjective += output.totalChi2;
            }
            else
            {
                ++holdoutRecords;
                holdoutObjective += output.totalChi2;
            }
        }
        if (trainingRecords == 0 || holdoutRecords == 0)
            throw std::runtime_error("both-splits-required");
        std::cout << "# RFAM_KIJ_FIT_OBJECTIVE_SYNTHETIC_TEST_ONLY_PASS"
                  << " training_records=" << trainingRecords
                  << " holdout_records=" << holdoutRecords
                  << " training_objective=" << trainingObjective
                  << " holdout_objective=" << holdoutObjective
                  << " fit_executed=0 runtime_kij_emitted=0"
                  << " binary_interaction_db_modified=0 production_ready=0\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "RFAM_KIJ_FIT_OBJECTIVE_PROBE_FAIL reason="
                  << error.what() << '\n';
        return 2;
    }
}
