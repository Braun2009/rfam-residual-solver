/*
    Offline probe for the low-concentration smoothing used by RFAM's
    irreversible fractional-order mass-action rates.

    This characterizes numerical behavior only. It does not validate the
    underlying reduced chemistry or authorize production chemistry.
*/

#include "ArrheniusRate.H"
#include "MassActionThermo.H"
#include "RhsCompositionPreparation.H"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr double transition = 1.0e-6;
constexpr double fuelExponent = 0.1;
constexpr double oxygenExponent = 1.65;
constexpr double fuelMolecularWeight = 0.04409562;
constexpr double reactionHeatMagnitude = 1.19502e6;

double parsePositive(const char* text, const std::string& label)
{
    char* end = nullptr;
    const double value = std::strtod(text, &end);
    if (end == text || *end != '\0' || !std::isfinite(value) || value <= 0.0)
        throw std::runtime_error("invalid-positive-number:" + label);
    return value;
}

double regularized(double concentration, double exponent)
{
    const double value = Foam::regularizedConcentrationPower(
        concentration, exponent);
    if (!std::isfinite(value) || value < 0.0)
        throw std::runtime_error("regularized-value-invalid");
    return value;
}

double stateRegularized(double value, double localTransition)
{
    const double result =
        Foam::RhsCompositionPreparation::regularizedNonnegative(
            value, localTransition);
    if (!std::isfinite(result) || result < 0.0)
        throw std::runtime_error("state-regularized-value-invalid");
    return result;
}

bool prepareRhsComposition(
    double rho,
    const Foam::SpeciesRegistry& registry,
    const Foam::SpeciesMapper& mapper,
    const std::vector<double>& chemistryMolecularWeights,
    const std::vector<double>& globalTemplate,
    std::vector<double>& chemistryY,
    std::vector<double>& globalY,
    double* preNormalizationMassFractionSum)
{
    return Foam::RhsCompositionPreparation::prepare(
        rho, registry, mapper, chemistryMolecularWeights, globalTemplate,
        chemistryY, globalY, preNormalizationMassFractionSum);
}

void writeNumber(double value)
{
    if (!std::isfinite(value))
        throw std::runtime_error("nonfinite-output");
    std::cout << value;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        std::cerr
            << "R2_LOW_CONCENTRATION_PROBE_FAIL "
            << "temperature-and-oxygen-bound-required\n";
        return 2;
    }

    try
    {
        const double temperatureBound = parsePositive(argv[1], "temperature");
        const double oxygenBound = parsePositive(argv[2], "oxygen");
        const std::vector<double> exponents = {0.1, 0.5, 1.0, 1.65};
        const std::vector<double> ratios = {
            0.0, 1.0e-12, 1.0e-9, 1.0e-6, 1.0e-3, 1.0e-2,
            0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999999,
            1.0, 1.000001, 10.0};

        const double derivativeStep = transition * 1.0e-7;
        const double zeroStep = transition * 1.0e-8;

        Foam::ArrheniusRate rate1(1.0e10, 0.0, 125520.0);
        const double k1 = rate1.rate(temperatureBound);
        const double maximumFuelFactor = regularized(
            transition, fuelExponent);
        const double progressBound = k1 * maximumFuelFactor
            * std::pow(oxygenBound, oxygenExponent);
        const double fuelMassSourceBound =
            progressBound * fuelMolecularWeight;
        const double heatReleaseBound =
            progressBound * reactionHeatMagnitude;
        const double fuelInventoryAtTransition =
            transition * fuelMolecularWeight;
        const double depletionTimeBound =
            fuelInventoryAtTransition / fuelMassSourceBound;

        std::cout << std::setprecision(17);
        std::cout
            << "{\"schema\":\"rfam.r2-low-concentration-regularization/v1\""
            << ",\"status\":\"numerical-characterization-only-production-blocked\""
            << ",\"transition_mol_m3\":";
        writeNumber(transition);
        std::cout << ",\"rows\":[";

        bool firstRow = true;
        for (double exponent : exponents)
        {
            for (double ratio : ratios)
            {
                const double concentration = transition * ratio;
                const double value = regularized(concentration, exponent);
                const double exact = concentration > 0.0
                    ? std::pow(concentration, exponent) : 0.0;
                if (!std::isfinite(exact) || exact < 0.0)
                    throw std::runtime_error("exact-value-invalid");
                if (!firstRow) std::cout << ',';
                firstRow = false;
                std::cout << "{\"exponent\":";
                writeNumber(exponent);
                std::cout << ",\"concentration_ratio\":";
                writeNumber(ratio);
                std::cout << ",\"concentration_mol_m3\":";
                writeNumber(concentration);
                std::cout << ",\"regularized\":";
                writeNumber(value);
                std::cout << ",\"exact\":";
                writeNumber(exact);
                std::cout << ",\"regularized_to_exact\":";
                if (concentration == 0.0)
                    std::cout << "null";
                else
                    writeNumber(value / exact);
                std::cout << '}';
            }
        }
        std::cout << "],\"state_rows\":[";

        for (std::size_t index = 0; index < ratios.size(); ++index)
        {
            const double ratio = ratios[index];
            const double concentration = transition * ratio;
            const double value = stateRegularized(concentration, transition);
            if (index != 0) std::cout << ',';
            std::cout << "{\"concentration_ratio\":";
            writeNumber(ratio);
            std::cout << ",\"input_concentration_mol_m3\":";
            writeNumber(concentration);
            std::cout << ",\"regularized_concentration_mol_m3\":";
            writeNumber(value);
            std::cout << ",\"regularized_to_input\":";
            if (concentration == 0.0)
                std::cout << "null";
            else
                writeNumber(value / concentration);
            std::cout << '}';
        }

        std::cout << "],\"pre_normalization_combined_rows\":[";
        firstRow = true;
        for (std::size_t exponentIndex = 0; exponentIndex < 2; ++exponentIndex)
        {
            const double exponent = exponents[exponentIndex];
            for (double ratio : ratios)
            {
                const double concentration = transition * ratio;
                const double stateValue = stateRegularized(
                    concentration, transition);
                const double combined = regularized(stateValue, exponent);
                const double exact = concentration > 0.0
                    ? std::pow(concentration, exponent) : 0.0;
                if (!firstRow) std::cout << ',';
                firstRow = false;
                std::cout << "{\"exponent\":";
                writeNumber(exponent);
                std::cout << ",\"concentration_ratio\":";
                writeNumber(ratio);
                std::cout << ",\"state_regularized_mol_m3\":";
                writeNumber(stateValue);
                std::cout << ",\"combined_factor\":";
                writeNumber(combined);
                std::cout << ",\"exact_factor\":";
                writeNumber(exact);
                std::cout << ",\"combined_to_exact\":";
                if (concentration == 0.0)
                    std::cout << "null";
                else
                    writeNumber(combined / exact);
                std::cout << '}';
            }
        }

        Foam::SpeciesRegistry registry;
        const std::vector<std::string> chemistryNames = {
            "C3H8", "O2", "CO", "CO2", "H2O", "N2"};
        Foam::SpeciesMapper mapper(registry, chemistryNames);
        std::vector<double> chemistryMolecularWeights;
        chemistryMolecularWeights.reserve(chemistryNames.size());
        for (int c = 0; c < mapper.nChem(); ++c)
            chemistryMolecularWeights.push_back(
                registry.molecularWeight(mapper.globalIndex(c)));

        struct RhsCase
        {
            const char* id;
            double rho;
            double fuelRatio;
            bool traceTemplate;
        };
        const std::vector<RhsCase> rhsCases = {
            {"aligned-rho1-ratio0.001", 1.0, 0.001, false},
            {"aligned-rho1-ratio0.5", 1.0, 0.5, false},
            {"trace-rho1-ratio0.5", 1.0, 0.5, true},
            {"trace-rho0.1-ratio0.5", 0.1, 0.5, true},
            {"trace-rho10-ratio0.5", 10.0, 0.5, true},
            {"trace-rho1-ratio0.9", 1.0, 0.9, true},
        };
        std::cout << "],\"rhs_composition_rows\":[";
        bool firstRhsRow = true;
        for (const RhsCase& rhsCase : rhsCases)
        {
            const double fuelTransition = transition
                *chemistryMolecularWeights[0]/rhsCase.rho;
            const double fuelMassFraction = fuelTransition*rhsCase.fuelRatio;
            const double remaining = 1.0 - fuelMassFraction;
            std::vector<double> chemistryY = {
                fuelMassFraction, 0.23*remaining, 0.0, 0.0, 0.0,
                0.77*remaining};
            std::vector<double> globalTemplate(9, 0.0);
            if (rhsCase.traceTemplate)
            {
                globalTemplate = {
                    0.20, 0.20, 0.01, 0.02, 0.03, 0.20, 0.10, 0.10, 0.14};
            }
            else
            {
                globalTemplate = {
                    0.20, 0.20, 0.00, 0.00, 0.00, 0.20, 0.10, 0.10, 0.20};
            }
            const double inputChemistrySum =
                std::accumulate(chemistryY.begin(), chemistryY.end(), 0.0);
            const double inputTemplateSum =
                std::accumulate(globalTemplate.begin(), globalTemplate.end(), 0.0);
            std::vector<double> normalizedGlobal(9, 0.0);
            double preNormalizationSum = 0.0;
            if (!prepareRhsComposition(
                    rhsCase.rho, registry, mapper, chemistryMolecularWeights,
                    globalTemplate, chemistryY, normalizedGlobal,
                    &preNormalizationSum))
                throw std::runtime_error("rhs-composition-preparation-failed");
            const double normalizedFuelConcentration =
                rhsCase.rho*chemistryY[0]/chemistryMolecularWeights[0];
            const double preNormalizedFuelConcentration =
                normalizedFuelConcentration*preNormalizationSum;
            const double shortcutFuelFactor = regularized(
                preNormalizedFuelConcentration, fuelExponent);
            const double rhsFuelFactor = regularized(
                normalizedFuelConcentration, fuelExponent);
            if (!firstRhsRow) std::cout << ',';
            firstRhsRow = false;
            std::cout << "{\"id\":\"" << rhsCase.id << "\""
                << ",\"rho_kg_m3\":";
            writeNumber(rhsCase.rho);
            std::cout << ",\"fuel_concentration_ratio\":";
            writeNumber(rhsCase.fuelRatio);
            std::cout << ",\"global_template_sum\":";
            writeNumber(inputTemplateSum);
            std::cout << ",\"input_chemistry_mass_fraction_sum\":";
            writeNumber(inputChemistrySum);
            std::cout << ",\"pre_normalization_mass_fraction_sum\":";
            writeNumber(preNormalizationSum);
            std::cout << ",\"normalized_global_mass_fractions\":[";
            for (std::size_t i = 0; i < normalizedGlobal.size(); ++i)
            {
                if (i != 0) std::cout << ',';
                writeNumber(normalizedGlobal[i]);
            }
            std::cout << "],\"normalized_chemistry_mass_fractions\":[";
            for (std::size_t i = 0; i < chemistryY.size(); ++i)
            {
                if (i != 0) std::cout << ',';
                writeNumber(chemistryY[i]);
            }
            std::cout << "],\"pre_normalized_fuel_concentration_mol_m3\":";
            writeNumber(preNormalizedFuelConcentration);
            std::cout << ",\"normalized_fuel_concentration_mol_m3\":";
            writeNumber(normalizedFuelConcentration);
            std::cout << ",\"shortcut_r1_fuel_factor\":";
            writeNumber(shortcutFuelFactor);
            std::cout << ",\"rhs_r1_fuel_factor\":";
            writeNumber(rhsFuelFactor);
            std::cout << '}';
        }

        std::cout << "],\"derivatives\":[";

        for (std::size_t index = 0; index < 2; ++index)
        {
            const double exponent = exponents[index];
            const double atTransition = regularized(transition, exponent);
            const double left = (atTransition - regularized(
                transition - derivativeStep, exponent)) / derivativeStep;
            const double right = (regularized(
                transition + derivativeStep, exponent) - atTransition)
                / derivativeStep;
            const double expected = exponent
                * std::pow(transition, exponent - 1.0);
            const double zeroRight = regularized(zeroStep, exponent) / zeroStep;
            const double zeroRightNormalized = zeroRight
                / std::pow(transition, exponent - 1.0);
            if (index != 0) std::cout << ',';
            std::cout << "{\"exponent\":";
            writeNumber(exponent);
            std::cout << ",\"transition_left\":";
            writeNumber(left);
            std::cout << ",\"transition_right\":";
            writeNumber(right);
            std::cout << ",\"transition_expected\":";
            writeNumber(expected);
            std::cout << ",\"zero_right_normalized\":";
            writeNumber(zeroRightNormalized);
            std::cout << '}';
        }

        const double stateAtTransition = stateRegularized(
            transition, transition);
        const double stateLeft = (stateAtTransition - stateRegularized(
            transition - derivativeStep, transition)) / derivativeStep;
        const double stateRight = (stateRegularized(
            transition + derivativeStep, transition) - stateAtTransition)
            / derivativeStep;
        const double stateZeroRight = stateRegularized(
            zeroStep, transition) / zeroStep;

        std::cout << "],\"state_derivative\":{";
        std::cout << "\"transition_left\":";
        writeNumber(stateLeft);
        std::cout << ",\"transition_right\":";
        writeNumber(stateRight);
        std::cout << ",\"transition_expected\":1.0";
        std::cout << ",\"zero_right_normalized\":";
        writeNumber(stateZeroRight);
        std::cout << "},\"source_bound\":{";
        std::cout << "\"temperature_k\":";
        writeNumber(temperatureBound);
        std::cout << ",\"oxygen_mol_m3\":";
        writeNumber(oxygenBound);
        std::cout << ",\"k1\":";
        writeNumber(k1);
        std::cout << ",\"maximum_fuel_factor\":";
        writeNumber(maximumFuelFactor);
        std::cout << ",\"progress_mol_m3_s\":";
        writeNumber(progressBound);
        std::cout << ",\"fuel_mass_source_kg_m3_s\":";
        writeNumber(fuelMassSourceBound);
        std::cout << ",\"heat_release_w_m3\":";
        writeNumber(heatReleaseBound);
        std::cout << ",\"fuel_inventory_kg_m3\":";
        writeNumber(fuelInventoryAtTransition);
        std::cout << ",\"depletion_time_s\":";
        writeNumber(depletionTimeBound);
        std::cout
            << "},\"r2_fractional_power_policy\":"
            << "\"reversible-exact-power-no-regularization\""
            << ",\"production_ready\":false}\n";
        std::cout << "R2_LOW_CONCENTRATION_PROBE_OK\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr
            << "R2_LOW_CONCENTRATION_PROBE_FAIL " << error.what() << '\n';
        return 2;
    }
}
