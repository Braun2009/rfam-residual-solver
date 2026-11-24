/*---------------------------------------------------------------------------*\
    RFAM TwoStepCO instantaneous-rate probe.

    The executable evaluates the exact production implementation at a caller-
    supplied T/rho/Y state.  It also evaluates the YAML-driven implementation
    and emits machine-readable progress, source, and formation-energy values.
    Artifact bytes are bound at CMake configure time and rechecked at runtime.
\*---------------------------------------------------------------------------*/

#include "src/combustion/FiniteRate/DynamicChemistry.H"
#include "src/combustion/FiniteRate/TwoStepCO.H"
#include "src/thermo/SpeciesRegistry.H"
#include "src/thermo/ThermoDataProvenance.H"

#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef RFAM_R2_GENERATOR_SHA256
#error "RFAM_R2_GENERATOR_SHA256 must be bound at build time"
#endif
#ifndef RFAM_R2_REFERENCE_SHA256
#error "RFAM_R2_REFERENCE_SHA256 must be bound at build time"
#endif
#ifndef RFAM_R2_CANTERA_MECHANISM_SHA256
#error "RFAM_R2_CANTERA_MECHANISM_SHA256 must be bound at build time"
#endif
#ifndef RFAM_R2_THERMO_SHA256
#error "RFAM_R2_THERMO_SHA256 must be bound at build time"
#endif
#ifndef RFAM_R2_RFAM_MECHANISM_SHA256
#error "RFAM_R2_RFAM_MECHANISM_SHA256 must be bound at build time"
#endif
#ifndef RFAM_R2_SPECIES_REGISTRY_SHA256
#error "RFAM_R2_SPECIES_REGISTRY_SHA256 must be bound at build time"
#endif
#ifndef RFAM_R2_MASS_ACTION_THERMO_SHA256
#error "RFAM_R2_MASS_ACTION_THERMO_SHA256 must be bound at build time"
#endif
#ifndef RFAM_R2_NASA7_THERMO_SHA256
#error "RFAM_R2_NASA7_THERMO_SHA256 must be bound at build time"
#endif
#ifndef RFAM_R2_MECHANISM_READER_SHA256
#error "RFAM_R2_MECHANISM_READER_SHA256 must be bound at build time"
#endif
#ifndef RFAM_R2_THERMO_PROVENANCE_HEADER_SHA256
#error "RFAM_R2_THERMO_PROVENANCE_HEADER_SHA256 must be bound at build time"
#endif
#ifndef RFAM_R2_THERMO_PROVENANCE_SOURCE_SHA256
#error "RFAM_R2_THERMO_PROVENANCE_SOURCE_SHA256 must be bound at build time"
#endif

namespace
{

using Foam::DynamicChemistry;
using Foam::SpeciesRegistry;
using Foam::TwoStepCO;

constexpr std::size_t nSpecies = 6;
const std::array<std::string, nSpecies> species = {
    "C3H8", "O2", "CO", "CO2", "H2O", "N2"};
const std::array<double, nSpecies> molecularWeights = {
    0.04409562, 0.0319988, 0.0280101,
    0.0440095, 0.01801528, 0.0280134};
const std::array<std::array<int, 4>, nSpecies> elements = {{
    {{3, 8, 0, 0}},
    {{0, 0, 2, 0}},
    {{1, 0, 1, 0}},
    {{1, 0, 2, 0}},
    {{0, 2, 1, 0}},
    {{0, 0, 0, 2}},
}};

std::string readBinaryFile(const std::string& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input.good())
        throw std::runtime_error("provenance-open-failed:" + path);
    return std::string(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>());
}

bool validSha256(const std::string& value)
{
    if (value.size() != 64) return false;
    for (const unsigned char character : value)
        if (!std::isdigit(character)
            && !(character >= 'a' && character <= 'f'))
            return false;
    return true;
}

void requireArtifact(
    const std::string& path,
    const std::string& expected,
    const std::string& label)
{
    if (!validSha256(expected))
        throw std::runtime_error(label + "-build-sha256-invalid");
    const std::string actual =
        Foam::thermoProvenance::sha256Hex(readBinaryFile(path));
    if (actual != expected)
        throw std::runtime_error(label + "-sha256-mismatch");
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
        throw std::runtime_error("invalid-number:" + label);
    }
    if (consumed != std::string(text).size() || !std::isfinite(value))
        throw std::runtime_error("invalid-number:" + label);
    return value;
}

void requireSpeciesContract(
    const TwoStepCO& twoStep,
    const DynamicChemistry& chemistry,
    const SpeciesRegistry& registry)
{
    if (twoStep.speciesNames().size() != nSpecies
        || chemistry.speciesNames().size() != nSpecies)
        throw std::runtime_error("rfam-mechanism-species-count-mismatch");
    for (std::size_t i = 0; i < nSpecies; ++i)
    {
        const int global = registry.index(species[i]);
        if (global < 0
            || twoStep.speciesNames()[i] != species[i]
            || chemistry.speciesNames()[i] != species[i])
            throw std::runtime_error("rfam-mechanism-species-order-mismatch");
        const double registryMw = registry.molecularWeight(global);
        if (std::abs(registryMw - molecularWeights[i]) > 1.0e-14
            || std::abs(twoStep.molecularWeights()[i] - registryMw) > 1.0e-14
            || std::abs(chemistry.molecularWeights()[i] - registryMw) > 1.0e-14)
            throw std::runtime_error("rfam-mechanism-mw-mismatch");
        for (std::size_t element = 0; element < 4; ++element)
        {
            const int registryCount[] = {
                registry.nC(global), registry.nH(global),
                registry.nO(global), registry.nN(global)};
            if (registryCount[element] != elements[i][element])
                throw std::runtime_error("rfam-species-element-mismatch");
        }
    }
}

void writeVector(const std::string& key, const std::vector<double>& values)
{
    std::cout << key << '=';
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (i != 0) std::cout << ',';
        std::cout << values[i];
    }
    std::cout << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    // Eleven artifacts, T, p, rho, and six mass fractions.
    if (argc != 21)
    {
        std::cerr
            << "R2_INSTANTANEOUS_PROBE_FAIL provenance-arguments-required\n";
        return 2;
    }

    try
    {
        requireArtifact(
            argv[1], RFAM_R2_GENERATOR_SHA256, "generator");
        requireArtifact(
            argv[2], RFAM_R2_REFERENCE_SHA256, "reference");
        requireArtifact(
            argv[3], RFAM_R2_CANTERA_MECHANISM_SHA256,
            "cantera-mechanism");
        requireArtifact(argv[4], RFAM_R2_THERMO_SHA256, "thermo");
        requireArtifact(
            argv[5], RFAM_R2_RFAM_MECHANISM_SHA256, "rfam-mechanism");
        requireArtifact(
            argv[6], RFAM_R2_SPECIES_REGISTRY_SHA256,
            "species-registry");
        requireArtifact(
            argv[7], RFAM_R2_MASS_ACTION_THERMO_SHA256,
            "mass-action-thermo");
        requireArtifact(
            argv[8], RFAM_R2_NASA7_THERMO_SHA256,
            "nasa7-thermo");
        requireArtifact(
            argv[9], RFAM_R2_MECHANISM_READER_SHA256,
            "mechanism-reader");
        requireArtifact(
            argv[10], RFAM_R2_THERMO_PROVENANCE_HEADER_SHA256,
            "thermo-provenance-header");
        requireArtifact(
            argv[11], RFAM_R2_THERMO_PROVENANCE_SOURCE_SHA256,
            "thermo-provenance-source");

        const double temperature = parseFinite(argv[12], "temperature");
        const double pressure = parseFinite(argv[13], "pressure");
        const double density = parseFinite(argv[14], "density");
        if (temperature <= 0.0 || pressure <= 0.0 || density <= 0.0)
            throw std::runtime_error("state-must-be-positive");

        std::vector<double> massFractions(nSpecies, 0.0);
        double massFractionSum = 0.0;
        for (std::size_t i = 0; i < nSpecies; ++i)
        {
            massFractions[i] = parseFinite(
                argv[15 + static_cast<int>(i)],
                "mass-fraction-" + std::to_string(i));
            if (massFractions[i] < 0.0)
                throw std::runtime_error("mass-fraction-negative");
            massFractionSum += massFractions[i];
        }
        if (std::abs(massFractionSum - 1.0) > 2.0e-12)
            throw std::runtime_error("mass-fraction-sum-mismatch");

        const SpeciesRegistry registry;
        TwoStepCO twoStep;
        const TwoStepCO::ReactionRates rates =
            twoStep.reactionRates(temperature, density, massFractions);
        std::vector<double> twoStepSource;
        double twoStepQdot = 0.0;
        twoStep.computeSourceTerms(
            temperature, pressure, density, massFractions,
            twoStepSource, twoStepQdot);

        DynamicChemistry dynamic(argv[5]);
        requireSpeciesContract(twoStep, dynamic, registry);
        std::vector<double> dynamicSource;
        double dynamicQdot = 0.0;
        dynamic.computeSourceTerms(
            temperature, pressure, density, massFractions,
            dynamicSource, dynamicQdot);

        std::array<long double, 4> elementResidual = {
            0.0L, 0.0L, 0.0L, 0.0L};
        long double massResidual = 0.0L;
        long double formationEnergyResidual =
            static_cast<long double>(twoStepQdot);
        for (std::size_t i = 0; i < nSpecies; ++i)
        {
            massResidual += static_cast<long double>(twoStepSource[i]);
            const long double molarRate =
                static_cast<long double>(twoStepSource[i])
                / static_cast<long double>(registry.molecularWeight(
                    registry.index(species[i])));
            for (std::size_t element = 0; element < 4; ++element)
                elementResidual[element] += molarRate * elements[i][element];
            formationEnergyResidual += molarRate
                * static_cast<long double>(registry.formationEnthalpy(
                    registry.index(species[i])));
        }

        std::cout << std::setprecision(17);
        std::cout << "schema=rfam.r2-instantaneous-probe/v1\n";
        std::cout << "r1_forward_mol_m3_s=" << rates.r1Forward << '\n';
        std::cout << "r2_forward_mol_m3_s=" << rates.r2Forward << '\n';
        std::cout << "r2_reverse_mol_m3_s=" << rates.r2Reverse << '\n';
        std::cout << "r2_net_mol_m3_s=" << rates.r2Net() << '\n';
        writeVector("twostep_mass_source_kg_m3_s", twoStepSource);
        std::cout << "twostep_qdot_w_m3=" << twoStepQdot << '\n';
        writeVector("dynamic_mass_source_kg_m3_s", dynamicSource);
        std::cout << "dynamic_qdot_w_m3=" << dynamicQdot << '\n';
        std::cout << "mass_residual_kg_m3_s="
                  << static_cast<double>(massResidual) << '\n';
        std::cout << "formation_energy_residual_w_m3="
                  << static_cast<double>(formationEnergyResidual) << '\n';
        std::cout << "element_residual_C_mol_m3_s="
                  << static_cast<double>(elementResidual[0]) << '\n';
        std::cout << "element_residual_H_mol_m3_s="
                  << static_cast<double>(elementResidual[1]) << '\n';
        std::cout << "element_residual_O_mol_m3_s="
                  << static_cast<double>(elementResidual[2]) << '\n';
        std::cout << "element_residual_N_mol_m3_s="
                  << static_cast<double>(elementResidual[3]) << '\n';
        std::cout << "R2_INSTANTANEOUS_PROBE_OK\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "R2_INSTANTANEOUS_PROBE_FAIL " << error.what() << '\n';
        return 2;
    }
}
