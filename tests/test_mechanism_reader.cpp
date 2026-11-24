#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "src/combustion/FiniteRate/DynamicChemistry.H"
#include "src/combustion/FiniteRate/MechanismReader.H"
#include "src/combustion/FiniteRate/TwoStepCO.H"
#include "src/thermo/NASA7Thermo.H"

using namespace Foam;

namespace
{

int passCount = 0;
int failCount = 0;

void check(const std::string& name, bool passed)
{
    std::cout << "  [" << (passed ? "PASS" : "FAIL") << "] " << name << "\n";
    passed ? ++passCount : ++failCount;
}

bool closeRelative(double a, double b, double tolerance = 2e-12)
{
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= tolerance * scale;
}

std::string formatConcentration(double concentration)
{
    std::ostringstream output;
    output << std::scientific << std::setprecision(6) << concentration;
    return output.str();
}

double independentKcR2(double T)
{
    static constexpr double Runiv = 8.314462618;
    static constexpr double pStd = 101325.0;
    const NASA7Thermo& thermo = nasa7Thermo();
    const double deltaH = thermo.h_mole("CO2", T)
        - thermo.h_mole("CO", T) - 0.5 * thermo.h_mole("O2", T);
    const double deltaS = thermo.s_mole("CO2", T)
        - thermo.s_mole("CO", T) - 0.5 * thermo.s_mole("O2", T);
    const double Kp = std::exp(-(deltaH - T * deltaS) / (Runiv * T));
    return Kp * std::sqrt(Runiv * T / pStd);
}

struct ChemistryState
{
    double rho = 0.0;
    std::vector<double> Y;
};

ChemistryState makeCoState(
    const ChemistryModel& chemistry,
    double concentrationO2,
    double concentrationCO,
    double concentrationCO2)
{
    ChemistryState state;
    state.Y.assign(chemistry.nSpecies(), 0.0);
    const int iO2 = chemistry.speciesIndex("O2");
    const int iCO = chemistry.speciesIndex("CO");
    const int iCO2 = chemistry.speciesIndex("CO2");
    const auto& MW = chemistry.molecularWeights();
    state.rho = concentrationO2 * MW[iO2]
        + concentrationCO * MW[iCO]
        + concentrationCO2 * MW[iCO2];
    state.Y[iO2] = concentrationO2 * MW[iO2] / state.rho;
    state.Y[iCO] = concentrationCO * MW[iCO] / state.rho;
    state.Y[iCO2] = concentrationCO2 * MW[iCO2] / state.rho;
    return state;
}

bool sourcesMatch(
    const ChemistryModel& reference,
    const ChemistryModel& candidate,
    double T,
    double p,
    const ChemistryState& state)
{
    std::vector<double> omegaReference;
    std::vector<double> omegaCandidate;
    double qdotReference = 0.0;
    double qdotCandidate = 0.0;
    reference.computeSourceTerms(
        T, p, state.rho, state.Y, omegaReference, qdotReference);
    candidate.computeSourceTerms(
        T, p, state.rho, state.Y, omegaCandidate, qdotCandidate);
    if (omegaReference.size() != omegaCandidate.size()) return false;
    for (std::size_t i = 0; i < omegaReference.size(); ++i)
        if (!closeRelative(omegaReference[i], omegaCandidate[i])) return false;
    return closeRelative(qdotReference, qdotCandidate);
}

bool constructionRejected(
    const MechanismData& mechanism,
    const std::string& expectedMessage)
{
    try
    {
        DynamicChemistry invalid(mechanism);
        (void)invalid;
    }
    catch (const std::runtime_error& error)
    {
        return std::string(error.what()).find(expectedMessage)
            != std::string::npos;
    }
    return false;
}

std::string readTextFile(const std::string& path)
{
    std::ifstream input(path);
    if (!input.good())
        throw std::runtime_error("cannot open " + path);
    return std::string(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>());
}

std::string replaceOnce(
    std::string text,
    const std::string& from,
    const std::string& to)
{
    const std::size_t position = text.find(from);
    if (position == std::string::npos)
        throw std::runtime_error("test fixture text not found: " + from);
    text.replace(position, from.size(), to);
    return text;
}

bool readerRejected(
    const std::string& contents,
    const std::string& expectedMessage)
{
    static int sequence = 0;
    const std::string path = "/tmp/rfam_invalid_mechanism_"
        + std::to_string(sequence++) + ".yaml";
    {
        std::ofstream output(path);
        if (!output.good()) return false;
        output << contents;
    }

    bool rejected = false;
    try
    {
        const MechanismData invalid = MechanismReader::read(path);
        (void)invalid;
    }
    catch (const std::runtime_error& error)
    {
        rejected = std::string(error.what()).find(expectedMessage)
            != std::string::npos;
    }
    std::remove(path.c_str());
    return rejected;
}

} // namespace

int main()
{
    std::cout << "=== Reversible TwoStepCO mechanism contract ===\n";

    MechanismData mechanism;
    try
    {
        mechanism = MechanismReader::read("tests/twostep_co.yaml");
    }
    catch (const std::exception& error)
    {
        std::cerr << "failed to read mechanism: " << error.what() << "\n";
        return 1;
    }

    check("mechanism name parsed", mechanism.name == "TwoStepCO");
    check("six species parsed", mechanism.species.size() == 6);
    check("two reactions parsed", mechanism.reactions.size() == 2);
    check("R1 remains irreversible", !mechanism.reactions[0].reversible);
    check("R2 is explicitly reversible", mechanism.reactions[1].reversible);
    check("R2 O2 exponent matches stoichiometry",
        closeRelative(mechanism.reactions[1].concExp.at("O2"), 0.5));

    TwoStepCO staticChemistry;
    DynamicChemistry dynamicChemistry(mechanism);

    ChemistryState reactants;
    reactants.rho = 15.0;
    reactants.Y = {0.216, 0.784, 0.0, 0.0, 0.0, 0.0};
    check("static and YAML mechanisms match for the reacting feed",
        sourcesMatch(
            staticChemistry, dynamicChemistry, 1500.0, 4.0e6, reactants));

    const double equilibriumTemperature = 2500.0;
    const double Kc = independentKcR2(equilibriumTemperature);
    const double concentrationCO = 1.0;
    const std::vector<double> oxygenConcentrations = {
        1e-12, 0.0625, 16.0, 256.0};
    for (double concentrationO2 : oxygenConcentrations)
    {
        const double concentrationCO2 =
            Kc * concentrationCO * std::sqrt(concentrationO2);
        const ChemistryState equilibrium = makeCoState(
            staticChemistry,
            concentrationO2,
            concentrationCO,
            concentrationCO2);
        std::vector<double> omega;
        double qdot = 0.0;
        staticChemistry.computeSourceTerms(
            equilibriumTemperature, 4.0e6, equilibrium.rho, equilibrium.Y,
            omega, qdot);

        const double kForward = 2.0e12
            * std::exp(-167360.0 / (8.314 * equilibriumTemperature));
        const double forwardRateScale = kForward * concentrationCO
            * std::sqrt(concentrationO2);
        const double netR2 = omega[staticChemistry.speciesIndex("CO2")]
            / staticChemistry.molecularWeights()[staticChemistry.speciesIndex("CO2")];
        const bool detailedBalance = std::abs(netR2)
            <= 2e-12 * std::max(1.0, forwardRateScale);
        const bool zeroHeat = std::abs(qdot)
            <= 2e-12 * std::max(1.0, forwardRateScale * 2.8298e5);
        check("R2 detailed balance at O2 concentration "
            + formatConcentration(concentrationO2), detailedBalance && zeroHeat);
        check("static and YAML mechanisms match at equilibrium O2 concentration "
            + formatConcentration(concentrationO2),
            sourcesMatch(
                staticChemistry, dynamicChemistry, equilibriumTemperature,
                4.0e6, equilibrium));
    }

    const double referenceO2 = 16.0;
    const double equilibriumCO2 = Kc * concentrationCO * std::sqrt(referenceO2);
    for (double quotientFactor : {0.5, 2.0})
    {
        const ChemistryState state = makeCoState(
            staticChemistry,
            referenceO2,
            concentrationCO,
            quotientFactor * equilibriumCO2);
        std::vector<double> omega;
        double qdot = 0.0;
        staticChemistry.computeSourceTerms(
            equilibriumTemperature, 4.0e6, state.rho, state.Y, omega, qdot);
        const double co2Source = omega[staticChemistry.speciesIndex("CO2")];
        const bool expectedDirection = quotientFactor < 1.0
            ? (co2Source > 0.0 && qdot > 0.0)
            : (co2Source < 0.0 && qdot < 0.0);
        check(quotientFactor < 1.0
            ? "R2 advances toward products below equilibrium"
            : "R2 dissociates products above equilibrium",
            expectedDirection);
        check("static and YAML mechanisms match away from equilibrium",
            sourcesMatch(
                staticChemistry, dynamicChemistry, equilibriumTemperature,
                4.0e6, state));
    }

    MechanismData invalidMechanism = mechanism;
    invalidMechanism.reactions[1].concExp["O2"] = 0.25;
    check("reversible empirical-order mismatch is rejected",
        constructionRejected(invalidMechanism, "exponent mismatch"));

    invalidMechanism = mechanism;
    invalidMechanism.species[invalidMechanism.speciesIndex("N2")].MW *= 2.0;
    check("nonreacting-species molecular weight drift is rejected",
        constructionRejected(invalidMechanism, "molecular weight mismatch"));

    invalidMechanism = mechanism;
    invalidMechanism.species[invalidMechanism.speciesIndex("O2")].MW *= 1.01;
    check("reacting-species molecular weight drift is rejected",
        constructionRejected(invalidMechanism, "molecular weight mismatch"));

    invalidMechanism = mechanism;
    invalidMechanism.species[invalidMechanism.speciesIndex("O2")].MW =
        std::numeric_limits<double>::quiet_NaN();
    check("nonfinite molecular weight is rejected",
        constructionRejected(invalidMechanism, "invalid molecular weight"));

    invalidMechanism = mechanism;
    invalidMechanism.species[invalidMechanism.speciesIndex("O2")].MW = 0.0;
    check("zero molecular weight is rejected",
        constructionRejected(invalidMechanism, "invalid molecular weight"));

    invalidMechanism = mechanism;
    invalidMechanism.species[invalidMechanism.speciesIndex("O2")].MW = -1.0;
    check("negative molecular weight is rejected",
        constructionRejected(invalidMechanism, "invalid molecular weight"));

    invalidMechanism = mechanism;
    invalidMechanism.species[invalidMechanism.speciesIndex("C3H8")].hf += 1.0;
    check("formation enthalpy drift is rejected",
        constructionRejected(invalidMechanism, "formation enthalpy mismatch"));

    invalidMechanism = mechanism;
    invalidMechanism.species[invalidMechanism.speciesIndex("C3H8")].hf =
        std::numeric_limits<double>::quiet_NaN();
    check("nonfinite formation enthalpy is rejected",
        constructionRejected(invalidMechanism, "invalid formation enthalpy"));

    invalidMechanism = mechanism;
    invalidMechanism.species[invalidMechanism.speciesIndex("N2")].name =
        "not-in-registry";
    check("unknown mechanism species is rejected",
        constructionRejected(invalidMechanism, "not present in SpeciesRegistry"));

    invalidMechanism = mechanism;
    invalidMechanism.species[invalidMechanism.speciesIndex("N2")].name = "CO2";
    check("duplicate mechanism species is rejected",
        constructionRejected(invalidMechanism, "duplicate species"));

    invalidMechanism = mechanism;
    invalidMechanism.reactions[1].DH = 0.0;
    check("reaction heat release must match authoritative formation enthalpies",
        constructionRejected(invalidMechanism, "reaction enthalpy mismatch"));

    invalidMechanism = mechanism;
    invalidMechanism.reactions[1].stoich["O2"] = -0.4;
    invalidMechanism.reactions[1].concExp["O2"] = 0.4;
    check("elementally unbalanced reaction is rejected",
        constructionRejected(invalidMechanism, "element balance mismatch"));

    invalidMechanism = mechanism;
    invalidMechanism.reactions[1].A = 0.0;
    check("zero forward pre-exponential is rejected",
        constructionRejected(invalidMechanism, "pre-exponential"));

    const std::string mechanismText = readTextFile("tests/twostep_co.yaml");
    check("omitted forward pre-exponential field is rejected",
        readerRejected(
            replaceOnce(mechanismText, "    A: 2.0e12", ""),
            "missing required field 'A'"));
    check("unknown reaction key is rejected",
        readerRejected(
            replaceOnce(
                mechanismText,
                "    A: 2.0e12",
                "    preexponential: 2.0e12"),
            "unknown reaction key 'preexponential'"));

    invalidMechanism = mechanism;
    for (auto& item : invalidMechanism.reactions[0].stoich)
        item.second *= 4.0e307;
    for (auto& item : invalidMechanism.reactions[0].concExp)
        item.second *= 4.0e307;
    check("irreversible reaction aggregate overflow is rejected",
        constructionRejected(invalidMechanism, "reaction aggregate"));

    invalidMechanism = mechanism;
    for (auto& item : invalidMechanism.reactions[1].stoich)
        item.second *= 4.0e307;
    for (auto& item : invalidMechanism.reactions[1].concExp)
        item.second *= 4.0e307;
    check("reversible reaction aggregate overflow is rejected",
        constructionRejected(invalidMechanism, "reaction aggregate"));

    std::cout << "\nRESULT: " << passCount << " passed, "
        << failCount << " failed\n";
    std::cout << (failCount == 0
        ? "TWOSTEP_CO_CONTRACT_PASS\n"
        : "TWOSTEP_CO_CONTRACT_FAIL\n");
    return failCount == 0 ? 0 : 1;
}
