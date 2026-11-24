/*---------------------------------------------------------------------------*\
    Parameterized wrapper for the offline RFAM R1/R2 0-D runner.

    The fixed runner remains byte-stable. This translation unit substitutes a
    parameterized chemistry object into that exact reactor/integrator source
    and forwards the original trajectory arguments unchanged.
\*---------------------------------------------------------------------------*/

#include "FiniteRate/TwoStepCO.H"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Foam
{
namespace
{

struct R2FitParameters
{
    double r1A = 1.0e10;
    double r1Ea = 125520.0;
    double r1FuelOrder = 0.1;
    double r1OxygenOrder = 1.65;
    double r2A = 2.0e12;
    double r2Ea = 167360.0;
};

R2FitParameters fitParameters;

class R2FitTwoStepCO : public ChemistryModel
{
    ArrheniusRate rate1_;
    ArrheniusRate rate2_;
    double r1FuelOrder_;
    double r1OxygenOrder_;

    static constexpr double MW_C3H8 = 0.04409562;
    static constexpr double MW_O2 = 0.0319988;
    static constexpr double MW_CO = 0.0280101;
    static constexpr double MW_CO2 = 0.0440095;
    static constexpr double MW_H2O = 0.01801528;
    static constexpr double MW_N2 = 0.0280134;
    static constexpr double DH1 = -1.19502e6;
    static constexpr double DH2 = -2.8298e5;

    static double r2EquilibriumConstant(double temperature)
    {
        static const std::vector<std::string> species = {"CO", "O2", "CO2"};
        static const std::vector<double> stoich = {-1.0, -0.5, 1.0};
        return concentrationEquilibriumConstant(temperature, species, stoich);
    }

public:
    R2FitTwoStepCO()
    :
        rate1_(fitParameters.r1A, 0.0, fitParameters.r1Ea),
        rate2_(fitParameters.r2A, 0.0, fitParameters.r2Ea),
        r1FuelOrder_(fitParameters.r1FuelOrder),
        r1OxygenOrder_(fitParameters.r1OxygenOrder)
    {
        nSpecies_ = 6;
        nReactions_ = 2;
        speciesNames_ = {"C3H8", "O2", "CO", "CO2", "H2O", "N2"};
        molecularWeights_ = {
            MW_C3H8, MW_O2, MW_CO, MW_CO2, MW_H2O, MW_N2};
    }

    void computeSourceTerms(
        double temperature,
        double pressure,
        double density,
        const std::vector<double>& massFractions,
        std::vector<double>& sources,
        double& heatRelease) const override
    {
        (void)pressure;
        if (!std::isfinite(temperature) || temperature <= 0.0
            || !std::isfinite(density) || density <= 0.0
            || massFractions.size() != 6U)
            throw std::invalid_argument("invalid-parameterized-chemistry-state");
        for (double value : massFractions)
            if (!std::isfinite(value) || value < 0.0)
                throw std::invalid_argument(
                    "invalid-parameterized-chemistry-composition");

        const double fuel = density * massFractions[0] / MW_C3H8;
        const double oxygen = density * massFractions[1] / MW_O2;
        const double carbonMonoxide = density * massFractions[2] / MW_CO;
        const double carbonDioxide = density * massFractions[3] / MW_CO2;

        double r1 = 0.0;
        if (fuel > 1.0e-20 && oxygen > 1.0e-20)
            r1 = rate1_.rate(temperature)
                * regularizedConcentrationPower(fuel, r1FuelOrder_)
                * regularizedConcentrationPower(oxygen, r1OxygenOrder_);

        const double r2ForwardRate = rate2_.rate(temperature);
        const double r2ReverseRate = r2ForwardRate
            / r2EquilibriumConstant(temperature);
        double r2Forward = 0.0;
        if (carbonMonoxide > 0.0 && oxygen > 0.0)
            r2Forward = r2ForwardRate * carbonMonoxide * std::sqrt(oxygen);
        const double r2 = r2Forward - r2ReverseRate * carbonDioxide;

        sources.assign(6U, 0.0);
        sources[0] -= r1 * MW_C3H8;
        sources[1] -= 3.5 * r1 * MW_O2;
        sources[2] += 3.0 * r1 * MW_CO;
        sources[4] += 4.0 * r1 * MW_H2O;
        sources[1] -= 0.5 * r2 * MW_O2;
        sources[2] -= r2 * MW_CO;
        sources[3] += r2 * MW_CO2;
        heatRelease = r1 * (-DH1) + r2 * (-DH2);
    }
};

double parseFitParameter(
    const char* text,
    const char* name,
    double lower,
    double upper)
{
    std::size_t consumed = 0U;
    double value = 0.0;
    try
    {
        value = std::stod(text, &consumed);
    }
    catch (const std::exception&)
    {
        throw std::invalid_argument(std::string("invalid-fit-parameter:") + name);
    }
    if (consumed != std::string(text).size() || !std::isfinite(value)
        || value < lower || value > upper)
        throw std::invalid_argument(std::string("invalid-fit-parameter:") + name);
    return value;
}

} // End anonymous namespace
} // End namespace Foam

#define TwoStepCO R2FitTwoStepCO
#define main rfamR2FixedRunnerMain
#include "r2_0d_rfam_runner.cpp"
#undef main
#undef TwoStepCO

int main(int argc, char** argv)
{
    try
    {
        if (argc < 20)
            throw std::invalid_argument("fit-runner-argument-count");

        Foam::fitParameters.r1A = Foam::parseFitParameter(
            argv[1], "R1_A_SI", 1.0e5, 1.0e20);
        Foam::fitParameters.r1Ea = Foam::parseFitParameter(
            argv[2], "R1_Ea_J_per_mol", 0.0, 3.0e5);
        Foam::fitParameters.r1FuelOrder = Foam::parseFitParameter(
            argv[3], "R1_order_C3H8", 0.0, 2.0);
        Foam::fitParameters.r1OxygenOrder = Foam::parseFitParameter(
            argv[4], "R1_order_O2", 0.0, 2.5);
        Foam::fitParameters.r2A = Foam::parseFitParameter(
            argv[5], "R2_A_SI", 1.0e5, 1.0e20);
        Foam::fitParameters.r2Ea = Foam::parseFitParameter(
            argv[6], "R2_Ea_J_per_mol", 0.0, 3.0e5);

        std::cout << std::setprecision(17)
                  << "fit_parameter_schema=rfam.r2-joint-fit-parameters/v1\n"
                  << "fit_R1_A_SI=" << Foam::fitParameters.r1A << '\n'
                  << "fit_R1_Ea_J_per_mol=" << Foam::fitParameters.r1Ea << '\n'
                  << "fit_R1_order_C3H8="
                  << Foam::fitParameters.r1FuelOrder << '\n'
                  << "fit_R1_order_O2="
                  << Foam::fitParameters.r1OxygenOrder << '\n'
                  << "fit_R2_A_SI=" << Foam::fitParameters.r2A << '\n'
                  << "fit_R2_Ea_J_per_mol=" << Foam::fitParameters.r2Ea << '\n';

        std::vector<char*> forwarded;
        forwarded.reserve(static_cast<std::size_t>(argc - 6));
        forwarded.push_back(argv[0]);
        for (int index = 7; index < argc; ++index)
            forwarded.push_back(argv[index]);
        return rfamR2FixedRunnerMain(
            static_cast<int>(forwarded.size()), forwarded.data());
    }
    catch (const std::exception& error)
    {
        std::cerr << "R2_FIT_RUNNER_ERROR " << error.what() << '\n';
        return 2;
    }
}
