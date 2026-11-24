#include "TwoStepCO.H"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>

namespace Foam
{

TwoStepCO::TwoStepCO(
    double A1, double n1, double Ea1,
    double A2, double n2, double Ea2,
    double a1, double b1
)
:
    rate1_(A1, n1, Ea1),
    rate2_(A2, n2, Ea2),
    a1_(a1),
    b1_(b1)
{
    nSpecies_ = 6;
    nReactions_ = 2;

    speciesNames_    = {"C3H8", "O2", "CO", "CO2", "H2O", "N2"};
    molecularWeights_ = {MW_C3H8, MW_O2, MW_CO, MW_CO2, MW_H2O, MW_N2};

    iC3H8_ = 0;
    iO2_   = 1;
    iCO_   = 2;
    iCO2_  = 3;
    iH2O_  = 4;
    iN2_   = 5;
}

TwoStepCO::ReactionRates TwoStepCO::reactionRates(
    double T,
    double rho,
    const std::vector<double>& Y
) const
{
    if (!std::isfinite(T) || T <= 0.0
        || !std::isfinite(rho) || rho <= 0.0)
        throw std::invalid_argument(
            "TwoStepCO::reactionRates requires finite positive T and rho");
    if (Y.size() != static_cast<std::size_t>(nSpecies_))
        throw std::invalid_argument(
            "TwoStepCO::reactionRates species dimension mismatch");
    for (double value : Y)
        if (!std::isfinite(value) || value < 0.0)
            throw std::invalid_argument(
                "TwoStepCO::reactionRates requires finite nonnegative mass fractions");

    // Molar concentrations [mol/m³]
    const double C_C3H8 = rho * Y[iC3H8_] / MW_C3H8;
    const double C_O2   = rho * Y[iO2_]   / MW_O2;
    const double C_CO   = rho * Y[iCO_]   / MW_CO;
    const double C_CO2  = rho * Y[iCO2_]  / MW_CO2;

    ReactionRates rates;

    // R1: C3H8 + 3.5 O2 → 3 CO + 4 H2O
    // r1 [mol_C3H8/m³/s]
    if (C_C3H8 > 1e-20 && C_O2 > 1e-20)
        rates.r1Forward = rate1_.rate(T)
            * regularizedConcentrationPower(C_C3H8, a1_)
            * regularizedConcentrationPower(C_O2, b1_);

    // R2 forward: CO + 0.5 O2 → CO2,  kf
    // R2 reverse: CO2 → CO + 0.5 O2,  kr = kf / Kc
    // Net rate r2_net = kf * C_CO * C_O2^0.5 - kr * C_CO2.
    // Matching the forward concentration exponents to the reaction
    // stoichiometry is required for kr=kf/Kc to satisfy detailed balance.
    const double kf = rate2_.rate(T);
    const double Kc = Kc_R2(T);
    const double kr = kf / Kc;
    if (C_CO > 0.0 && C_O2 > 0.0)
        rates.r2Forward = kf * C_CO * std::sqrt(C_O2);
    rates.r2Reverse = kr * C_CO2;
    return rates;
}

void TwoStepCO::computeSourceTerms(
    double T, double p, double rho,
    const std::vector<double>& Y,
    std::vector<double>& omega_i,
    double& qdot
) const
{
    (void)p;
    omega_i.assign(nSpecies_, 0.0);

    const ReactionRates rates = reactionRates(T, rho, Y);
    const double r1 = rates.r1Forward;
    const double r2 = rates.r2Net();

    // Mass source terms [kg/m³/s]
    // R1 stoichiometry: -1 C3H8, -3.5 O2, +3 CO, +4 H2O
    omega_i[iC3H8_] += -1.0 * r1 * MW_C3H8;
    omega_i[iO2_]   += -3.5 * r1 * MW_O2;
    omega_i[iCO_]   += +3.0 * r1 * MW_CO;
    omega_i[iH2O_]  += +4.0 * r1 * MW_H2O;

    // R2 stoichiometry: -1 CO, -0.5 O2, +1 CO2
    omega_i[iO2_]   += -0.5 * r2 * MW_O2;
    omega_i[iCO_]   += -1.0 * r2 * MW_CO;
    omega_i[iCO2_]  += +1.0 * r2 * MW_CO2;

    // Heat release [W/m³]
    qdot = r1 * (-DH1) + r2 * (-DH2);
}

void TwoStepCO::print() const
{
    std::cout << "2-Step CO Mechanism:" << std::endl;
    std::cout << "  R1: C3H8 + 3.5 O2 -> 3 CO + 4 H2O" << std::endl;
    std::cout << "  R2: CO + 0.5 O2 <=> CO2" << std::endl;
    std::cout << "  Species: " << nSpecies_ << ", Reactions: " << nReactions_ << std::endl;
}

} // End namespace Foam
