#ifndef GlobalReaction1Step_H
#define GlobalReaction1Step_H

#include "ChemistryModel.H"
#include "ArrheniusRate.H"
#include <cmath>
#include <algorithm>

namespace Foam
{

/*---------------------------------------------------------------------------*\
    GlobalReaction1Step - One-step global reaction for C3H8/O2
    
    Reaction:
        C3H8 + 5 O2 → 3 CO2 + 4 H2O
    
    Rate:
        r = k(T) * [C3H8]^a * [O2]^b   [mol_fuel/m³/s]
        k(T) = A * T^n * exp(-Ea/RT)
    
    Species source terms (stoichiometry):
        omega_mol_C3H8 = -1 * r
        omega_mol_O2   = -5 * r
        omega_mol_CO2  = +3 * r
        omega_mol_H2O  = +4 * r
    
    Mass source terms:
        omega_i = omega_mol_i * MW_i   [kg/m³/s]
    
    Heat release:
        qdot = r * (-ΔH)   [W/m³]
        ΔH = -2.046e6 J/mol_fuel (exothermic)
\*---------------------------------------------------------------------------*/

class GlobalReaction1Step : public ChemistryModel
{
private:
    ArrheniusRate rate_;
    
    // Species indices
    int iC3H8_, iO2_, iCO2_, iH2O_, iN2_;
    
    // Stoichiometric coefficients (positive for products)
    double nu_C3H8_, nu_O2_, nu_CO2_, nu_H2O_;
    
    // Molecular weights [kg/mol]
    static constexpr double MW_C3H8 = 0.04409562;
    static constexpr double MW_O2   = 0.0319988;
    static constexpr double MW_CO2  = 0.0440095;
    static constexpr double MW_H2O  = 0.01801528;
    static constexpr double MW_N2   = 0.0280134;
    
    // Heat of reaction [J/mol_fuel]
    // LHV_C3H8 = 46.4 MJ/kg * 0.044 kg/mol = 2.046 MJ/mol
    static constexpr double DELTA_H = -2.046e6;  // Negative = exothermic
    
    // Reaction order exponents
    double a_C3H8_;  // Fuel exponent
    double b_O2_;    // Oxidizer exponent

public:
    GlobalReaction1Step(
        double A = 8.6e11,      // Pre-exponential [depends on order]
        double n = 0.0,         // Temperature exponent
        double Ea = 125520.0,   // Activation energy [J/mol]
        double a = 1.0,         // Fuel order
        double b = 1.5          // Oxidizer order
    );
    
    void computeSourceTerms(
        double T,
        double p,
        double rho,
        const std::vector<double>& Y,
        std::vector<double>& omega_i,
        double& qdot
    ) const override;
    
    // Test utilities
    void testConservation(const std::vector<double>& omega_i, double qdot) const;
};

} // End namespace Foam

#endif
