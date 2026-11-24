#include "GlobalReaction1Step.H"
#include <iostream>
#include <iomanip>

namespace Foam
{

GlobalReaction1Step::GlobalReaction1Step(
    double A, double n, double Ea,
    double a, double b
)
:
    rate_(A, n, Ea),
    a_C3H8_(a),
    b_O2_(b)
{
    nSpecies_ = 5;
    nReactions_ = 1;
    
    speciesNames_ = {"C3H8", "O2", "CO2", "H2O", "N2"};
    molecularWeights_ = {MW_C3H8, MW_O2, MW_CO2, MW_H2O, MW_N2};
    
    iC3H8_ = 0;
    iO2_   = 1;
    iCO2_  = 2;
    iH2O_  = 3;
    iN2_   = 4;
    
    nu_C3H8_ = -1.0;
    nu_O2_   = -5.0;
    nu_CO2_  = +3.0;
    nu_H2O_  = +4.0;
}

void GlobalReaction1Step::computeSourceTerms(
    double T,
    double p,
    double rho,
    const std::vector<double>& Y,
    std::vector<double>& omega_i,
    double& qdot
) const
{
    omega_i.resize(nSpecies_, 0.0);
    qdot = 0.0;
    
    double C_C3H8 = rho * Y[iC3H8_] / MW_C3H8;
    double C_O2   = rho * Y[iO2_]   / MW_O2;
    
    if (C_C3H8 < 1e-20 || C_O2 < 1e-20) return;
    
    double k = rate_.rate(T);
    double r_fuel = k * std::pow(C_C3H8, a_C3H8_) * std::pow(C_O2, b_O2_);
    
    double omega_mol_C3H8 = nu_C3H8_ * r_fuel;
    double omega_mol_O2   = nu_O2_   * r_fuel;
    double omega_mol_CO2  = nu_CO2_  * r_fuel;
    double omega_mol_H2O  = nu_H2O_  * r_fuel;
    
    omega_i[iC3H8_] = omega_mol_C3H8 * MW_C3H8;
    omega_i[iO2_]   = omega_mol_O2   * MW_O2;
    omega_i[iCO2_]  = omega_mol_CO2  * MW_CO2;
    omega_i[iH2O_]  = omega_mol_H2O  * MW_H2O;
    omega_i[iN2_]   = 0.0;
    
    qdot = r_fuel * (-DELTA_H);
}

void GlobalReaction1Step::testConservation(
    const std::vector<double>& omega_i,
    double qdot
) const
{
    using namespace std;
    
    cout << "\n=== Conservation Tests ===" << endl;
    
    // Test 1: Mass conservation (relative)
    double sum_omega = 0.0;
    for (int i = 0; i < nSpecies_; i++)
    {
        sum_omega += omega_i[i];
    }
    
    double rel_mass_error = abs(sum_omega) / (abs(omega_i[iC3H8_]) + 1e-20);
    
    cout << "1. Mass conservation:" << endl;
    cout << "   Σ omega_i = " << scientific << sum_omega << " kg/m³/s" << endl;
    cout << "   Relative error = " << rel_mass_error << endl;
    cout << "   Status: " << (rel_mass_error < 1e-4 ? "✅ PASS" : "❌ FAIL") << endl;
    
    // Test 2: Element conservation (relative)
    double r_C3H8 = omega_i[iC3H8_] / MW_C3H8;
    double r_O2   = omega_i[iO2_]   / MW_O2;
    double r_CO2  = omega_i[iCO2_]  / MW_CO2;
    double r_H2O  = omega_i[iH2O_]  / MW_H2O;
    
    double C_balance = -3.0 * r_C3H8 + 3.0 * r_CO2;
    double H_balance = -8.0 * r_C3H8 + 8.0 * r_H2O;
    double O_balance = -10.0 * r_O2 + 6.0 * r_CO2 + 4.0 * r_H2O;
    
    double r_fuel = abs(r_C3H8);
    double rel_C_error = abs(C_balance) / (3.0 * r_fuel + 1e-20);
    double rel_H_error = abs(H_balance) / (8.0 * r_fuel + 1e-20);
    double rel_O_error = abs(O_balance) / (10.0 * r_fuel + 1e-20);
    
    cout << "\n2. Element conservation:" << endl;
    cout << "   C balance: " << scientific << C_balance << " mol/m³/s (rel: " << rel_C_error << ")" << endl;
    cout << "   H balance: " << H_balance << " mol/m³/s (rel: " << rel_H_error << ")" << endl;
    cout << "   O balance: " << O_balance << " mol/m³/s (rel: " << rel_O_error << ")" << endl;
    
    bool elem_ok = rel_C_error < 1e-4 && rel_H_error < 1e-4 && rel_O_error < 1e-4;
    cout << "   Status: " << (elem_ok ? "✅ PASS" : "❌ FAIL") << endl;
    
    // Test 3: Heat release consistency
    if (abs(omega_i[iC3H8_]) > 1e-20)
    {
        double omega_C3H8_mass = -omega_i[iC3H8_];
        double specific_heat_release = qdot / omega_C3H8_mass;
        double LHV_C3H8 = 46.4e6;
        double error = abs(specific_heat_release - LHV_C3H8) / LHV_C3H8;
        
        cout << "\n3. Heat release consistency:" << endl;
        cout << "   qdot / (-omega_C3H8) = " << fixed << setprecision(2) 
             << specific_heat_release/1e6 << " MJ/kg" << endl;
        cout << "   Expected LHV         = " << LHV_C3H8/1e6 << " MJ/kg" << endl;
        cout << "   Error                = " << setprecision(1) << error*100 << "%" << endl;
        cout << "   Status: " << (error < 0.01 ? "✅ PASS" : "❌ FAIL") << endl;
    }
    
    cout << "==========================\n" << endl;
}

} // End namespace Foam
