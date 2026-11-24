#ifndef TwoStepCO_H
#define TwoStepCO_H

#include "ChemistryModel.H"
#include "ArrheniusRate.H"
#include "MassActionThermo.H"
#include <cmath>

namespace Foam
{

/*---------------------------------------------------------------------------*\
    TwoStepCO - 2-step mechanism for C3H8/O2

    R1: C3H8 + 3.5 O2  → 3 CO + 4 H2O   (fuel combustion)
    R2: CO   + 0.5 O2  → CO2              (CO oxidation)

    References:
    - Westbrook & Dryer (1981)
    - Jones & Lindstedt (1988)

    R1 retains the global empirical concentration orders. R2 uses a
    stoichiometric reversible mass-action closure so that kr=kf/Kc satisfies
    detailed balance. Its A2/Ea2 values remain a candidate kinetic closure
    pending independent rate validation; this contract does not enable any
    production capability.
\*---------------------------------------------------------------------------*/

class TwoStepCO : public ChemistryModel
{
private:
    ArrheniusRate rate1_;   // R1: C3H8 → CO
    ArrheniusRate rate2_;   // R2: CO   → CO2
    double a1_, b1_;        // R1 concentration exponents (calibrated)

    // Species indices
    int iC3H8_, iO2_, iCO_, iCO2_, iH2O_, iN2_;

    // MW [kg/mol]
    static constexpr double MW_C3H8 = 0.04409562;
    static constexpr double MW_O2   = 0.0319988;
    static constexpr double MW_CO   = 0.0280101;
    static constexpr double MW_CO2  = 0.0440095;
    static constexpr double MW_H2O  = 0.01801528;
    static constexpr double MW_N2   = 0.0280134;

    // Heat of reaction [J/mol_fuel or J/mol_CO]
    // R1: C3H8 + 3.5 O2 → 3 CO + 4 H2O
    //   DH1 = 3*(-110.53e3) + 4*(-241.82e3) - (-103.85e3) = -1195.02e3 J/mol_C3H8
    //   (fixed 2026-06-22: was -1.5498e6, inconsistent with the formation
    //    enthalpies above -> R1 heat release ~30% too high. Caught by the
    //    Phase 6 CVODE energy-consistency test [test_phase6_cvode Test 2].)
    static constexpr double DH1 = -1.19502e6;

    // R2: CO + 0.5 O2 → CO2
    //   DH2 = -393.51e3 - (-110.53e3) = -282.98e3 J/mol_CO
    static constexpr double DH2 = -2.8298e5;

    // Equilibrium constant Kc for R2 (CO + 0.5O2 <=> CO2) at temperature T.
    // Kp = exp(-ΔG°/(RT)) = exp((ΔS° - ΔH°/T)/R)
    // Using NASA-7 absolute enthalpies and entropies (GRI-Mech 3.0 coefficients
    // already in NASA7Thermo, consistent with the rest of the RFAM thermo core).
    // Kc = Kp / (RT/p°)^Δn   where Δn = 1 - 1 - 0.5 = -0.5  (product - reactant moles)
    // Standard pressure p° = 101325 Pa.
    //
    // This makes R2 REVERSIBLE, which enables CO2 dissociation at high T and
    // presses the A-box terminal temperature toward the physical adiabatic value
    // (~3800K) rather than the no-dissociation limit (~7000K).
    static double Kc_R2(double T)
    {
        static const std::vector<std::string> species = {"CO", "O2", "CO2"};
        static const std::vector<double> stoich = {-1.0, -0.5, 1.0};
        return concentrationEquilibriumConstant(T, species, stoich);
    }

public:
    struct ReactionRates
    {
        double r1Forward = 0.0;  // mol/m3/s, C3H8 oxidation progress
        double r2Forward = 0.0;  // mol/m3/s, CO oxidation progress
        double r2Reverse = 0.0;  // mol/m3/s, CO2 dissociation progress

        double r2Net() const { return r2Forward - r2Reverse; }
    };

    // Concentration exponents (2026-08-07 GRI-30 calibration):
    //   a1 = 1.303 (C3H8), b1 = 0.263 (O2) replace the Westbrook-Dryer
    //   atmospheric-air values (0.1 / 1.65) whose O2^1.65 explodes in the
    //   rich-oxygen chamber (qdot 1e14-1e16 vs physical 1e10-1e12).
    TwoStepCO(
        double A1  = 7.77e5,     double n1  = 0.0, double Ea1 = 125520.0,
        double A2  = 1.55e8,     double n2  = 0.0, double Ea2 = 167360.0,
        double a1  = 1.303,      double b1  = 0.263
    );

    // Exact progress rates used by computeSourceTerms.  This diagnostic API
    // exists for same-mechanism verification; it is not solver-owned state.
    ReactionRates reactionRates(
        double T,
        double rho,
        const std::vector<double>& Y
    ) const;

    void computeSourceTerms(
        double T, double p, double rho,
        const std::vector<double>& Y,
        std::vector<double>& omega_i,
        double& qdot
    ) const override;

    void print() const override;
};

} // End namespace Foam

#endif
