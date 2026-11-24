/*---------------------------------------------------------------------------*\
  realFluidAllMach - Conserved Variables Structure (FPV Extended)
\*---------------------------------------------------------------------------*/

#ifndef ConservedVars_H
#define ConservedVars_H

#include <vector>
#include <string>

namespace Foam
{
class RealFluidState;

// Conserved variables: rho, rhoU, rhoE, rhoY, rhoZ, rhoC
struct ConservedVars
{
    double rho;                     // Density [kg/m³]
    std::vector<double> rhoU;       // Momentum [kg/(m²·s)]
    double rhoE;                    // Total energy [J/m³]
    std::vector<double> rhoY;       // Species densities [kg/m³]
    
    // FPV variables
    double rhoZ;                    // Mixture fraction density [kg/m³]
    double rhoC;                    // Progress variable density [kg/m³]

    ConservedVars() : rho(0.0), rhoE(0.0), rhoZ(0.0), rhoC(0.0)
    {
        rhoU.resize(3, 0.0);
        rhoY.resize(2, 0.0);
    }

    ConservedVars(int nSpecies) : rho(0.0), rhoE(0.0), rhoZ(0.0), rhoC(0.0)
    {
        rhoU.resize(3, 0.0);
        rhoY.resize(nSpecies, 0.0);
    }
};

// Primitive variables: p, T, U, Y, Z, C
struct PrimitiveVars
{
    double p;                       // Pressure [Pa]
    double T;                       // Temperature [K]
    std::vector<double> U;          // Velocity [m/s]
    std::vector<double> Y;          // Mass fractions [-]
    
    // FPV variables
    double Z;                       // Mixture fraction [-]
    double C;                       // Progress variable [-]

    PrimitiveVars() : p(101325.0), T(300.0), Z(0.0), C(0.0)
    {
        U.resize(3, 0.0);
        Y.resize(2, 0.0);
        Y[0] = 1.0;
    }

    PrimitiveVars(int nSpecies) : p(101325.0), T(300.0), Z(0.0), C(0.0)
    {
        U.resize(3, 0.0);
        Y.resize(nSpecies, 0.0);
    }
};

// State conversion functions
class StateConversion
{
public:
    // species: nullptr => LEGACY default EOS species path ({O2,C3H8}); the
    //   2-species MVP/precond reback stays byte-for-byte unchanged.
    // non-null => use this EXACT species order, which MUST match cons.rhoY /
    //   prim.Y (the global SpeciesRegistry order for reacting flow). The
    //   function never guesses species identity — the caller (Domain) owns it.

    // Convert conserved to primitive using RealFluidState
    static bool conservedToPrimitive(
        const ConservedVars& cons,
        PrimitiveVars& prim,
        const std::vector<std::string>* species = nullptr
    );
    static bool conservedToPrimitive(
        const ConservedVars& cons,PrimitiveVars& prim,
        RealFluidState& state,double initialTemperature);

    // Convert primitive to conserved using RealFluidState
    static bool primitiveToConserved(
        const PrimitiveVars& prim,
        ConservedVars& cons,
        const std::vector<std::string>* species = nullptr
    );
    static bool primitiveToConserved(
        const PrimitiveVars& prim,ConservedVars& cons,RealFluidState& state);
};

} // End namespace Foam

#endif
