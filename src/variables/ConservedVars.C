/*---------------------------------------------------------------------------*\
  realFluidAllMach - State Conversion Implementation (FPV Extended)
\*---------------------------------------------------------------------------*/

#include "ConservedVars.H"
#include "RealFluidState.H"
#include <cmath>
#include <algorithm>
#include <utility>

namespace Foam
{

bool StateConversion::conservedToPrimitive(
    const ConservedVars& cons,
    PrimitiveVars& prim,
    const std::vector<std::string>* species
)
{
    RealFluidState state;
    if(species) state.setSpecies(*species);
    return conservedToPrimitive(cons,prim,state,species?prim.T:300.0);
}

bool StateConversion::conservedToPrimitive(
    const ConservedVars& cons,PrimitiveVars& prim,
    RealFluidState& state,double initialTemperature)
{
    if (!std::isfinite(cons.rho) || !std::isfinite(cons.rhoE) ||
        cons.rho <= 1e-10 || cons.rhoU.size() != 3 ||
        cons.rhoY.size() != prim.Y.size() || prim.U.size() != 3)
        return false;

    PrimitiveVars candidate = prim;

    // Extract velocity
    for (size_t i = 0; i < 3; i++)
    {
        if (!std::isfinite(cons.rhoU[i])) return false;
        candidate.U[i] = cons.rhoU[i] / cons.rho;
        if (!std::isfinite(candidate.U[i])) return false;
    }

    // Extract species mass fractions
    double speciesMass = 0.0;
    for (size_t i = 0; i < cons.rhoY.size(); i++)
    {
        if (!std::isfinite(cons.rhoY[i]) || cons.rhoY[i] < 0.0)
            return false;
        speciesMass += cons.rhoY[i];
        candidate.Y[i] = cons.rhoY[i] / cons.rho;
    }
    if (!std::isfinite(speciesMass) ||
        std::abs(speciesMass - cons.rho) > 1e-6*cons.rho)
        return false;

    // Extract FPV variables
    if (!std::isfinite(cons.rhoZ) || !std::isfinite(cons.rhoC)) return false;
    candidate.Z = cons.rhoZ / cons.rho;
    candidate.C = cons.rhoC / cons.rho;
    if (!std::isfinite(candidate.Z) || candidate.Z < 0.0 || candidate.Z > 1.0
        || !std::isfinite(candidate.C) || candidate.C < 0.0
        || candidate.C > 1.0)
        return false;

    // Compute internal energy
    double Ek = 0.0;
    for (size_t i = 0; i < 3; i++)
    {
        Ek += candidate.U[i] * candidate.U[i];
    }
    Ek *= 0.5;

    double e_internal = cons.rhoE / cons.rho - Ek;
    if (!std::isfinite(e_internal)) return false;

    // Use RealFluidState to get (p, T) from (rho, e, Y)
    auto status = state.updateRhoEY(
        cons.rho,e_internal,candidate.Y,initialTemperature);

    if (status != RealFluidState::Status::SUCCESS ||
        !std::isfinite(state.p()) || !std::isfinite(state.T()))
        return false;

    candidate.p = state.p();
    candidate.T = state.T();
    prim = std::move(candidate);
    return true;
}

bool StateConversion::primitiveToConserved(
    const PrimitiveVars& prim,
    ConservedVars& cons,
    const std::vector<std::string>* species
)
{
    RealFluidState state;
    if(species) state.setSpecies(*species);
    return primitiveToConserved(prim,cons,state);
}

bool StateConversion::primitiveToConserved(
    const PrimitiveVars& prim,ConservedVars& cons,RealFluidState& state)
{
    if (!std::isfinite(prim.p) || !std::isfinite(prim.T) ||
        prim.p <= 0.0 || prim.T <= 0.0 || prim.U.size() != 3 ||
        prim.Y.empty() || !std::isfinite(prim.Z) || prim.Z < 0.0
        || prim.Z > 1.0 || !std::isfinite(prim.C) || prim.C < 0.0
        || prim.C > 1.0)
        return false;

    double sumY = 0.0;
    for (double y : prim.Y)
    {
        if (!std::isfinite(y) || y < 0.0) return false;
        sumY += y;
    }
    if (!std::isfinite(sumY) || std::abs(sumY - 1.0) > 1e-6)
        return false;

    if (state.updatePTY(prim.p, prim.T, prim.Y)
        != RealFluidState::Status::SUCCESS)
        return false;

    ConservedVars candidate(static_cast<int>(prim.Y.size()));

    candidate.rho = state.rho();
    if (!std::isfinite(candidate.rho) || candidate.rho <= 0.0 ||
        !std::isfinite(state.e()))
        return false;

    // Momentum
    for (size_t i = 0; i < 3; i++)
    {
        if (!std::isfinite(prim.U[i])) return false;
        candidate.rhoU[i] = candidate.rho * prim.U[i];
    }

    // Kinetic energy
    double Ek = 0.0;
    for (size_t i = 0; i < 3; i++)
    {
        Ek += prim.U[i] * prim.U[i];
    }
    Ek *= 0.5;

    // Total energy
    candidate.rhoE = candidate.rho * (state.e() + Ek);
    if (!std::isfinite(candidate.rhoE)) return false;

    // Species
    for (size_t i = 0; i < prim.Y.size(); i++)
    {
        candidate.rhoY[i] = candidate.rho * prim.Y[i];
    }

    candidate.rhoZ = candidate.rho * prim.Z;
    candidate.rhoC = candidate.rho * prim.C;
    if (!std::isfinite(candidate.rhoZ) || !std::isfinite(candidate.rhoC))
        return false;

    cons = std::move(candidate);
    return true;
}

} // End namespace Foam
