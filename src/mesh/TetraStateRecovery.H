#ifndef TetraStateRecovery_H
#define TetraStateRecovery_H

#include "CompleteCellState.H"
#include "TetraAcousticAdmission.H"
#include "StateAcceptance.H"
#include "ThermoCounters.H"
#include <algorithm>
#include <cmath>

namespace Foam
{
// Tetra solver recovery is stricter than the legacy diagnostic HEM contract:
// no mass/energy projection, no significant composition clipping, no acoustic
// fallback. Only floating-point composition noise (<1e-12*rho) is repaired.
inline CompleteCellConversionResult recoverTetraState(
    const ConservedVars& input,const PrimitiveVars& previous,
    const thermo::EquilibriumState& previousQ,
    const thermo::EquilibriumThermoModel& model)
{
    auto fail=[](const char* text)
    {
        CompleteCellConversionResult r;
        r.status=thermo::EquilibriumStatus::RESIDUAL_CONTRACT_FAILED;
        r.diagnostics.backendStatus=text;return r;
    };
    if(!std::isfinite(input.rho)||input.rho<=1e-10)
        return fail("TETRA_INVALID_DENSITY");
    ConservedVars repaired=input;
    double sum=0.;std::size_t largest=0;
    double repairedMass=0.;
    const double repairTolerance=speciesRepairRel();
    for(std::size_t s=0;s<repaired.rhoY.size();++s)
    {
        double& value=repaired.rhoY[s];
        if(!std::isfinite(value)||value < -repairTolerance*input.rho
           ||value>(1.+repairTolerance)*input.rho)
            return fail("TETRA_SPECIES_OUT_OF_BOUNDS");
        if(value<0.){repairedMass+=-value;}
        else if(value>input.rho){repairedMass+=value-input.rho;}
        value=std::clamp(value,0.,input.rho);sum+=value;
        if(value>repaired.rhoY[largest])largest=s;
    }
    if(repairedMass>0.)
    {
        // Mass-conservative positivity fix-up: the deficit is added back to the
        // largest species, so the cell keeps its total mass. The repaired mass
        // is recorded so the size of the transient can be audited per step.
        thermo::thermoCount(thermo::thermoCounters().speciesRepairEvents);
        thermo::thermoMax(thermo::thermoCounters().maxSpeciesRepairRel, repairedMass/input.rho);
    }
    if(repaired.rhoY.empty()||std::abs(sum-input.rho)>repairTolerance*input.rho)
        return fail("TETRA_SPECIES_MASS_MISMATCH");
    repaired.rhoY[largest]+=input.rho-sum;
    for(double* value:{&repaired.rhoZ,&repaired.rhoC})
    {
        if(!std::isfinite(*value)||*value < -1e-12*input.rho||*value>(1.+1e-12)*input.rho)
            return fail("TETRA_SCALAR_OUT_OF_BOUNDS");
        *value=std::clamp(*value,0.,input.rho);
    }
    thermo::EquilibriumSeed seed;
    seed.temperature=previousQ.T;seed.previousTemperature=previousQ.T;
    seed.pressure=previousQ.p;seed.previousRegion=previousQ.region;
    seed.previousVaporMassFraction=previousQ.vaporMassFraction;
    auto result=CompleteStateConversion::fromConserved(repaired,previous,model,&seed);
    if(!result.success())return result;
    const auto& q=result.state->Q;
    double kinetic=0.;for(double m:input.rhoU)kinetic+=.5*(m/input.rho)*(m/input.rho);
    const double e=input.rhoE/input.rho-kinetic;
    const bool closureOk=
        std::abs(q.rho-input.rho)<=stateAcceptanceRel()*std::max(input.rho,1.)
        &&std::abs(q.e-e)<=stateAcceptanceRel()*std::max({std::abs(e),std::abs(q.e),1.});
    if(q.conservedStateProjected||!tetraAdmitsAcoustics(q)
       ||(!closureOk&&!projectionAcceptEnabled()))
        return fail("TETRA_PROJECTED_OR_UNADMITTED_STATE");
    return result;
}
}
#endif
