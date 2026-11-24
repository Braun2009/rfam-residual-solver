#ifndef TetraBinaryEulerJacobian_H
#define TetraBinaryEulerJacobian_H
#include "FrozenEulerFluxJacobian.H"
#include "TetraAcousticAdmission.H"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Foam
{
// Diagnostic binary extension of the frozen-composition Euler block. This is
// NOT a derivative of the AUSM numerical flux. Pressure derivatives are built
// at the same stable single-phase branch, with Y1 = 1-Y0. No phase projection.
// At fixed p,T:
//   0 = p_rho|e,Y * rho_Y|p,T + p_e|rho,Y * e_Y|p,T + p_Y|rho,e.
// The missing term is p_Y*(d(rhoY0)-Y0*d(rho))/rho. Including it restores
// independence of species-dependent reference enthalpy offsets.
class TetraBinaryEulerJacobian
{
    static bool singlePhase(const thermo::EquilibriumState& q)
    {
        return tetraAdmitsAcoustics(q)&&q.singlePhase.has_value()
            &&q.thermodynamics&&q.thermodynamics->valid()
            &&q.acoustic.closure==thermo::AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE;
    }
public:
    struct Coefficient
    {
        double dpdY_rhoE=0.;
        double relativeRefinementError=0.;
    };
    static Coefficient prepare(const thermo::EquilibriumState& q,
                               const thermo::EquilibriumThermoModel& model)
    {
        if(q.overallMassFractions.size()!=2||!singlePhase(q))
            throw std::runtime_error("binary Jacobian needs a valid single-phase derivative payload");
        const double y=q.overallMassFractions[0];const auto& td=*q.thermodynamics;
        const double pe=td.dpdT_rho/td.dedT_rho;
        const double pr=td.dpdrho_T-pe*td.dedrho_T;
        const auto derivative=[&](double step)
        {
            const double lower=std::max(0.,y-step),upper=std::min(1.,y+step);
            const auto sample=[&](double composition)
            {
                if(composition==y)return q;
                const auto r=model.initializePTY(q.p,q.T,{composition,1.-composition});
                if(!r.success()||!singlePhase(*r.state))
                    throw std::runtime_error("binary Jacobian composition perturbation failed/crossed phase boundary: "+r.diagnostics.backendStatus);
                return *r.state;
            };
            const auto lo=sample(lower),hi=sample(upper);
            return -pr*(hi.rho-lo.rho)/(upper-lower)-pe*(hi.e-lo.e)/(upper-lower);
        };
        const double coarse=derivative(2e-5),fine=derivative(1e-5);
        const double scale=std::max({std::abs(coarse),std::abs(fine),q.p,1.});
        const double error=std::abs(fine-coarse)/scale;
        if(!std::isfinite(fine)||!std::isfinite(error)||error>1e-4)
            throw std::runtime_error("binary Jacobian derivative did not converge under refinement");
        return {fine,error};
    }
    static bool apply(const PrimitiveVars& w,const thermo::EquilibriumState& q,
                      const ConservedVars& d,const std::array<double,3>& n,
                      const Coefficient& coefficient,ConservedVars& product)
    {
        if(w.Y.size()!=2||d.rhoY.size()!=2||!singlePhase(q)
           ||!std::isfinite(coefficient.dpdY_rhoE)
           ||!FrozenEulerFluxJacobian::apply(w,q,d,n,product))return false;
        const double dy=(d.rhoY[0]-w.Y[0]*d.rho)/q.rho;
        const double dp=coefficient.dpdY_rhoE*dy;
        double un=0.;for(int i=0;i<3;++i){product.rhoU[i]+=dp*n[i];un+=w.U[i]*n[i];}
        product.rhoE+=dp*un;
        return std::isfinite(dp)&&std::isfinite(product.rhoE)
            &&std::all_of(product.rhoU.begin(),product.rhoU.end(),[](double v){return std::isfinite(v);});
    }
};
}
#endif
