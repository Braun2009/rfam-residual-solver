#ifndef SectorNonReflectingOutlet_H
#define SectorNonReflectingOutlet_H

#include "Domain.H"
#include "SectorMeshMetrics.H"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Foam
{

class SectorNonReflectingOutlet
{
    const SectorMeshMetrics& mesh_;
    double pInf_,sigma_;

    static bool supportsCharacteristicState(
        const thermo::EquilibriumState& state)
    {
        return state.valid()
            && state.acoustic.closure
                == thermo::AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE;
    }

public:
    SectorNonReflectingOutlet(const SectorMeshMetrics& mesh,double pInf,
                              double sigma=0.25)
        :mesh_(mesh),pInf_(pInf),sigma_(sigma) {}

    void apply(Domain& d,double dt) const
    {
        if(!d.stateReady()||d.tripped()||d.stepTransactionActive())
            throw std::logic_error(
                "SectorNonReflectingOutlet requires a ready Domain outside a step transaction");
        if(!d.thermoModel().supportsPrimitivePTYInitialization())
            throw std::runtime_error(
                "SectorNonReflectingOutlet phase-aware characteristic closure is not implemented");
        const int ib=d.Nx()-1,ii=d.Nx()-2;
        const double dx=mesh_.xFace(ib)-mesh_.xFace(ii);

        struct BoundaryCandidate
        {
            PrimitiveVars primitive;
            thermo::EquilibriumState equilibrium;
            bool ready=false;
            std::string error;

            explicit BoundaryCandidate(int nSpecies)
            : primitive(nSpecies)
            {}
        };
        const int radialCount=std::max(0,d.Ny()-2);
        const int azimuthalCount=std::max(0,d.Nz()-2);
        std::vector<BoundaryCandidate> candidates;
        candidates.reserve(
            static_cast<std::size_t>(radialCount*azimuthalCount));
        for(int n=0;n<radialCount*azimuthalCount;++n)
            candidates.emplace_back(d.nSpecies());
        const auto candidateIndex=[&](int j,int k)
        {
            return static_cast<std::size_t>(
                (j-1)*azimuthalCount+(k-1));
        };

        auto computeCandidate=[&](int j,int k)
        {
            BoundaryCandidate& candidate=candidates[candidateIndex(j,k)];
            auto initializeState=[&](
                const PrimitiveVars& primitive,
                thermo::EquilibriumState& equilibrium)
            {
                thermo::EquilibriumResult result=
                    d.thermoModel().initializePTY(
                        primitive.p,primitive.T,primitive.Y);
                if(!result.success()
                   ||!supportsCharacteristicState(*result.state))
                {
                    candidate.error=
                        "PTY recovery failed or returned an unsupported acoustic closure";
                    return false;
                }
                equilibrium=std::move(*result.state);
                return true;
            };

            const PrimitiveVars wb=d.W(ib,j,k);
            const PrimitiveVars wi=d.W(ii,j,k);
            const thermo::EquilibriumState& qb=d.Q(ib,j,k);
            const thermo::EquilibriumState& qi=d.Q(ii,j,k);
            if(!supportsCharacteristicState(qb)
               ||!supportsCharacteristicState(qi))
            {
                candidate.error=
                    "frozen-composition single-phase characteristic state required";
                return;
            }
            const double rho=qb.rho;
            const double a=qb.acoustic.soundSpeed();
            const double u=wb.U[0];

            if(u>=a)
            {
                candidate.primitive=wi;
                candidate.equilibrium=qi;
                candidate.ready=true;
                return;
            }
            if(u<=0.0)
            {
                PrimitiveVars farField=wi;
                farField.p=pInf_;
                farField.U[0]=std::max(0.0,wi.U[0]);
                candidate.primitive=farField;
                candidate.ready=initializeState(
                    candidate.primitive,candidate.equilibrium);
                return;
            }

            const double dpdx=(wb.p-wi.p)/dx;
            const double dudx=(wb.U[0]-wi.U[0])/dx;
            const double drhodx=(rho-qi.rho)/dx;
            const double l5=(u+a)*(dpdx+rho*a*dudx);
            const double l2=u*(a*a*drhodx-dpdx);
            const double mach=u/a;
            const double relax=sigma_*(1.0-mach*mach)*a/mesh_.length();
            const double l1=relax*(wb.p-pInf_);

            const double pNew=std::max(1e3,wb.p-0.5*dt*(l5+l1));
            const double uNew=wb.U[0]-dt*(l5-l1)/(2.0*rho*a);
            const double rhoNew=std::max(
                1e-6,rho-dt*(l2+0.5*(l5+l1))/(a*a));
            PrimitiveVars out=wb;
            out.p=pNew;
            out.U[0]=uNew;

            auto densityError=[&](double temperature,double& error)
            {
                PrimitiveVars trial=out;
                trial.T=temperature;
                thermo::EquilibriumState equilibrium;
                if(!initializeState(trial,equilibrium)) return false;
                error=equilibrium.rho-rhoNew;
                return true;
            };
            double lo=50.0,hi=6000.0,flo=0.0,fhi=0.0;
            if(!densityError(lo,flo)||!densityError(hi,fhi)) return;
            if(flo*fhi<0.0)
            {
                double temperature=std::clamp(wb.T,lo,hi);
                for(int it=0;it<12;++it)
                {
                    double residual=0.0;
                    if(!densityError(temperature,residual)) return;
                    if(residual>0.0) lo=temperature; else hi=temperature;
                    if(std::abs(residual)<1e-9*std::max(rhoNew,1.0)) break;
                    const double delta=std::max(0.01,1e-3*temperature);
                    const double upper=std::min(temperature+delta,hi);
                    const double lower=std::max(temperature-delta,lo);
                    if(!(upper>lower)) break;
                    double upperError=0.0,lowerError=0.0;
                    if(!densityError(upper,upperError)
                       ||!densityError(lower,lowerError)) return;
                    const double derivative=(upperError-lowerError)/(upper-lower);
                    double next=temperature-residual/derivative;
                    if(!std::isfinite(next)||next<=lo||next>=hi)
                        next=0.5*(lo+hi);
                    if(std::abs(next-temperature)<1e-4)
                    {
                        temperature=next;
                        break;
                    }
                    temperature=next;
                }
                out.T=temperature;
            }
            candidate.primitive=out;
            candidate.ready=initializeState(
                candidate.primitive,candidate.equilibrium);
        };

        #pragma omp parallel for collapse(2) schedule(static)
        for(int j=1;j<d.Ny()-1;++j) for(int k=1;k<d.Nz()-1;++k)
            computeCandidate(j,k);

        for(int j=1;j<d.Ny()-1;++j) for(int k=1;k<d.Nz()-1;++k)
        {
            BoundaryCandidate& candidate=candidates[candidateIndex(j,k)];
            if(!candidate.ready)
                throw std::runtime_error(
                    "SectorNonReflectingOutlet candidate recovery failed at ("
                    +std::to_string(j)+","+std::to_string(k)+"): "
                    +candidate.error);
        }
        for(int j=1;j<d.Ny()-1;++j) for(int k=1;k<d.Nz()-1;++k)
        {
            BoundaryCandidate& candidate=candidates[candidateIndex(j,k)];
            if(!d.setCell(
                    ib,j,k,candidate.primitive,candidate.equilibrium))
                throw std::runtime_error(
                    "SectorNonReflectingOutlet failed to publish a complete state: "
                    +d.failure().detail);
        }
    }

    // Steady/pseudo-time outlet with no physical-time characteristic update.
    // pressureContinuation=0 gives pure extrapolation for startup. A positive
    // value advances the checkpointed boundary pressure geometrically toward
    // p_inf while velocity/temperature/composition come from the adjacent
    // active cell. Supersonic outflow remains pure extrapolation.
    void applySteadyImpl(
        Domain& d,double pressureContinuation,
        bool inStepTransaction) const
    {
        if(!d.stateReady()||d.tripped()
            ||d.stepTransactionActive()!=inStepTransaction)
            throw std::logic_error(
                inStepTransaction
                    ?"SectorNonReflectingOutlet steady inner closure requires an active transaction"
                    :"SectorNonReflectingOutlet steady closure requires a ready Domain outside a step transaction");
        if(!d.thermoModel().supportsPrimitivePTYInitialization())
            throw std::runtime_error(
                "SectorNonReflectingOutlet steady phase-aware pressure "
                "closure is not implemented");
        pressureContinuation=std::clamp(
            pressureContinuation,0.0,1.0);
        const int ib=d.Nx()-1,ii=d.Nx()-2;
        std::vector<Domain::CompleteCellReplacement> replacements;
        replacements.reserve(static_cast<std::size_t>(
            std::max(0,d.Ny()-2)*std::max(0,d.Nz()-2)));
        for(int j=1;j<d.Ny()-1;++j)
        for(int k=1;k<d.Nz()-1;++k)
        {
            PrimitiveVars out=d.W(ii,j,k);
            const thermo::EquilibriumState& interior=d.Q(ii,j,k);
            if(out.U[0]<interior.acoustic.soundSpeed())
            {
                if(pressureContinuation>0.0)
                {
                    const double previousPressure=d.W(ib,j,k).p;
                    out.p=previousPressure+pressureContinuation
                        *(pInf_-previousPressure);
                }
                out.U[0]=std::max(0.0,out.U[0]);
                const thermo::EquilibriumResult initialized=
                    d.thermoModel().initializePTY(out.p,out.T,out.Y);
                if(!initialized.success()
                   ||!supportsCharacteristicState(*initialized.state))
                    throw std::runtime_error(
                        "SectorNonReflectingOutlet steady PTY recovery failed");
                replacements.push_back(Domain::CompleteCellReplacement{
                    ib,j,k,out,*initialized.state});
            }
            else
                replacements.push_back(Domain::CompleteCellReplacement{
                    ib,j,k,out,interior});
        }
        const bool published=inStepTransaction
            ?d.replaceCompleteCellsInStepTransaction(replacements)
            :d.replaceCompleteCells(replacements);
        if(!published)
            throw std::runtime_error(
                "SectorNonReflectingOutlet steady publication failed: "
                +d.failure().detail);
    }

    void applySteady(
        Domain& d,double pressureContinuation=0.0) const
    { applySteadyImpl(d,pressureContinuation,false); }

    void applySteadyInStepTransaction(
        Domain& d,double pressureContinuation=0.0) const
    { applySteadyImpl(d,pressureContinuation,true); }

    double farFieldPressure() const { return pInf_; }
};

} // End namespace Foam
#endif
