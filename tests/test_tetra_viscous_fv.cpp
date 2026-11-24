#include "TetraFiniteVolume.H"
#include "CompleteCellState.H"
#include "CoolPropThermoModel.H"
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace Foam;
int main(int argc,char**argv)
{
    if(argc!=2)return 2;
    if(!thermo::coolPropThermoAvailable())return 77;
    try
    {
        TetraMesh mesh;mesh.read(argv[1]);
        auto model=thermo::makeCoolPropThermoModel();
        PrimitiveVars w(2);w.p=4.11e6;w.T=400.;w.Y={.784,.216};
        w.U={4.,2.,11.};
        const auto q=model->initializePTY(w.p,w.T,w.Y);
        if(!q.success())throw std::runtime_error("test EOS state failed");
        TetraFiniteVolume fv(mesh,{"O2","C3H8"},
            std::unique_ptr<FluxScheme>(new AUSMFluxAdapter({"O2","C3H8"})));
        for(std::size_t c=0;c<mesh.cellCount();++c)fv.setPrimitive(c,w,*q.state);
        const TetraBoundaryProvider boundary=[&](std::size_t,std::size_t,
            PrimitiveVars& ghost,thermo::EquilibriumState& eq)
        {ghost=w;eq=*q.state;};
        std::vector<ConservedVars> baseline,after;
        std::vector<double> inviscidDt,viscousDt;
        fv.computeResidual(boundary,baseline);
        fv.localTimeStep(.5,inviscidDt);
        for(std::size_t c=0;c<mesh.cellCount();++c)
            if(fv.molecularSpectral(c)!=0.)
                throw std::runtime_error("default inviscid path has molecular spectral diffusion");
        bool rejection=false;std::size_t calls=0;
        fv.setMolecularTransportProvider([&](const PrimitiveVars&,
                                               const thermo::EquilibriumState& qstate)
        {
            ++calls;
            if(rejection)throw std::runtime_error("transport rejected");
            if(qstate.region==thermo::EquilibriumRegion::TWO_PHASE)
                throw std::runtime_error("two-phase transport not admitted");
            return std::make_pair(1.8e-5,0.03);
        });
        fv.computeResidual(boundary,after);
        fv.localTimeStep(.5,viscousDt);
        if(calls!=1)throw std::runtime_error("identical cell states were reflashed");
        for(std::size_t c=0;c<mesh.cellCount();++c)
            if(!(fv.molecularSpectral(c)>0.)||!(viscousDt[c]<inviscidDt[c]))
                throw std::runtime_error("molecular spectral did not reduce local pseudo-time");
        double totalEnergy=0.,work=0.,maxMass=0.,maxSpecies=0.,maxEnergy=0.;
        for(std::size_t c=0;c<mesh.cellCount();++c)
        {
            maxMass=std::max(maxMass,std::abs(after[c].rho-baseline[c].rho));
            maxEnergy=std::max(maxEnergy,std::abs(after[c].rhoE-baseline[c].rhoE));
            for(std::size_t s=0;s<2;++s)
                maxSpecies=std::max(maxSpecies,std::abs(after[c].rhoY[s]-baseline[c].rhoY[s]));
            for(int d=0;d<3;++d)
                work+=w.U[d]*(after[c].rhoU[d]-baseline[c].rhoU[d]);
            totalEnergy+=after[c].rhoE-baseline[c].rhoE;
        }
        if(maxMass!=0.||maxSpecies!=0.||!(work>0.)||std::abs(totalEnergy)>1e-9
           ||!std::isfinite(maxEnergy))
            throw std::runtime_error("FV molecular residual sign/conservation failed");
        rejection=true;bool rejected=false;
        auto modified=w;modified.T+=1.;
        auto modifiedQ=model->initializePTY(modified.p,modified.T,modified.Y);
        if(!modifiedQ.success())throw std::runtime_error("cache rejection state failed");
        fv.setPrimitive(0,modified,*modifiedQ.state);
        try{fv.computeResidual(boundary,after);}catch(const std::runtime_error&){rejected=true;}
        if(!rejected)throw std::runtime_error("invalid molecular transport did not fail closed");
        std::cout<<"TETRA_VISCOUS_FV_PASS cells="<<mesh.cellCount()<<'\n';
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
