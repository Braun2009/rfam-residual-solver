// Nonuniform states expose neighbour-index and flux-cancellation defects which
// a uniform freestream cannot detect. Also checks real momentum advancement
// from rest and exact recovery of the updated conservative state.
#include "TetraFiniteVolume.H"
#include "TetraLUSGS.H"
#include "TetraSteadyMonitor.H"
#include <fstream>
#include <iostream>
#include <cstdio>
#include <filesystem>
#include <chrono>
using namespace Foam;
int main()
{
    const auto path=std::filesystem::temp_directory_path()/
        ("rfam-two-tetra-"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count())+".rfamtet");
    try
    {
        {std::ofstream out(path,std::ios::binary);auto put=[&](const auto& x){out.write(reinterpret_cast<const char*>(&x),sizeof(x));};
        out.write("RFAMTET1",8);std::uint64_t np=5,nc=2,nf=7;put(np);put(nc);put(nf);
        double points[15]={0,0,0, 1,0,0, 0,1,0, 0,0,1, 0,0,-1};put(points);
        std::int32_t cells[8]={0,1,2,3, 0,2,1,4};put(cells);double vols[2]={1./6,1./6};put(vols);
        std::int32_t faces[42]={0,2,1,0,1,0, 0,1,3,0,-1,0, 1,2,3,0,-1,0, 2,0,3,0,-1,0,
                              0,2,4,1,-1,0, 2,1,4,1,-1,0, 1,0,4,1,-1,0};put(faces);}
        TetraMesh mesh;mesh.read(path.string());std::vector<std::string> names{"O2","C3H8"};
        auto model=thermo::makeLegacySinglePhaseThermoModel(names);
        TetraFiniteVolume fv(mesh,names,std::make_unique<AUSMFluxAdapter>(names));
        PrimitiveVars w(2);w.p=5e5;w.T=450.;w.Y={.6,.4};w.U={30.,-10.,50.};
        auto q=model->initializePTY(w.p,w.T,w.Y);fv.setPrimitive(0,w,*q.state);
        w.p=7e5;w.T=520.;w.Y={.2,.8};w.U={-15.,35.,-20.};q=model->initializePTY(w.p,w.T,w.Y);fv.setPrimitive(1,w,*q.state);
        const TetraBoundaryProvider extrapolate=[&](std::size_t c,std::size_t,PrimitiveVars& ghost,thermo::EquilibriumState& thermo)
        {ghost=fv.W(c);thermo=fv.Q(c);};
        std::vector<ConservedVars> residual;fv.computeResidual(extrapolate,residual);
        auto values=[](const ConservedVars& u){std::vector<double> x{u.rho,u.rhoE,u.rhoZ,u.rhoC};x.insert(x.end(),u.rhoU.begin(),u.rhoU.end());x.insert(x.end(),u.rhoY.begin(),u.rhoY.end());return x;};
        std::vector<double> total(9),boundary(9),magnitude(9);
        for(auto& r:residual){auto x=values(r);for(int i=0;i<9;++i)total[i]+=x[i];}
        for(std::size_t f=1;f<mesh.faceCount();++f){auto x=values(fv.lastFaceFlux(f));for(int i=0;i<9;++i){boundary[i]+=x[i]*fv.faceArea(f);magnitude[i]+=std::abs(x[i])*fv.faceArea(f);}}
        for(int i=0;i<9;++i)if(std::abs(total[i]-boundary[i])>1e-12*std::max(magnitude[i],1.))throw std::runtime_error("shared face does not cancel");
        // A prescribed physical flux must be used identically in the residual
        // and reported boundary mass flow, while interior fluxes still cancel.
        const auto baseline=residual;
        const double oldFlux=fv.lastFaceFlux(1).rho;
        fv.setBoundaryFluxProvider([](std::size_t,std::size_t f,ConservedVars& flux)
        {
            if(f!=1)return false;
            flux=ConservedVars(2);flux.rho=-2.;flux.rhoE=3e5;
            flux.rhoY={-.4,-1.6};return true;
        });
        fv.computeResidual(extrapolate,residual);
        if(std::abs(residual[0].rho-baseline[0].rho-(-2.-oldFlux)*fv.faceArea(1))>1e-10
            ||values(residual[1])!=values(baseline[1]))throw std::runtime_error("prescribed flux residual mismatch");
        if(std::abs(fv.faceMassFlow(1,fv.W(0),fv.W(0),fv.Q(0),fv.Q(0))+2.*fv.faceArea(1))>1e-14)
            throw std::runtime_error("prescribed boundary mass monitor mismatch");
        fv.setBoundaryFluxProvider({});
        // Same composition, unequal pressure, zero initial velocity.
        for(int c=0;c<2;++c){w.p=(5.+c*.1)*1e5;w.T=450.;w.Y={.6,.4};w.U={0.,0.,0.};q=model->initializePTY(w.p,w.T,w.Y);fv.setPrimitive(c,w,*q.state);}
        fv.computeResidual(extrapolate,residual);
        TetraLUSGS stepper(fv);stepper.step(extrapolate,.05,1.,.5,*model);
        double speed=0.;for(int c=0;c<2;++c){for(double v:fv.W(c).U)speed+=v*v;auto state=CompleteStateConversion::fromConserved(fv.U(c),fv.W(c),*model);if(!state.success())throw std::runtime_error("updated conserved state cannot recover");}
        if(speed<1e-12)throw std::runtime_error("LU-SGS freezes initially stationary momentum");
        const auto& d=stepper.diagonal();
        for(int c=0;c<2;++c)if(!std::isfinite(d[c])||d[c]<=0)throw std::runtime_error("invalid implicit diagonal");
        // Failed recovery must not publish a partly advanced field.
        auto saved0=values(fv.U(0)),saved1=values(fv.U(1));stepper.setGuard(1000,1100,1e4,2e7);
        bool rejected=false;try{stepper.step(extrapolate,.05,1.,.5,*model);}catch(const std::runtime_error&){rejected=true;}
        if(!rejected||saved0!=values(fv.U(0))||saved1!=values(fv.U(1)))throw std::runtime_error("failed update was not atomic");
        stepper.setGuard(80.,4000.,1e4,2e7);
        rejected=false;bool sawTrial=false;
        try
        {
            stepper.step(extrapolate,.05,1.,.5,*model,[&]()
            {
                sawTrial=values(fv.U(0))!=saved0;
                throw std::runtime_error("deliberate boundary-state rejection");
            });
        }
        catch(const std::runtime_error&){rejected=true;}
        if(!rejected||!sawTrial||saved0!=values(fv.U(0))||saved1!=values(fv.U(1)))
            throw std::runtime_error("boundary failure did not roll back complete cell update");
        // Primitive/thermo boundary mismatch must fail, never use cached fallback.
        const TetraBoundaryProvider bad=[&](std::size_t c,std::size_t,PrimitiveVars& ghost,thermo::EquilibriumState& thermo)
        {ghost=fv.W(c);ghost.p*=2.;thermo=fv.Q(c);};
        TetraKernelFailure failure;fv.computeResidual(bad,residual,&failure);
        if(!failure.failed||!residual.empty())throw std::runtime_error("inconsistent boundary state was accepted");
        const TetraBoundaryProvider nonfinite=[&](std::size_t c,std::size_t,PrimitiveVars& ghost,thermo::EquilibriumState& thermo)
        {ghost=fv.W(c);ghost.T=std::numeric_limits<double>::quiet_NaN();thermo=fv.Q(c);};
        fv.computeResidual(nonfinite,residual,&failure);
        if(!failure.failed||!residual.empty())throw std::runtime_error("NaN cached boundary temperature accepted");
        // Global cancellation alone must not pass the steady gate. Start with
        // an actual uniform moving field, then inject opposing local defects.
        w.p=5e5;w.T=450.;w.Y={.6,.4};w.U={3.,5.,7.};
        q=model->initializePTY(w.p,w.T,w.Y);
        for(int c=0;c<2;++c)fv.setPrimitive(c,w,*q.state);
        fv.computeResidual(extrapolate,residual);
        const auto uniformResidual=residual;
        TetraSteadyMonitor monitor(1.,300.,3e5,.6,.4);
        if(monitor.evaluate(fv,residual,0.,0,0).converged)
            throw std::runtime_error("initial field marked steady without iterations");
        for(int iteration=1;iteration<=5;++iteration)
        {
            const auto status=monitor.evaluate(fv,residual,0.,0,iteration);
            if(status.converged!=(iteration==5))throw std::runtime_error("steady sample window incorrect");
        }
        residual[0].rhoE+=3e5;residual[1].rhoE-=3e5;
        auto status=monitor.evaluate(fv,residual,0.,0,6);
        if(status.converged||status.stableSamples!=0||status.local[4]<1.9)
            throw std::runtime_error("cancelling local energy defects marked converged");
        residual=uniformResidual;
        residual[0].rhoY[0]+=.1;residual[0].rhoY[1]-=.1;
        residual[1].rhoY[0]-=.1;residual[1].rhoY[1]+=.1;
        status=monitor.evaluate(fv,residual,0.,0,7);
        if(status.converged||status.stableSamples!=0||status.local[5]<.3)
            throw std::runtime_error("opposing species defects hidden by total mass balance");
        residual=uniformResidual;
        if(monitor.evaluate(fv,residual,.01,0,8).stableSamples!=0
            ||monitor.evaluate(fv,residual,0.,1,9).stableSamples!=0)
            throw std::runtime_error("changing or limited field admitted as steady");
        bool repeatedRejected=false;
        try{monitor.evaluate(fv,residual,0.,0,9);}catch(const std::invalid_argument&){repeatedRejected=true;}
        if(!repeatedRejected)throw std::runtime_error("duplicate monitor sample accepted");
        std::vector<double> cellDt,boundaryDt;
        fv.localTimeStep(.5,cellDt);
        const TetraBoundaryProvider fastBoundary=[&](std::size_t c,std::size_t,PrimitiveVars& ghost,thermo::EquilibriumState& thermo)
        {ghost=fv.W(c);ghost.U[0]=1500.;thermo=fv.Q(c);};
        fv.localTimeStep(.5,boundaryDt,fastBoundary);
        if(!(boundaryDt[0]<cellDt[0]&&boundaryDt[1]<cellDt[1]))
            throw std::runtime_error("physical CFL ignores faster boundary trace");
        // A mass-only inner gate would accept this nonzero energy equation.
        w.U={0.,0.,0.};
        for(int c=0;c<2;++c)fv.setPrimitive(c,w,*q.state);
        fv.setBoundaryFluxProvider([&](std::size_t,std::size_t f,ConservedVars& flux)
        {
            flux=ConservedVars(2);
            const auto& n=fv.faceNormal(f);
            for(int d=0;d<3;++d)flux.rhoU[d]=w.p*n[d];
            if(f==1)flux.rhoE=5000.;
            return true;
        });
        std::vector<ConservedVars> reference{fv.U(0),fv.U(1)};
        TetraLUSGS physical(fv);
        const auto energyDefect=physical.measurePhysicalResidual(extrapolate,reference,1e-5);
        if(!(energyDefect.first>0.&&energyDefect.second>0.)
           ||std::abs(physical.residual()[0].rho)>1e-12
           ||std::abs(physical.residual()[1].rho)>1e-12)
            throw std::runtime_error("energy-only physical residual hidden by mass norm");
        const auto physicalReport=physical.step(extrapolate,.5,1.,.5,*model,{},1.,1.,&reference,1e-5);
        if(!(physicalReport.residualNorm>0.&&physicalReport.residualMaxNorm>0.)
           ||!(physical.diagonal()[0]>mesh.volume(0)/1e-5))
            throw std::runtime_error("physical-time diagonal or residual not assembled");
        const auto accepted0=values(fv.U(0)),accepted1=values(fv.U(1));
        physical.setGuard(1000.,1100.,1e4,2e7);
        rejected=false;
        try{physical.step(extrapolate,.5,1.,.5,*model,{},1.,1.,&reference,1e-5);}
        catch(const std::runtime_error&){rejected=true;}
        if(!rejected||accepted0!=values(fv.U(0))||accepted1!=values(fv.U(1)))
            throw std::runtime_error("rejected physical inner iteration changed accepted field");
        std::filesystem::remove(path);std::cout<<"TETRA_CONSERVATION_PASS speed2="<<speed<<"\n";return 0;
    }
    catch(const std::exception& e){std::filesystem::remove(path);std::cerr<<e.what()<<"\n";return 1;}
}
