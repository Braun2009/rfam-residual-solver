#include "TetraPressureOutlet.H"
#include "CoolPropThermoModel.H"
#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>
using namespace Foam;
int main()
{
    if(!thermo::coolPropThermoAvailable())return 77;
    try
    {
        auto model=thermo::makeCoolPropThermoModel();
        PrimitiveVars wi(2),wo(2);wi.p=4e6;wi.T=400.;wi.Y={.73,.27};wi.U={7.,-9.,20.};
        auto initial=model->initializePTY(wi.p,wi.T,wi.Y);
        if(!initial.success())throw std::runtime_error("initial test state failed");
        const auto qi=*initial.state;thermo::EquilibriumState qo;
        const std::vector<double> n{0.,0.,1.};
        TetraPressureOutlet matched(*model,wi.p,300.,{1.,0.});matched.apply(wi,qi,n,wo,qo);
        if(std::abs(wo.U[2]-wi.U[2])>1e-8||std::abs(wo.T-wi.T)>1e-7||std::abs(qo.rho-qi.rho)>1e-7)
            throw std::runtime_error("matched pressure outlet changes freestream");
        // Compare every conserved flux component to the unmodified PSY route.
        auto local=wi;local.p=qi.p;
        for(double pressure:{std::nextafter(qi.p,std::numeric_limits<double>::infinity()),
                             std::nextafter(qi.p,0.)})
        {
            TetraPressureOutlet full(*model,pressure,400.,{.73,.27});
            TetraPressureOutlet fast(*model,pressure,400.,{.73,.27},true);
            PrimitiveVars exactW(2),fastW(2);
            thermo::EquilibriumState exactQ,fastQ;
            full.apply(local,qi,n,exactW,exactQ);
            fast.apply(local,qi,n,fastW,fastQ);
            const auto exactFlux=TetraPressureOutlet::physicalFlux(exactW,exactQ,n);
            const auto fastFlux=TetraPressureOutlet::physicalFlux(fastW,fastQ,n);
            const auto close=[](double a,double b)
            {return std::abs(a-b)<=2e-10*std::max({1.,std::abs(a),std::abs(b)});};
            if(!close(exactFlux.rho,fastFlux.rho)
               ||!close(exactFlux.rhoE,fastFlux.rhoE)
               ||!close(exactFlux.rhoU[0],fastFlux.rhoU[0])
               ||!close(exactFlux.rhoU[1],fastFlux.rhoU[1])
               ||!close(exactFlux.rhoU[2],fastFlux.rhoU[2])
               ||!close(exactFlux.rhoY[0],fastFlux.rhoY[0])
               ||!close(exactFlux.rhoY[1],fastFlux.rhoY[1]))
                throw std::runtime_error("one-ULP outlet shortcut changed conserved flux");
            if(fast.matchHits()!=1)
                throw std::runtime_error("one-ULP outlet path was not exercised");
        }
        {
            // The cold A outlet has a negative roundoff-scale normal speed
            // and an owner pressure one ULP above the prescribed pressure.
            TetraPressureOutlet full(*model,std::nextafter(qi.p,0.),400.,{.73,.27});
            TetraPressureOutlet fastNear(*model,std::nextafter(qi.p,0.),400.,{.73,.27},true);
            local.U[2]=-4e-16;
            PrimitiveVars exactW(2),fastW(2);
            thermo::EquilibriumState exactQ,fastQ;
            full.apply(local,qi,n,exactW,exactQ);
            fastNear.apply(local,qi,n,fastW,fastQ);
            const auto exactFlux=TetraPressureOutlet::physicalFlux(exactW,exactQ,n);
            const auto fastFlux=TetraPressureOutlet::physicalFlux(fastW,fastQ,n);
            if(fastNear.matchHits()!=1||exactFlux.rho*fastFlux.rho<0.
               ||std::abs(exactFlux.rho-fastFlux.rho)>1e-10
               ||std::abs(exactFlux.rhoE-fastFlux.rhoE)>1e-8
               ||std::abs(exactFlux.rhoU[2]-fastFlux.rhoU[2])>1e-7
               ||std::abs(exactFlux.rhoY[0]-fastFlux.rhoY[0])>1e-10
               ||std::abs(exactFlux.rhoY[1]-fastFlux.rhoY[1])>1e-10)
                throw std::runtime_error("near-zero one-ULP outlet flux differs from full PSY");
            local.U[2]=20.;
        }
        {
            const double pressure=std::nextafter(qi.p,std::numeric_limits<double>::infinity());
            TetraPressureOutlet full(*model,pressure,310.,{0.,1.});
            TetraPressureOutlet fastReverse(*model,pressure,310.,{0.,1.},true);
            local.U[2]=-1e-10;
            PrimitiveVars exactW(2),fastW(2);
            thermo::EquilibriumState exactQ,fastQ;
            full.apply(local,qi,n,exactW,exactQ);
            fastReverse.apply(local,qi,n,fastW,fastQ);
            const auto exactFlux=TetraPressureOutlet::physicalFlux(exactW,exactQ,n);
            const auto fastFlux=TetraPressureOutlet::physicalFlux(fastW,fastQ,n);
            if(fastReverse.matchHits()!=1||exactFlux.rho>=0.||fastFlux.rho>=0.
               ||std::abs(exactFlux.rho-fastFlux.rho)>1e-10
               ||std::abs(exactFlux.rhoE-fastFlux.rhoE)>1e-8
               ||std::abs(exactFlux.rhoU[2]-fastFlux.rhoU[2])>1e-7
               ||std::abs(exactFlux.rhoY[0]-fastFlux.rhoY[0])>1e-10
               ||std::abs(exactFlux.rhoY[1]-fastFlux.rhoY[1])>1e-10)
                throw std::runtime_error("one-ULP reservoir backflow differs from full PSY");
            local.U[2]=20.;
        }
        TetraPressureOutlet fast(*model,qi.p+2.,400.,{.73,.27},true);
        fast.apply(local,qi,n,wo,qo);
        local.U[2]=-1.;
        TetraPressureOutlet reverse(*model,qi.p,400.,{.73,.27},true);
        reverse.apply(local,qi,n,wo,qo);
        if(fast.matchHits()!=0||reverse.matchHits()!=0)
            throw std::runtime_error("nonroundoff or reverse flow used shortcut");
        TetraPressureOutlet expanded(*model,3.8e6,300.,{1.,0.});expanded.apply(wi,qi,n,wo,qo);
        if(std::abs(qo.p-3.8e6)>1.||std::abs(qo.s-qi.s)>1e-7||wo.Y!=wi.Y
           ||wo.U[0]!=wi.U[0]||wo.U[1]!=wi.U[1]||wo.U[2]<=wi.U[2])
            throw std::runtime_error("isentropic pressure outlet contract failed");
        auto cell=CompleteStateConversion::fromPrimitiveAndEquilibrium(wo,qo,*model);
        if(!cell.success())throw std::runtime_error("pressure outlet W/Q inconsistent");
        auto flux=TetraPressureOutlet::physicalFlux(wo,qo,n);
        if(std::abs(flux.rho-qo.rho*wo.U[2])>1e-10
           ||std::abs(flux.rhoU[2]-flux.rho*wo.U[2]-3.8e6)>1e-5
           ||std::abs(flux.rhoY[0]+flux.rhoY[1]-flux.rho)>1e-10)
            throw std::runtime_error("pressure trace flux is not conservative");
        const double energyFlux=flux.rho*(qo.h+.5*(wo.U[0]*wo.U[0]+wo.U[1]*wo.U[1]+wo.U[2]*wo.U[2]));
        if(std::abs(flux.rhoE-energyFlux)>1e-8*std::max(std::abs(energyFlux),1.))
            throw std::runtime_error("outlet energy flux mismatch");
        auto mismatch=wo;mismatch.T+=1.;bool mismatchRejected=false;
        try{TetraPressureOutlet::physicalFlux(mismatch,qo,n);}catch(const std::runtime_error&){mismatchRejected=true;}
        if(!mismatchRejected)throw std::runtime_error("mismatched outlet trace was accepted");
        // Python CoolProp 7.2.0 + scipy.integrate.quad of dp/(rho*a),
        // at fixed entropy from O2 10 kPa / 300 K, down to 9 kPa.
        wi.p=1e4;wi.T=300.;wi.Y={1.,0.};wi.U={7.,-9.,20.};
        initial=model->initializePTY(wi.p,wi.T,wi.Y);
        TetraPressureOutlet dilute(*model,9e3,300.,{1.,0.});dilute.apply(wi,*initial.state,n,wo,qo);
        const double expected=20.+24.71935841718838;
        if(std::abs(wo.U[2]-expected)>2e-5)
            throw std::runtime_error("dilute Riemann invariant mismatch");
        TetraPressureOutlet choked(*model,1000.,300.,{1.,0.});choked.apply(wi,*initial.state,n,wo,qo);
        if(qo.p<=1000.||std::abs(wo.U[2]/qo.acoustic.soundSpeed()-1.)>2e-8)
            throw std::runtime_error("pressure outlet sonic choke failed");
        wi.U[2]=2*initial.state->acoustic.soundSpeed();dilute.apply(wi,*initial.state,n,wo,qo);
        if(wo.U!=wi.U||qo.p!=initial.state->p)throw std::runtime_error("supersonic outlet responds to exterior pressure");
        wi.U[2]=-1.;TetraPressureOutlet backflow(*model,1.01e4,310.,{0.,1.});backflow.apply(wi,*initial.state,n,wo,qo);
        if(wo.U[2]>=0.||wo.Y!=std::vector<double>({0.,1.})||std::abs(qo.p-1.01e4)>1.||wo.T!=310.)
            throw std::runtime_error("reservoir backflow conditions failed");
        flux=TetraPressureOutlet::physicalFlux(wo,qo,n);
        if(flux.rho>=0.||flux.rhoY[0]!=0.||flux.rhoY[1]!=flux.rho)
            throw std::runtime_error("backflow flux uses wrong reservoir species");
        auto invalid=*initial.state;invalid.acoustic.evidence=thermo::AcousticEvidence::UNVERIFIED;
        bool rejected=false;try{backflow.apply(wi,invalid,n,wo,qo);}catch(const std::runtime_error&){rejected=true;}
        if(!rejected)throw std::runtime_error("unverified acoustic outlet state accepted");
        std::cout<<"TETRA_PRESSURE_OUTLET_PASS\n";return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
