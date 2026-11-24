// Tetrahedral finite-volume kernel checks on the audited nozzle mesh.
//
// Case 1 (freestream): a uniform state against a matching ghost must leave
// every cell residual at roundoff, because per-cell outward area vectors sum
// to zero and internal faces cancel exactly between owner and neighbour.
// This validates face orientation, neighbour bookkeeping and residual scaling
// before any real boundary physics is attached.
#include "TetraMesh.H"
#include "TetraFiniteVolume.H"
#include "EquilibriumThermoModel.H"
#include "FluxScheme.H"
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

using namespace Foam;

int main(int argc,char**argv)
{
    if(argc!=2){std::cerr<<"usage: test_tetra_freestream <mesh.rfamtet>\n";return 2;}
    try
    {
        TetraMesh mesh;mesh.read(argv[1]);
        std::vector<std::string> species{"O2","C3H8"};
        auto model=thermo::makeLegacySinglePhaseThermoModel(species);
        PrimitiveVars w(species.size());
        w.p=4.0e6;w.T=700.0;w.U={0.0,0.0,0.0};w.Y={0.8,0.2};
        const auto initialized=model->initializePTY(w.p,w.T,w.Y);
        if(!initialized.success())
        {std::cerr<<"uniform state rejected: "<<thermo::equilibriumStatusName(initialized.status)<<"\n";return 1;}
        const thermo::EquilibriumState q=*initialized.state;
        if(!solverConsumerAdmitsAcousticState(q))
        {std::cerr<<"uniform state acoustic closure not admitted\n";return 1;}

        TetraFiniteVolume fv(mesh,species,
            std::unique_ptr<FluxScheme>(new AUSMFluxAdapter(species)));
        for(std::size_t c=0;c<mesh.cellCount();++c)fv.setPrimitive(c,w,q);
        // A quiescent uniform state must be its own ghost state.
        const TetraBoundaryProvider boundary=[&](std::size_t,std::size_t,
                                                 PrimitiveVars& gw,thermo::EquilibriumState& gq)
        {gw=w;gq=q;};
        std::vector<ConservedVars> residual;
        fv.computeResidual(boundary,residual);
        if(residual.size()!=mesh.cellCount()){std::cerr<<"residual size mismatch\n";return 1;}

        // Scale from the state itself, not from an absolute tolerance.
        // A stationary uniform state has zero mass flux, so the only scale
        // that matters for momentum is the pressure force p*|A| per cell.
        const double rhoScale=q.rho;
        const double energyScale=q.rho*q.h;
        std::vector<double> pressureForce(mesh.cellCount(),0.0);
        for(std::size_t c=0;c<mesh.cellCount();++c)
            for(std::int32_t fi:mesh.cellFaces(c))
            {
                const auto a=mesh.faceAreaVector(fi);
                pressureForce[c]+=w.p*std::sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
            }
        double worstMass=0.0,worstMomentum=0.0,worstEnergy=0.0,worstSpecies=0.0;
        std::size_t worstMomentumCell=0;
        for(std::size_t c=0;c<residual.size();++c)
        {
            const auto& r=residual[c];
            worstMass=std::max(worstMass,std::abs(r.rho));
            for(double m:r.rhoU)
                if(std::abs(m)>worstMomentum){worstMomentum=std::abs(m);worstMomentumCell=c;}
            worstEnergy=std::max(worstEnergy,std::abs(r.rhoE));
            for(double y:r.rhoY)worstSpecies=std::max(worstSpecies,std::abs(y));
        }
        const double massRel=worstMass/mesh.volume(0)/rhoScale;
        const double momentumRel=worstMomentum/std::max(pressureForce[worstMomentumCell],1e-30);
        std::cout<<"worst residual per volume: mass="<<worstMass/mesh.volume(0)
                 <<" momentum="<<worstMomentum/mesh.volume(0)
                 <<" energy="<<worstEnergy/mesh.volume(0)
                 <<" species="<<worstSpecies/mesh.volume(0)<<"\n";
        std::cout<<"relative mass residual="<<massRel
                 <<" momentum vs pressure force="<<momentumRel
                 <<" (cell "<<worstMomentumCell<<")\n";
        if(massRel>1e-12)
        {std::cerr<<"freestream mass residual too large\n";return 1;}
        if(momentumRel>1e-12)
        {std::cerr<<"freestream momentum residual too large\n";return 1;}
        if(std::abs(worstEnergy/mesh.volume(0))/std::max(energyScale,1e-30)>1e-12)
        {std::cerr<<"freestream energy residual too large\n";return 1;}
        if(std::abs(worstSpecies)/mesh.volume(0)/std::max(rhoScale,1e-30)>1e-12)
        {std::cerr<<"freestream species residual too large\n";return 1;}

        std::cout<<"TETRA_FREESTREAM_PASS cells="<<mesh.cellCount()
                 <<" faces="<<mesh.faceCount()<<"\n";
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
