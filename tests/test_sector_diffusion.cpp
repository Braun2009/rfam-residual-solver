#include <algorithm>
#include <cmath>
#include <iostream>
#include "src/framework/Domain.H"
#include "src/framework/SectorDiffusionOperator.H"
#include "src/mesh/SectorMeshMetrics.H"
#include "src/thermo/SpeciesRegistry.H"

using namespace Foam;

int main()
{
    const int nx=14,nr=10,nt=10;
    SectorMeshMetrics mesh(vortexProfilePath(),nx,nr,nt,120.0);
    SpeciesRegistry registry;
    const std::vector<std::string> names={"O2","C3H8"};
    Domain d(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2); d.setSpecies(names);
    for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
    {
        PrimitiveVars w(2); w.p=4e6; w.T=700; w.U={0.0,0.0,0.0}; w.Y={0.7,0.3};
        if(j==nr-1) w.Y={1.0,0.0};
        if(i==nx/2&&j==nr/2&&k==nt/2) { w.T=900; w.Y={0.2,0.8}; }
        d.stagedW(i,j,k)=w;
    }
    d.primitiveToConserved(); d.beginStep();
    TransportProps::Point localTransport;
    const auto& referenceCell=d.W(nx/2,nr/2,nt/2);
    if(!TransportProps::evaluateSinglePhaseGas(
        names,referenceCell.Y,referenceCell.p,referenceCell.T,
        d.U(nx/2,nr/2,nt/2).rho,registry,localTransport)) return 1;
    RealFluidState referenceThermo; referenceThermo.setSpecies(names);
    if(referenceThermo.updatePTY(
        referenceCell.p,referenceCell.T,referenceCell.Y)
       !=RealFluidState::Status::SUCCESS ||
       std::abs(localTransport.Cv-referenceThermo.Cv())>
           1e-12*referenceThermo.Cv() ||
       std::abs(localTransport.Cp-referenceThermo.Cp())>
           1e-12*referenceThermo.Cp())
    { std::cerr<<"sector transport does not use real-fluid Cp/Cv\n"; return 1; }
    SectorDiffusionOperator diffusion(mesh,registry);
    const double stable=diffusion.maxStableDt(d,0.2);
    if(!std::isfinite(stable)||stable<=0.0) return 1;
    diffusion.apply(d,std::min(stable,1e-8));

    double speciesChange[2]={0.0,0.0},energyChange=0.0,maxMassClosure=0.0;
    for(int i=1;i<nx-1;++i) for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
    {
        const double v=mesh.volume(i,j,k);
        if(!std::isfinite(d.Unew(i,j,k).rhoE)) return 1;
        double local=0.0;
        for(int s=0;s<2;++s)
        {
            if(!std::isfinite(d.Unew(i,j,k).rhoY[s])) return 1;
            const double change=d.Unew(i,j,k).rhoY[s]-d.U(i,j,k).rhoY[s];
            speciesChange[s]+=v*change; local+=change;
        }
        maxMassClosure=std::max(maxMassClosure,std::abs(local));
        energyChange+=v*(d.Unew(i,j,k).rhoE-d.U(i,j,k).rhoE);
    }
    for(double change:speciesChange) if(std::abs(change)>1e-14) return 1;
    if(std::abs(energyChange)>1e-9||maxMassClosure>1e-12) return 1;

    Domain wallDomain(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2);
    wallDomain.setSpecies(names);
    for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
    {
        PrimitiveVars w(2); w.p=4e6; w.T=j==nr-1?700.0:1000.0;
        if(j==nr-1) w.U={0.0,0.0,0.0}; else w.U={80.0,20.0,0.0};
        w.Y={0.7,0.3}; wallDomain.stagedW(i,j,k)=w;
    }
    wallDomain.primitiveToConserved(); wallDomain.beginStep();
    SectorTurbulenceOperator turbulence(mesh,registry,0.05,0.001);
    SectorDiffusionOperator ransDiffusion(mesh,registry,&turbulence);
    double wallPower=0.0;
    for(int i=1;i<nx-1;++i) for(int k=1;k<nt-1;++k)
    {
        const auto q=SectorWallHeatFlux::evaluate(
            wallDomain,mesh,registry,&turbulence,i,k);
        wallPower+=q.totalFlux*q.area;
    }
    if(!std::isfinite(wallPower)||wallPower<=0.0) return 1;
    const double wallDt=std::min(1e-9,ransDiffusion.maxStableDt(wallDomain,0.2));
    ransDiffusion.apply(wallDomain,wallDt);
    double wallEnergyChange=0.0,wallSpeciesChange[2]={0.0,0.0};
    for(int i=1;i<nx-1;++i) for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
    {
        const double volume=mesh.volume(i,j,k);
        wallEnergyChange+=volume*(wallDomain.Unew(i,j,k).rhoE
            -wallDomain.U(i,j,k).rhoE);
        for(int s=0;s<2;++s) wallSpeciesChange[s]+=volume*(
            wallDomain.Unew(i,j,k).rhoY[s]-wallDomain.U(i,j,k).rhoY[s]);
    }
    const double wallEnergyTolerance=5e-7*std::abs(wallDt*wallPower)+1e-14;
    if(std::abs(wallEnergyChange+wallDt*wallPower)>wallEnergyTolerance)
    {
        std::cerr<<"wall-function energy flux is not conservative: change="
                 <<wallEnergyChange<<" expected="<<-wallDt*wallPower
                 <<" residual="<<wallEnergyChange+wallDt*wallPower<<"\n";
        return 1;
    }
    for(double change:wallSpeciesChange) if(std::abs(change)>1e-14) return 1;

    Domain viscous(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2);
    viscous.setSpecies(names);
    for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
    {
        PrimitiveVars w(2); w.p=4e6; w.T=700.0;
        w.U={0.0,0.0,0.0}; w.Y={0.7,0.3};
        if(i==nx/2&&j==nr/2&&k==nt/2) w.U={60.0,-15.0,8.0};
        viscous.stagedW(i,j,k)=w;
    }
    viscous.primitiveToConserved(); viscous.beginStep();
    SectorDiffusionOperator viscousDiffusion(mesh,registry);
    viscousDiffusion.apply(viscous,1e-9);
    double viscousMomentum[3]={0.0,0.0,0.0};
    double viscousEnergy=0.0,maxLocalMomentum=0.0;
    for(int i=1;i<nx-1;++i) for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
    {
        const double volume=mesh.volume(i,j,k);
        for(int m=0;m<3;++m)
        {
            const double change=viscous.Unew(i,j,k).rhoU[m]
                               -viscous.U(i,j,k).rhoU[m];
            viscousMomentum[m]+=volume*change;
            maxLocalMomentum=std::max(maxLocalMomentum,std::abs(change));
        }
        viscousEnergy+=volume*(
            viscous.Unew(i,j,k).rhoE-viscous.U(i,j,k).rhoE);
    }
    if(!std::isfinite(maxLocalMomentum)||maxLocalMomentum<=0.0)
    { std::cerr<<"molecular stress did not reach momentum equation\n"; return 1; }
    for(double change:viscousMomentum)
        if(!std::isfinite(change)||std::abs(change)>1e-14)
        { std::cerr<<"molecular momentum flux is not conservative\n"; return 1; }
    if(!std::isfinite(viscousEnergy)||std::abs(viscousEnergy)>1e-12)
    { std::cerr<<"molecular stress work is not conservative\n"; return 1; }
    std::cout<<"SECTOR_DIFFUSION_PASS dt="<<stable<<"\n";
    return 0;
}
