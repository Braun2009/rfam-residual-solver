#include <cmath>
#include <iostream>

#include "src/framework/Domain.H"
#include "src/mesh/SectorMeshMetrics.H"
#include "src/turbulence/SectorTurbulenceOperator.H"
#include "src/turbulence/WallFunction.H"

using namespace Foam;

int main()
{
    {
        const double y=1e-8,rho=5.0,cp=2000.0,conductivity=0.1;
        const double q=WallFunction::thermalHeatFlux(
            y,1e-4,1.5e-5,rho,cp,conductivity,800.0,700.0);
        const double fourier=conductivity*(800.0-700.0)/y;
        if(!std::isfinite(q)||std::abs(q-fourier)>1e-8*fourier)
        { std::cerr<<"thermal wall law misses viscous limit\n"; return 1; }
        const double stagnant=WallFunction::thermalHeatFlux(
            1e-3,0.0,1.5e-5,rho,cp,conductivity,800.0,700.0);
        if(std::abs(stagnant-conductivity*(800.0-700.0)/1e-3)>1e-12*stagnant)
        { std::cerr<<"thermal wall law misses stagnant limit\n"; return 1; }

        const double logLayer=WallFunction::thermalHeatFlux(
            1e-3,100.0,1.5e-5,6.5,2000.0,0.13,2200.0,700.0);
        const double reverse=WallFunction::thermalHeatFlux(
            1e-3,100.0,1.5e-5,6.5,2000.0,0.13,700.0,2200.0);
        if(!std::isfinite(logLayer)||logLayer<=0.0||logLayer>=1e7
           ||std::abs(logLayer+reverse)>1e-12*logLayer)
        { std::cerr<<"invalid log-layer thermal wall flux\n"; return 1; }
    }
    const int nx=16,nr=10,nt=10;
    SectorMeshMetrics mesh(vortexProfilePath(),nx,nr,nt,120.0);
    SpeciesRegistry registry;
    const std::vector<std::string> species={"O2","C3H8"};
    Domain d(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2);
    d.setSpecies(species);
    for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
    {
        PrimitiveVars w(2); w.p=4e6; w.T=700;
        w.U={20.0+80.0*j/(nr-1.0),0.0,0.0}; w.Y[0]=0.8; w.Y[1]=0.2;
        d.stagedW(i,j,k)=w;
    }
    d.primitiveToConserved(); d.beginStep();
    SectorTurbulenceOperator turbulence(mesh,registry,0.05,0.001);
    const double wallNuForBoundary=turbulence.molecularKinematicViscosity(
        nx/2,nr-2,nt/2,d);
    const double dt=std::min(1e-8,turbulence.maxStableDt(d,0.2));
    turbulence.apply(d,dt); d.commit();
    if(d.tripped()) return 1;
    double maxNut=0;
    for(int i=1;i<nx-1;++i) for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
        maxNut=std::max(maxNut,turbulence.turbulentViscosity(i,j,k,d));
    const double centerNu=turbulence.molecularKinematicViscosity(
        nx/2,nr/2,nt/2,d);
    if(!std::isfinite(centerNu)||centerNu<=0.0||
       std::abs(centerNu-1.5e-5)<1e-8*1.5e-5)
    { std::cerr<<"sector SST still uses the fixed legacy nu\n"; return 1; }
    if(!std::isfinite(maxNut)||maxNut<=centerNu) return 1;
    const auto& kField=turbulence.kField();const auto& omegaField=turbulence.omegaField();
    auto fieldIndex=[&](int i,int j,int k){return k*nx*nr+j*nx+i;};
    for(size_t id=0;id<kField.size();++id)
        if(!std::isfinite(kField[id])||!std::isfinite(omegaField[id])
           ||kField[id]<KOmegaSST::k_min||omegaField[id]<KOmegaSST::om_min) return 1;
    for(int i=0;i<nx;++i) for(int j=0;j<nr;++j)
    {
        if(kField[fieldIndex(i,j,0)]!=kField[fieldIndex(i,j,nt-2)]) return 1;
        if(kField[fieldIndex(i,j,nt-1)]!=kField[fieldIndex(i,j,1)]) return 1;
    }
    for(int i=0;i<nx;++i) for(int k=0;k<nt;++k)
    {
        if(kField[fieldIndex(i,0,k)]!=kField[fieldIndex(i,1,k)]) return 1;
        if(kField[fieldIndex(i,nr-1,k)]!=KOmegaSST::k_min) return 1;
    }
    {
        const int i=nx/2,k=nt/2;
        const auto area=mesh.radialAreaVector(i,nr-1,k,+1);
        const auto cell=mesh.cellCenter(i,nr-2,k),ghost=mesh.cellCenter(i,nr-1,k);
        const double areaMagnitude=std::sqrt(
            area.x*area.x+area.y*area.y+area.z*area.z);
        const double wallDistance=0.5*std::abs(
            ((ghost.x-cell.x)*area.x+(ghost.y-cell.y)*area.y
             +(ghost.z-cell.z)*area.z)/areaMagnitude);
        const double expectedOmega=std::max(
            60.0*wallNuForBoundary/(KOmegaSST::beta1*wallDistance*wallDistance),
            KOmegaSST::om_min);
        const double actualOmega=omegaField[fieldIndex(i,nr-1,k)];
        if(std::abs(actualOmega-expectedOmega)>1e-12*expectedOmega)
        {
            std::cerr<<"wall omega uses incorrect halo-adjusted distance: "
                     <<actualOmega<<" vs "<<expectedOmega<<"\n";
            return 1;
        }
    }

    if(TransportProps::supportsTwoPhase||TransportProps::supportsDenseCryogenic)
    { std::cerr<<"gas transport closure overclaims phase capability\n"; return 1; }
    {
        SectorTurbulenceOperator gradientTurbulence(mesh,registry,0.05,0.001);
        std::vector<double> gradientK(d.N()),gradientOmega(d.N());
        for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
        {
            const auto c=mesh.cellCenter(i,j,k);
            const size_t id=static_cast<size_t>(k*nx*nr+j*nx+i);
            gradientK[id]=2.0+20.0*c.x+5.0*c.y+3.0*c.z;
            gradientOmega[id]=1000.0+200.0*c.x+50.0*c.y+30.0*c.z;
        }
        gradientTurbulence.restoreFields(gradientK,gradientOmega);
        const double cd=gradientTurbulence.crossDiffusionSource(
            nx/2,nr/2,nt/2);
        if(!std::isfinite(cd)||cd<=0.0)
        { std::cerr<<"mapped SST cross diffusion was not evaluated\n"; return 1; }
    }

    Domain wallFlow(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2);
    wallFlow.setSpecies(species);
    for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
    {
        PrimitiveVars w(2); w.p=4e6; w.T=700;
        w.U=j==nr-1?std::vector<double>{0.0,0.0,0.0}
                    :std::vector<double>{50.0,0.0,0.0};
        w.Y[0]=0.8; w.Y[1]=0.2; wallFlow.stagedW(i,j,k)=w;
    }
    wallFlow.primitiveToConserved(); wallFlow.beginStep();
    SectorTurbulenceOperator wallTurbulence(mesh,registry,0.05,0.001);
    wallTurbulence.apply(wallFlow,1e-10);
    double wallMomentumChange=0.0;
    for(int i=1;i<nx-1;++i) for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
        wallMomentumChange+=mesh.volume(i,j,k)
            *(wallFlow.Unew(i,j,k).rhoU[0]-wallFlow.U(i,j,k).rhoU[0]);
    if(!(wallMomentumChange<0.0))
    { std::cerr<<"wall function did not oppose tangential momentum\n"; return 1; }

    // A uniform Cartesian velocity has no shear on the curved mapped mesh.
    Domain uniform(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2);
    uniform.setSpecies(species);
    for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
    {
        PrimitiveVars w(2); w.p=4e6; w.T=700; w.U={30.0,-7.0,4.0};
        w.Y[0]=0.8; w.Y[1]=0.2; uniform.stagedW(i,j,k)=w;
    }
    uniform.primitiveToConserved(); uniform.beginStep();
    SectorTurbulenceOperator uniformTurbulence(mesh,registry,0.05,0.001);
    uniformTurbulence.apply(uniform,1e-10);
    double maxMomentumChange=0.0,uniformIntegratedChange[3]={0.0,0.0,0.0};
    double uniformIntegratedEnergyChange=0.0;
    for(int i=1;i<nx-1;++i) for(int j=1;j<nr-2;++j) for(int k=1;k<nt-1;++k)
        for(int m=0;m<3;++m)
            maxMomentumChange=std::max(maxMomentumChange,
                std::abs(uniform.Unew(i,j,k).rhoU[m]-uniform.U(i,j,k).rhoU[m]));
    if(!std::isfinite(maxMomentumChange)||maxMomentumChange>1e-8) return 1;
    for(int i=1;i<nx-1;++i) for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
    {
        for(int m=0;m<3;++m)
            uniformIntegratedChange[m]+=mesh.volume(i,j,k)
                *(uniform.Unew(i,j,k).rhoU[m]-uniform.U(i,j,k).rhoU[m]);
        uniformIntegratedEnergyChange+=mesh.volume(i,j,k)
            *(uniform.Unew(i,j,k).rhoE-uniform.U(i,j,k).rhoE);
    }

    // A compact velocity perturbation may redistribute momentum, but mapped
    // internal stress faces must not create net Cartesian momentum.
    Domain conservative(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2);
    conservative.setSpecies(species);
    for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
    {
        PrimitiveVars w(2); w.p=4e6; w.T=700; w.U={30.0,-7.0,4.0};
        if(i==nx/2 && j==nr/2 && k==nt/2) w.U={80.0,12.0,-9.0};
        w.Y[0]=0.8; w.Y[1]=0.2; conservative.stagedW(i,j,k)=w;
    }
    conservative.primitiveToConserved(); conservative.beginStep();
    SectorTurbulenceOperator conservativeTurbulence(mesh,registry,0.05,0.001);
    conservativeTurbulence.apply(conservative,1e-10);
    double integratedMomentumChange[3]={0.0,0.0,0.0};
    double integratedEnergyChange=0.0;
    for(int i=1;i<nx-1;++i) for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
    {
        for(int m=0;m<3;++m)
            integratedMomentumChange[m]+=mesh.volume(i,j,k)
                *(conservative.Unew(i,j,k).rhoU[m]-conservative.U(i,j,k).rhoU[m]);
        integratedEnergyChange+=mesh.volume(i,j,k)
            *(conservative.Unew(i,j,k).rhoE-conservative.U(i,j,k).rhoE);
    }
    for(int m=0;m<3;++m)
    {
        const double change=integratedMomentumChange[m]-uniformIntegratedChange[m];
        if(!std::isfinite(change)||std::abs(change)>1e-12) return 1;
    }
    if(!std::isfinite(integratedEnergyChange)||
       std::abs(integratedEnergyChange-uniformIntegratedEnergyChange)>1e-9)
    { std::cerr<<"turbulent stress work is not face conservative\n"; return 1; }
    std::cout<<"SECTOR_TURBULENCE_PASS nut="<<maxNut<<"\n";
    return 0;
}
