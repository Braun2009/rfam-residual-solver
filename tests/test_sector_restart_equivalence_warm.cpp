#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "Domain.H"
#include "SectorCheckpoint.H"
#include "SectorMeshMetrics.H"
#include "SectorTurbulenceOperator.H"
#include "SpeciesRegistry.H"

using namespace Foam;

namespace
{

void initialize(
    Domain& domain,const SectorMeshMetrics& mesh,
    const SpeciesRegistry& registry)
{
    domain.setSpecies(registry.names());
    std::vector<double> composition(registry.size(),0.0);
    composition[registry.index("O2")]=0.8;
    composition[registry.index("C3H8")]=0.2;
    for(int i=0;i<domain.Nx();++i)
    for(int j=0;j<domain.Ny();++j)
    for(int k=0;k<domain.Nz();++k)
    {
        PrimitiveVars& state=domain.stagedW(i,j,k);
        const auto center=mesh.cellCenter(i,j,k);
        state.p=4.0e6;
        state.T=700.0;
        state.U={30.0+2.0*center.x,-3.0+center.y,0.5*center.z};
        state.Y=composition;
    }
    domain.primitiveToConserved();
    if(domain.tripped())
        throw std::runtime_error("warm restart initial state failed");
}

void advance(
    Domain& domain,SectorTurbulenceOperator& turbulence,double dt)
{
    domain.beginStep();
    turbulence.apply(domain,dt);
    domain.commit();
    if(domain.tripped())
        throw std::runtime_error("warm restart advancement tripped");
}

bool close(double lhs,double rhs)
{
    return std::abs(lhs-rhs)
        <=1.0e-13*std::max({std::abs(lhs),std::abs(rhs),1.0});
}

} // namespace

int main()
{
    const std::string path="test_sector_restart_equivalence_warm.chk";
    SpeciesRegistry registry;
    const int nx=8,ny=6,nz=6;
    SectorMeshMetrics mesh(vortexProfilePath(),nx,ny,nz,120.0);

    Domain continuous(
        nx,ny,nz,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
    initialize(continuous,mesh,registry);
    SectorTurbulenceOperator continuousTurbulence(mesh,registry);

    Domain split(
        nx,ny,nz,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
    initialize(split,mesh,registry);
    SectorTurbulenceOperator splitTurbulence(mesh,registry);

    const double dt=1.0e-10;
    advance(continuous,continuousTurbulence,dt);
    advance(continuous,continuousTurbulence,dt);

    advance(split,splitTurbulence,dt);
    SectorCheckpointMetadata metadata{
        dt,1,1,"warm-gas-rans","restart-test",12345};
    SectorCheckpoint::write(path,split,metadata,&splitTurbulence);

    Domain restarted(
        nx,ny,nz,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
    restarted.setSpecies(registry.names());
    SectorTurbulenceOperator restartedTurbulence(mesh,registry);
    const auto loaded=SectorCheckpoint::read(
        path,restarted,&restartedTurbulence);
    std::remove(path.c_str());
    if(loaded.iteration!=1 || !close(loaded.time,dt)) return 1;
    advance(restarted,restartedTurbulence,dt);

    for(int i=0;i<nx;++i) for(int j=0;j<ny;++j) for(int k=0;k<nz;++k)
    {
        const auto& lhs=continuous.U(i,j,k);
        const auto& rhs=restarted.U(i,j,k);
        if(!close(lhs.rho,rhs.rho) || !close(lhs.rhoE,rhs.rhoE) ||
           !close(lhs.rhoZ,rhs.rhoZ) || !close(lhs.rhoC,rhs.rhoC))
            return 1;
        for(int m=0;m<3;++m)
            if(!close(lhs.rhoU[m],rhs.rhoU[m])) return 1;
        for(size_t s=0;s<lhs.rhoY.size();++s)
            if(!close(lhs.rhoY[s],rhs.rhoY[s])) return 1;
        if(!close(continuous.W(i,j,k).T,restarted.W(i,j,k).T)) return 1;
    }

    const auto& continuousK=continuousTurbulence.kField();
    const auto& restartedK=restartedTurbulence.kField();
    const auto& continuousOmega=continuousTurbulence.omegaField();
    const auto& restartedOmega=restartedTurbulence.omegaField();
    if(continuousK.size()!=restartedK.size() ||
       continuousOmega.size()!=restartedOmega.size()) return 1;
    for(size_t index=0;index<continuousK.size();++index)
        if(!close(continuousK[index],restartedK[index]) ||
           !close(continuousOmega[index],restartedOmega[index])) return 1;

    std::cout<<"SECTOR_RESTART_EQUIVALENCE_WARM_PASS\n";
    return 0;
}
