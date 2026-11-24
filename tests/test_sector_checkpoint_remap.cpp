#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "src/framework/Domain.H"
#include "src/framework/SectorCheckpointRemap.H"
#include "src/mesh/SectorMeshMetrics.H"
#include "src/thermo/SpeciesRegistry.H"
#include "src/turbulence/SectorTurbulenceOperator.H"

using namespace Foam;

int main()
{
    SpeciesRegistry registry;
    const std::string profile=vortexProfilePath();
    SectorMeshMetrics coarseMesh(profile,10,8,6,120.0);
    SectorMeshMetrics fineMesh(profile,18,12,10,120.0);
    Domain coarse(10,8,6,coarseMesh.length(),0.065,coarseMesh.sectorAngle(),registry.size());
    Domain fine(18,12,10,fineMesh.length(),0.065,fineMesh.sectorAngle(),registry.size());
    coarse.setSpecies(registry.names());fine.setSpecies(registry.names());
    std::vector<double> composition(registry.size(),0.0);
    composition[registry.index("CO2")]=0.65;
    composition[registry.index("H2O")]=0.35;
    for(int i=0;i<coarse.Nx();++i) for(int j=0;j<coarse.Ny();++j)
    for(int k=0;k<coarse.Nz();++k)
    {
        PrimitiveVars& w=coarse.stagedW(i,j,k);
        w.p=4.0e6;w.T=1200.0;w.U={40.0,-2.0,3.0};w.Y=composition;
    }
    coarse.primitiveToConserved();
    SectorTurbulenceOperator coarseTurbulence(coarseMesh,registry);
    SectorTurbulenceOperator fineTurbulence(fineMesh,registry);
    std::vector<double> kField(coarse.N()),omegaField(coarse.N());
    for(size_t n=0;n<kField.size();++n)
    {kField[n]=0.1+1e-5*n;omegaField[n]=100.0+1e-3*n;}
    coarseTurbulence.restoreFields(kField,omegaField);

    SectorCheckpointRemap::apply(coarse,&coarseTurbulence,fine,&fineTurbulence);
    for(int i=0;i<fine.Nx();++i) for(int j=0;j<fine.Ny();++j)
    for(int k=0;k<fine.Nz();++k)
    {
        const PrimitiveVars& w=fine.W(i,j,k);
        if(std::abs(w.p-4.0e6)>1e-8*4.0e6||std::abs(w.T-1200.0)>1e-8*1200.0
           ||std::abs(w.U[0]-40.0)>1e-10||std::abs(w.U[1]+2.0)>1e-10
           ||std::abs(w.U[2]-3.0)>1e-10) return 1;
        for(size_t s=0;s<composition.size();++s)
            if(std::abs(w.Y[s]-composition[s])>1e-10) return 1;
    }
    if(fineTurbulence.kField().size()!=static_cast<size_t>(fine.N())
       ||fineTurbulence.omegaField().size()!=static_cast<size_t>(fine.N())
       ||*std::min_element(fineTurbulence.kField().begin(),fineTurbulence.kField().end())<=0.0
       ||*std::min_element(fineTurbulence.omegaField().begin(),fineTurbulence.omegaField().end())<=0.0)
        return 1;
    std::cout << "SECTOR_CHECKPOINT_REMAP_PASS\n";
    return 0;
}
