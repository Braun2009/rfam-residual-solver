#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/SectorConvectionOperator.H"
#include "src/thermo/SpeciesRegistry.H"

using namespace Foam;

int main()
{
    SpeciesRegistry registry;
    const int nx=18,nr=10,nt=10;
    SectorMeshMetrics mesh(vortexProfilePath(),nx,nr,nt,120.0);
    Domain d(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
    d.setSpecies(registry.names());
    std::vector<double> y(registry.size(),0.0);
    y[registry.index("CO2")]=0.6468; y[registry.index("H2O")]=0.3532;
    for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
    {
        PrimitiveVars& w=d.stagedW(i,j,k); w.p=4e6; w.T=700; w.U={0.0,0.0,0.0}; w.Y=y;
    }
    d.primitiveToConserved(); d.beginStep();
    SectorConvectionOperator convection(
        std::make_unique<AUSMFluxAdapter>(registry.names()),mesh);
    convection.apply(d,1e-9); d.commit();
    if(d.tripped()) return 1;
    double maxVelocity=0.0,maxPressureError=0.0;
    for(int i=1;i<nx-1;++i) for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
    {
        const auto& w=d.W(i,j,k);
        maxVelocity=std::max(maxVelocity,std::sqrt(w.U[0]*w.U[0]+w.U[1]*w.U[1]+w.U[2]*w.U[2]));
        maxPressureError=std::max(maxPressureError,std::abs(w.p-4e6));
    }
    if(maxVelocity>1e-7||maxPressureError>1e-3) return 1;
    std::cout<<"SECTOR_REACTING_FREESTREAM_PASS\n"; return 0;
}
