#include <cmath>
#include <iostream>
#include "src/framework/Domain.H"
#include "src/framework/SectorPerformanceDiagnostics.H"
#include "src/mesh/SectorMeshMetrics.H"
#include "src/thermo/SpeciesRegistry.H"

using namespace Foam;

int main()
{
    SpeciesRegistry registry; const int nx=20,nr=12,nt=10;
    SectorMeshMetrics mesh(vortexProfilePath(),nx,nr,nt,120.0);
    Domain d(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
    d.setSpecies(registry.names());
    std::vector<double> y(registry.size(),0.0);
    y[registry.index("CO2")]=0.6468; y[registry.index("H2O")]=0.3532;
    for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
    {
        PrimitiveVars& w=d.stagedW(i,j,k); w.p=4e6; w.T=(j==nr-1?700.0:1000.0);
        w.U={50.0,0.0,0.0}; w.Y=y;
    }
    d.primitiveToConserved();
    const SectorPerformance p=SectorPerformanceDiagnostics(mesh,registry).evaluate(d);
    if(!std::isfinite(p.cStar)||p.cStar<=0.0||p.throatMassFlow<=0.0) return 1;
    if(p.wallHeatFluxMaximum<=0.0||p.wallHeatFluxAverage<=0.0||p.wallHeatPower<=0.0) return 1;
    if(std::abs(p.chamberPressure-4e6)>1e-6) return 1;
    std::cout<<"SECTOR_PERFORMANCE_PASS cstar="<<p.cStar
             <<" qmax="<<p.wallHeatFluxMaximum<<"\n"; return 0;
}
