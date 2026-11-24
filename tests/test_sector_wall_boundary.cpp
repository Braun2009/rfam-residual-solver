#include <algorithm>
#include <cmath>
#include <iostream>

#include "src/framework/Domain.H"
#include "src/framework/SectorDiagnostics.H"
#include "src/framework/SectorWallBoundary.H"
#include "src/mesh/SectorMeshMetrics.H"

using namespace Foam;

namespace
{

bool close(double lhs,double rhs,double relativeTolerance=1.0e-8)
{
    return std::isfinite(lhs)&&std::isfinite(rhs)
        &&std::abs(lhs-rhs)<=relativeTolerance
            *std::max({std::abs(lhs),std::abs(rhs),1.0});
}

}

int main()
{
    const int nx=14,nr=10,nt=10;
    SectorMeshMetrics mesh(vortexProfilePath(),nx,nr,nt,120.0);
    Domain d(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2);
    for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
    {
        PrimitiveVars& w=d.stagedW(i,j,k);
        w.p=3.9e6+1000.0*i; w.T=900.0+2.0*j;
        w.U={80.0+i,5.0-j,3.0+k}; w.Y={0.7,0.3};
    }
    d.primitiveToConserved();
    SectorWallBoundary wall(700.0); wall.apply(d);

    for(int i=0;i<nx;++i) for(int k=0;k<nt;++k)
    {
        const PrimitiveVars& inside=d.W(i,nr-2,k);
        const PrimitiveVars& ghost=d.W(i,nr-1,k);
        if(!close(ghost.p,inside.p)||!close(ghost.T,700.0)
           ||ghost.Y.size()!=inside.Y.size()) return 1;
        for(std::size_t species=0;species<ghost.Y.size();++species)
            if(!close(ghost.Y[species],inside.Y[species])) return 1;
        for(int m=0;m<3;++m)
            if(!close(ghost.U[m],-inside.U[m])) return 1;
    }

    const SectorBalance balance=SectorDiagnostics(mesh).evaluate(d);
    if(std::abs(balance.wallOutwardMassFlow)>1e-14)
    {
        std::cerr << "wall mass flux is not zero: "
                  << balance.wallOutwardMassFlow << "\n";
        return 1;
    }
    std::cout << "SECTOR_WALL_BOUNDARY_PASS\n";
    return 0;
}
