#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>

#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/SectorConvectionOperator.H"

using namespace Foam;

int main()
{
    const int nx=24,nr=10,nt=12;
    SectorMeshMetrics mesh(vortexProfilePath(),nx,nr,nt,120.0);
    Domain d(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2);
    for(int i=0;i<nx;++i)
    for(int j=0;j<nr;++j)
    for(int k=0;k<nt;++k)
    {
        PrimitiveVars w(2);
        w.p=4e6; w.T=700; w.U={0,0,0}; w.Y[0]=0.8; w.Y[1]=0.2;
        d.stagedW(i,j,k)=w;
    }
    d.primitiveToConserved();
    d.beginStep();
    SectorConvectionOperator op(std::unique_ptr<FluxScheme>(new AUSMFluxAdapter()),mesh);
    op.apply(d,1e-9);
    d.commit();
    if(d.tripped())
    {
        std::cerr<<"freestream domain tripped during commit\n";
        return 1;
    }

    double maxVelocity=0,maxDp=0;
    for(int i=1;i<nx-1;++i)
    for(int j=1;j<nr-1;++j)
    for(int k=1;k<nt-1;++k)
    {
        const auto& w=d.W(i,j,k);
        maxVelocity=std::max(maxVelocity,std::sqrt(w.U[0]*w.U[0]+w.U[1]*w.U[1]+w.U[2]*w.U[2]));
        maxDp=std::max(maxDp,std::abs(w.p-4e6));
    }
    const double pressureTolerance = 1.0e-9*4e6;
    if(maxVelocity>1e-7 || maxDp>pressureTolerance)
    {
        std::cerr<<"freestream drift velocity="<<maxVelocity<<" dp="<<maxDp<<"\n";
        return 1;
    }
    std::cout<<"SECTOR_FREESTREAM_PASS\n";
    return 0;
}
