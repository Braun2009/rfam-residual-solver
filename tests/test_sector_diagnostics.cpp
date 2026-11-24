#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include "src/framework/Domain.H"
#include "src/framework/SectorConvectionOperator.H"
#include "src/framework/SectorDiagnostics.H"
#include "src/mesh/SectorMeshMetrics.H"

using namespace Foam;

int main()
{
    const int nx=12,nr=9,nt=10;
    SectorMeshMetrics mesh(vortexProfilePath(),nx,nr,nt,120.0);
    Domain d(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2);
    for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
    {
        const double theta=mesh.thetaCenter(k);
        PrimitiveVars w(2); w.p=4e6; w.T=700;
        w.U={25.0,-5.0*std::sin(theta),5.0*std::cos(theta)}; w.Y={0.8,0.2};
        d.stagedW(i,j,k)=w;
    }
    d.primitiveToConserved(); d.applyCyclicZ(120.0);
    SectorDiagnostics diagnostics(mesh);
    const SectorBalance b=diagnostics.evaluate(d);

    double expectedMass=0.0;
    for(int i=1;i<nx-1;++i) for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
        expectedMass+=d.U(i,j,k).rho*mesh.volume(i,j,k);
    if(std::abs(b.massInventory-expectedMass)>1e-13*expectedMass) return 1;
    if(std::abs(b.periodicNetMassFlow)>1e-14||b.periodicRelativeMismatch>1e-14) return 1;

    double expectedInlet=0.0,expectedOutlet=0.0;
    for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
    {
        expectedInlet+=d.U(1,j,k).rho*25.0*mesh.axialArea(1,j,k);
        expectedOutlet+=d.U(nx-2,j,k).rho*25.0*mesh.axialArea(nx-1,j,k);
    }
    if(std::abs(b.inletMassFlow-expectedInlet)>1e-13*expectedInlet) return 1;
    if(std::abs(b.outletMassFlow-expectedOutlet)>1e-13*expectedOutlet) return 1;
    if(std::abs(b.axisOutwardMassFlow)>1e-14) return 1;

    const double dt=1e-10;
    SectorConvectionOperator convection(
        std::make_unique<AUSMFluxAdapter>(),mesh);
    d.beginStep(); convection.apply(d,dt); d.commit();
    if(d.tripped()) return 1;
    const SectorBalance after=diagnostics.evaluate(d);
    const double expectedChange=dt*(b.inletMassFlow-b.outletMassFlow
        -b.axisOutwardMassFlow-b.wallOutwardMassFlow-b.periodicNetMassFlow);
    const double actualChange=after.massInventory-b.massInventory;
    const double tolerance=std::max(1e-8*std::abs(expectedChange),
                                    1e-13*b.massInventory);
    if(std::abs(actualChange-expectedChange)>tolerance)
    {
        std::cerr << "sector mass closure failed: actual/expected="
                  << actualChange << "/" << expectedChange << "\n";
        return 1;
    }
    std::cout<<"SECTOR_DIAGNOSTICS_PASS mass="<<b.massInventory<<"\n";
    return 0;
}
