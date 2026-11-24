#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../src/framework/SectorConservativeGridTransfer.H"
#include "../src/mesh/SectorMeshMetrics.H"
#include "../src/thermo/SpeciesRegistry.H"

using namespace Foam;

namespace
{

void addScaled(
    ConservedVars& target,const ConservedVars& source,double scale)
{
    target.rho+=scale*source.rho;
    target.rhoE+=scale*source.rhoE;
    for(std::size_t component=0;component<target.rhoU.size();++component)
        target.rhoU[component]+=scale*source.rhoU[component];
    for(std::size_t component=0;component<target.rhoY.size();++component)
        target.rhoY[component]+=scale*source.rhoY[component];
    target.rhoZ+=scale*source.rhoZ;
    target.rhoC+=scale*source.rhoC;
}

ConservedVars inventory(
    const Domain& domain,const SectorMeshMetrics& mesh)
{
    ConservedVars total(domain.nSpecies());
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
        addScaled(total,domain.U(i,j,k),mesh.volume(i,j,k));
    return total;
}

double maximumRelativeDifference(
    const ConservedVars& left,const ConservedVars& right)
{
    double result=0.0;
    auto compare=[&](double a,double b)
    {
        result=std::max(result,
            std::abs(a-b)/std::max({std::abs(a),std::abs(b),1.0e-300}));
    };
    compare(left.rho,right.rho);
    compare(left.rhoE,right.rhoE);
    for(std::size_t component=0;component<left.rhoU.size();++component)
        compare(left.rhoU[component],right.rhoU[component]);
    for(std::size_t component=0;component<left.rhoY.size();++component)
        compare(left.rhoY[component],right.rhoY[component]);
    compare(left.rhoZ,right.rhoZ);
    compare(left.rhoC,right.rhoC);
    return result;
}

bool initializeNonuniform(Domain& domain,const SpeciesRegistry& registry)
{
    const int oxygen=registry.index("O2");
    const int nitrogen=registry.index("N2");
    const int carbonDioxide=registry.index("CO2");
    const int water=registry.index("H2O");
    for(int i=0;i<domain.Nx();++i)
    for(int j=0;j<domain.Ny();++j)
    for(int k=0;k<domain.Nz();++k)
    {
        PrimitiveVars& primitive=domain.stagedW(i,j,k);
        const double x=static_cast<double>(i)/std::max(domain.Nx()-1,1);
        const double r=static_cast<double>(j)/std::max(domain.Ny()-1,1);
        const double theta=static_cast<double>(k)/std::max(domain.Nz()-1,1);
        primitive.p=8.0e5+2.0e5*x+5.0e4*r;
        primitive.T=500.0+300.0*x+40.0*theta;
        primitive.U={20.0+15.0*x,-2.0+3.0*r,1.0-2.0*theta};
        primitive.Y.assign(static_cast<std::size_t>(registry.size()),0.0);
        primitive.Y[static_cast<std::size_t>(oxygen)]=0.15+0.05*x;
        primitive.Y[static_cast<std::size_t>(nitrogen)]=0.25+0.05*r;
        primitive.Y[static_cast<std::size_t>(carbonDioxide)]=0.35-0.05*x;
        primitive.Y[static_cast<std::size_t>(water)]=
            1.0-primitive.Y[static_cast<std::size_t>(oxygen)]
               -primitive.Y[static_cast<std::size_t>(nitrogen)]
               -primitive.Y[static_cast<std::size_t>(carbonDioxide)];
        primitive.Z=0.1+0.2*x;
        primitive.C=0.05+0.1*r;
    }
    return domain.primitiveToConserved();
}

} // namespace

int main()
{
    SpeciesRegistry registry;
    const std::string profile=vortexProfilePath();
    // Active fine=8x6x4, coarse=4x3x2: exact integer 2x2x2 topology.
    SectorMeshMetrics fineMesh(profile,10,8,6,120.0);
    SectorMeshMetrics coarseMesh(profile,6,5,4,120.0);
    Domain fine(10,8,6,fineMesh.length(),0.065,
        fineMesh.sectorAngle(),registry.size());
    Domain coarse(6,5,4,coarseMesh.length(),0.065,
        coarseMesh.sectorAngle(),registry.size());
    Domain roundTrip(10,8,6,fineMesh.length(),0.065,
        fineMesh.sectorAngle(),registry.size());
    fine.setSpecies(registry.names());
    coarse.setSpecies(registry.names());
    roundTrip.setSpecies(registry.names());
    if(!initializeNonuniform(fine,registry)) return 1;

    const auto restrictionPlan=
        SectorConservativeGridTransfer::plan(fine,coarse);
    if(restrictionPlan.direction
        !=SectorConservativeGridTransfer::Direction::RESTRICT
        ||restrictionPlan.ratioI!=2
        ||restrictionPlan.ratioJ!=2
        ||restrictionPlan.ratioK!=2) return 1;
    SectorConservativeGridTransfer::apply(
        fine,fineMesh,nullptr,coarse,coarseMesh,nullptr);
    if(coarse.tripped()||!coarse.stateReady()) return 1;

    const ConservedVars fineInventory=inventory(fine,fineMesh);
    const ConservedVars coarseInventory=inventory(coarse,coarseMesh);
    const double restrictionError=maximumRelativeDifference(
        fineInventory,coarseInventory);
    if(restrictionError>2.0e-12) return 1;

    const auto prolongationPlan=
        SectorConservativeGridTransfer::plan(coarse,roundTrip);
    if(prolongationPlan.direction
        !=SectorConservativeGridTransfer::Direction::PROLONG
        ||prolongationPlan.ratioI!=2
        ||prolongationPlan.ratioJ!=2
        ||prolongationPlan.ratioK!=2) return 1;
    SectorConservativeGridTransfer::apply(
        coarse,coarseMesh,nullptr,roundTrip,fineMesh,nullptr,4);
    if(roundTrip.tripped()||!roundTrip.stateReady()) return 1;
    const ConservedVars roundTripInventory=inventory(roundTrip,fineMesh);
    const double prolongationError=maximumRelativeDifference(
        coarseInventory,roundTripInventory);
    const double roundTripError=maximumRelativeDifference(
        fineInventory,roundTripInventory);
    if(prolongationError>2.0e-12||roundTripError>3.0e-12) return 1;

    std::cout << "restriction/prolongation/roundtrip="
              << restrictionError << "/" << prolongationError
              << "/" << roundTripError << "\n"
              << "SECTOR_CONSERVATIVE_GRID_TRANSFER_PASS\n";
    return 0;
}
