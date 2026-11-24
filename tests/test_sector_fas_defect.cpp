#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../src/framework/FixedConservedSourceOperator.H"
#include "../src/framework/SectorConservativeGridTransfer.H"
#include "../src/framework/SectorFASDefect.H"
#include "../src/mesh/SectorMeshMetrics.H"
#include "../src/thermo/SpeciesRegistry.H"

using namespace Foam;

namespace
{

int index(const Domain& domain,int i,int j,int k)
{ return k*domain.Nx()*domain.Ny()+j*domain.Nx()+i; }

bool initialize(Domain& domain,const SpeciesRegistry& registry)
{
    const int oxygen=registry.index("O2");
    const int nitrogen=registry.index("N2");
    for(int i=0;i<domain.Nx();++i)
    for(int j=0;j<domain.Ny();++j)
    for(int k=0;k<domain.Nz();++k)
    {
        PrimitiveVars& primitive=domain.stagedW(i,j,k);
        primitive.p=1.0e6;
        primitive.T=600.0+2.0*i+j;
        primitive.U={40.0,0.0,0.0};
        primitive.Y.assign(static_cast<std::size_t>(registry.size()),0.0);
        primitive.Y[static_cast<std::size_t>(oxygen)]=0.3;
        primitive.Y[static_cast<std::size_t>(nitrogen)]=0.7;
        primitive.Z=0.2;
        primitive.C=0.1;
    }
    return domain.primitiveToConserved();
}

class MeshRateOperator : public Operator
{
    double coefficient_;

public:
    explicit MeshRateOperator(double coefficient)
        :coefficient_(coefficient)
    {}

    void apply(Domain& domain,double dt) const override
    {
        for(int i=0;i<domain.Nx();++i)
        for(int j=0;j<domain.Ny();++j)
        for(int k=0;k<domain.Nz();++k)
        {
            const double spatial=1.0+0.1*i+0.03*j+0.007*k;
            ConservedVars& candidate=domain.Unew(i,j,k);
            candidate.rhoZ+=dt*coefficient_*spatial;
            candidate.rhoC+=dt*(0.25*coefficient_+0.02*i-0.01*j);
        }
    }

    std::string name() const override { return "mesh-dependent rate"; }
};

bool sameField(
    const std::vector<ConservedVars>& left,
    const std::vector<ConservedVars>& right,
    const Domain& domain,double tolerance)
{
    if(left.size()!=right.size()) return false;
    auto compare=[&](double a,double b)
    {
        return std::abs(a-b)
            <=tolerance*std::max({std::abs(a),std::abs(b),1.0});
    };
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const std::size_t cell=static_cast<std::size_t>(
            index(domain,i,j,k));
        const ConservedVars& a=left[cell];
        const ConservedVars& b=right[cell];
        if(!compare(a.rho,b.rho)||!compare(a.rhoE,b.rhoE)
            ||!compare(a.rhoZ,b.rhoZ)||!compare(a.rhoC,b.rhoC))
            return false;
        for(std::size_t component=0;component<a.rhoU.size();++component)
            if(!compare(a.rhoU[component],b.rhoU[component])) return false;
        for(std::size_t component=0;component<a.rhoY.size();++component)
            if(!compare(a.rhoY[component],b.rhoY[component])) return false;
    }
    return true;
}

} // namespace

int main()
{
    SpeciesRegistry registry;
    const std::string profile=vortexProfilePath();
    SectorMeshMetrics fineMesh(profile,10,8,6,120.0);
    SectorMeshMetrics coarseMesh(profile,6,5,4,120.0);
    Domain fine(10,8,6,fineMesh.length(),0.065,
        fineMesh.sectorAngle(),registry.size());
    Domain coarse(6,5,4,coarseMesh.length(),0.065,
        coarseMesh.sectorAngle(),registry.size());
    fine.setSpecies(registry.names());
    coarse.setSpecies(registry.names());
    if(!initialize(fine,registry)) return 1;
    SectorConservativeGridTransfer::apply(
        fine,fineMesh,nullptr,coarse,coarseMesh,nullptr);

    std::vector<std::unique_ptr<Operator>> fineOperators;
    fineOperators.push_back(std::make_unique<MeshRateOperator>(1.0));
    std::vector<std::unique_ptr<Operator>> coarseOperators;
    coarseOperators.push_back(std::make_unique<MeshRateOperator>(4.0));
    SteadyResidualEvaluator::Result fineResidual,coarseResidual;
    if(!SteadyResidualEvaluator::evaluate(
            fine,fineOperators,1.0e-6,fineResidual,{},true)
        ||!SteadyResidualEvaluator::evaluate(
            coarse,coarseOperators,1.0e-6,coarseResidual,{},true))
        return 1;

    const SectorFASDefect::Result defect=SectorFASDefect::build(
        fine,fineMesh,fineResidual,coarse,coarseMesh,coarseResidual);
    if(defect.maximumInitialConsistencyError>2.0e-12
        ||defect.maximumRestrictionInventoryError>2.0e-12)
    {
        std::cerr << "FAS construction error consistency/inventory="
                  << defect.maximumInitialConsistencyError << "/"
                  << defect.maximumRestrictionInventoryError << "\n";
        return 1;
    }

    std::vector<std::unique_ptr<Operator>> forcedCoarseOperators;
    forcedCoarseOperators.push_back(std::make_unique<MeshRateOperator>(4.0));
    forcedCoarseOperators.push_back(
        std::make_unique<FixedConservedSourceOperator>(
            defect.coarseSource));
    SteadyResidualEvaluator::Result forcedResidual;
    if(!SteadyResidualEvaluator::evaluate(
            coarse,forcedCoarseOperators,1.0e-6,
            forcedResidual,{},true))
        return 1;
    if(!sameField(
            forcedResidual.field,defect.restrictedFineResidual,
            coarse,3.0e-10))
    {
        std::cerr << "forced coarse residual does not equal restricted fine residual\n";
        return 1;
    }

    bool nonzeroTau=false;
    for(int i=1;i<coarse.Nx()-1;++i)
    for(int j=1;j<coarse.Ny()-1;++j)
    for(int k=1;k<coarse.Nz()-1;++k)
    {
        const ConservedVars& tau=defect.tau[static_cast<std::size_t>(
            index(coarse,i,j,k))];
        if(std::abs(tau.rhoZ)>1.0e-12||std::abs(tau.rhoC)>1.0e-12)
            nonzeroTau=true;
    }
    if(!nonzeroTau)
    {
        std::cerr << "FAS tau unexpectedly vanished for unequal fine/coarse operators\n";
        return 1;
    }

    std::cout << "consistency_error="
              << defect.maximumInitialConsistencyError
              << " restriction_inventory_error="
              << defect.maximumRestrictionInventoryError << "\n"
              << "SECTOR_FAS_DEFECT_PASS\n";
    return 0;
}
