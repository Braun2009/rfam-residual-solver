#ifndef SectorDiagnostics_H
#define SectorDiagnostics_H

#include "Domain.H"
#include "FluxScheme.H"
#include "SectorMeshMetrics.H"
#include "SpeciesRegistry.H"
#include <algorithm>
#include <cmath>
#include <vector>

namespace Foam
{

struct SectorBalance
{
    double massInventory;
    double inletMassFlow;
    double outletMassFlow;
    double axisOutwardMassFlow;
    double wallOutwardMassFlow;
    double periodicNetMassFlow;
    double periodicRelativeMismatch;
};

struct SectorSpeciesBalance
{
    std::vector<double> minimumMassFraction;
    double maximumMassFractionSumError;
    double elementInventory[4];
};

class SectorDiagnostics
{
    const SectorMeshMetrics& mesh_;

    static double magnitude(const FaceAreaVector& area)
    {
        return std::sqrt(area.x*area.x+area.y*area.y+area.z*area.z);
    }

    static double numericalMassFlux(const FluxScheme& flux,
                                    const PrimitiveVars& inside,
                                    const thermo::EquilibriumState& insideThermo,
                                    const PrimitiveVars& outside,
                                    const thermo::EquilibriumState& outsideThermo,
                                    const FaceAreaVector& area)
    {
        const double scale=magnitude(area);
        if(scale<=1e-30) return 0.0;
        ConservedVars face(static_cast<int>(inside.Y.size()));
        flux.computeFluxFromEquilibrium(
            inside,insideThermo,outside,outsideThermo,
            {area.x/scale,area.y/scale,area.z/scale},face);
        return scale*face.rho;
    }

public:
    explicit SectorDiagnostics(const SectorMeshMetrics& mesh):mesh_(mesh) {}

    SectorBalance evaluate(const Domain& d) const
    {
        SectorBalance result{0.0,0.0,0.0,0.0,0.0,0.0,0.0};
        AUSMFluxAdapter flux(d.speciesNames());
        for(int i=1;i<d.Nx()-1;++i)
        for(int j=1;j<d.Ny()-1;++j)
        for(int k=1;k<d.Nz()-1;++k)
            result.massInventory+=d.U(i,j,k).rho*mesh_.volume(i,j,k);

        for(int j=1;j<d.Ny()-1;++j)
        for(int k=1;k<d.Nz()-1;++k)
        {
            result.inletMassFlow-=numericalMassFlux(
                flux,d.W(1,j,k),d.Q(1,j,k),d.W(0,j,k),d.Q(0,j,k),
                mesh_.axialAreaVector(1,j,k,-1));
            result.outletMassFlow+=numericalMassFlux(
                flux,d.W(d.Nx()-2,j,k),d.Q(d.Nx()-2,j,k),
                d.W(d.Nx()-1,j,k),d.Q(d.Nx()-1,j,k),
                mesh_.axialAreaVector(d.Nx()-1,j,k,+1));
        }

        for(int i=1;i<d.Nx()-1;++i)
        for(int k=1;k<d.Nz()-1;++k)
        {
            result.axisOutwardMassFlow+=numericalMassFlux(
                flux,d.W(i,1,k),d.Q(i,1,k),d.W(i,0,k),d.Q(i,0,k),
                mesh_.radialAreaVector(i,1,k,-1));
            result.wallOutwardMassFlow+=numericalMassFlux(
                flux,d.W(i,d.Ny()-2,k),d.Q(i,d.Ny()-2,k),
                d.W(i,d.Ny()-1,k),d.Q(i,d.Ny()-1,k),
                mesh_.radialAreaVector(i,d.Ny()-1,k,+1));
        }

        double periodicAbsolute=0.0;
        for(int i=1;i<d.Nx()-1;++i)
        for(int j=1;j<d.Ny()-1;++j)
        {
            const double lower=numericalMassFlux(
                flux,d.W(i,j,1),d.Q(i,j,1),d.W(i,j,0),d.Q(i,j,0),
                mesh_.azimuthalAreaVector(i,j,1,-1));
            const double upper=numericalMassFlux(
                flux,d.W(i,j,d.Nz()-2),d.Q(i,j,d.Nz()-2),
                d.W(i,j,d.Nz()-1),d.Q(i,j,d.Nz()-1),
                mesh_.azimuthalAreaVector(i,j,d.Nz()-1,+1));
            result.periodicNetMassFlow+=lower+upper;
            periodicAbsolute+=std::abs(lower)+std::abs(upper);
        }
        result.periodicRelativeMismatch=std::abs(result.periodicNetMassFlow)
            /std::max(periodicAbsolute,1e-30);
        return result;
    }

    SectorSpeciesBalance evaluateSpecies(const Domain& d,
                                         const SpeciesRegistry& registry) const
    {
        SectorSpeciesBalance result;
        result.minimumMassFraction.assign(d.nSpecies(),1e30);
        result.maximumMassFractionSumError=0.0;
        for(double& value:result.elementInventory) value=0.0;
        for(int i=1;i<d.Nx()-1;++i) for(int j=1;j<d.Ny()-1;++j)
        for(int k=1;k<d.Nz()-1;++k)
        {
            const PrimitiveVars& w=d.W(i,j,k);
            double sum=0.0;
            for(int s=0;s<d.nSpecies();++s)
            {
                result.minimumMassFraction[s]=std::min(result.minimumMassFraction[s],w.Y[s]);
                sum+=w.Y[s];
                const int g=registry.index(d.speciesNames()[s]);
                const double moles=d.U(i,j,k).rho*w.Y[s]*mesh_.volume(i,j,k)
                    /registry.molecularWeight(g);
                result.elementInventory[0]+=moles*registry.nC(g);
                result.elementInventory[1]+=moles*registry.nH(g);
                result.elementInventory[2]+=moles*registry.nO(g);
                result.elementInventory[3]+=moles*registry.nN(g);
            }
            result.maximumMassFractionSumError=std::max(
                result.maximumMassFractionSumError,std::abs(sum-1.0));
        }
        return result;
    }
};

} // End namespace Foam
#endif
