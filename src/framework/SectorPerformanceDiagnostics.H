#ifndef SectorPerformanceDiagnostics_H
#define SectorPerformanceDiagnostics_H

#include "Domain.H"
#include "SectorMeshMetrics.H"
#include "SectorWallHeatFlux.H"
#include "SpeciesRegistry.H"
#include "TransportProps.H"
#include "SectorTurbulenceOperator.H"
#include <algorithm>
#include <cmath>

namespace Foam
{

struct SectorPerformance
{
    double chamberPressure,throatMassFlow,throatArea,cStar;
    double wallHeatFluxMaximum,wallHeatFluxAverage,wallHeatPower;
};

class SectorPerformanceDiagnostics
{
    const SectorMeshMetrics& mesh_;
    const SpeciesRegistry& registry_;

    static double magnitude(const FaceAreaVector& a)
    { return std::sqrt(a.x*a.x+a.y*a.y+a.z*a.z); }

public:
    SectorPerformanceDiagnostics(const SectorMeshMetrics& mesh,
                                 const SpeciesRegistry& registry)
        :mesh_(mesh),registry_(registry) {}

    SectorPerformance evaluate(const Domain& d,
        const SectorTurbulenceOperator* turbulence=nullptr) const
    {
        SectorPerformance r{0.0,0.0,0.0,0.0,0.0,0.0,0.0};
        double chamberVolume=0.0;
        for(int i=1;i<d.Nx()-1;++i)
        {
            if(0.5*(mesh_.xFace(i)+mesh_.xFace(i+1))>0.273) continue;
            for(int j=1;j<d.Ny()-1;++j) for(int k=1;k<d.Nz()-1;++k)
            {
                const double volume=mesh_.volume(i,j,k);
                r.chamberPressure+=d.W(i,j,k).p*volume; chamberVolume+=volume;
            }
        }
        r.chamberPressure/=std::max(chamberVolume,1e-30);

        int throatFace=1;
        for(int face=1;face<d.Nx();++face)
            if(std::abs(mesh_.xFace(face)-0.2765)<std::abs(mesh_.xFace(throatFace)-0.2765))
                throatFace=face;
        const int left=std::clamp(throatFace-1,1,d.Nx()-2);
        const int right=std::clamp(throatFace,1,d.Nx()-2);
        for(int j=1;j<d.Ny()-1;++j) for(int k=1;k<d.Nz()-1;++k)
        {
            const double area=mesh_.axialArea(throatFace,j,k);
            r.throatArea+=area;
            r.throatMassFlow+=0.5*(d.U(left,j,k).rho*d.W(left,j,k).U[0]
                +d.U(right,j,k).rho*d.W(right,j,k).U[0])*area;
        }
        if(r.throatMassFlow>1e-6)
            r.cStar=r.chamberPressure*r.throatArea/r.throatMassFlow;

        double wallArea=0.0;
        for(int i=1;i<d.Nx()-1;++i) for(int k=1;k<d.Nz()-1;++k)
        {
            const SectorWallHeatFluxResult q=SectorWallHeatFlux::evaluate(
                d,mesh_,registry_,turbulence,i,k,
                turbulence?turbulence->diagnosticMode():false);
            if(q.area<=1e-30) continue;
            r.wallHeatFluxMaximum=std::max(r.wallHeatFluxMaximum,q.totalFlux);
            r.wallHeatPower+=q.totalFlux*q.area; wallArea+=q.area;
        }
        r.wallHeatFluxAverage=r.wallHeatPower/std::max(wallArea,1e-30);
        return r;
    }
};

} // End namespace Foam
#endif
