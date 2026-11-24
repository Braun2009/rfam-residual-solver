#ifndef SectorWallHeatFlux_H
#define SectorWallHeatFlux_H

#include "Domain.H"
#include "SectorMeshMetrics.H"
#include "SpeciesRegistry.H"
#include "TransportProps.H"
#include "SectorTurbulenceOperator.H"
#include "WallFunction.H"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace Foam
{

struct SectorWallHeatFluxResult
{
    double area=0.0,distance=0.0;
    double laminarConductivity=0.0,turbulentConductivity=0.0;
    double turbulentViscosity=0.0,turbulentK=0.0,turbulentOmega=0.0;
    double frictionVelocity=0.0,yPlus=0.0,Prandtl=0.0;
    double laminarFlux=0.0,turbulentFlux=0.0,totalFlux=0.0;
};

class SectorWallHeatFlux
{
    static double magnitude(const FaceAreaVector& a)
    { return std::sqrt(a.x*a.x+a.y*a.y+a.z*a.z); }

    static void validateLayout(
        const Domain& d,const SectorMeshMetrics& mesh,int i,int k)
    {
        if(d.Nx()<4 || d.Ny()<4 || d.Nz()<4)
            throw std::runtime_error(
                "sector wall transport cell -1 status=INVALID_INPUT "
                "detail=INSUFFICIENT_HALO_DIMENSIONS");
        if(mesh.Nx()!=d.Nx() || mesh.Nr()!=d.Ny()
           ||mesh.Ntheta()!=d.Nz())
            throw std::runtime_error(
                "sector wall transport cell -1 status=INVALID_INPUT "
                "detail=DOMAIN_MESH_MISMATCH");
        if(i<0 || i>=d.Nx() || k<0 || k>=d.Nz())
            throw std::runtime_error(
                "sector wall transport cell -1 status=INVALID_INPUT "
                "detail=WALL_INDEX_OUT_OF_RANGE");
    }

public:
    static SectorWallHeatFluxResult evaluateCached(
        const Domain& d,const SectorMeshMetrics& mesh,int i,int k,
        const TransportProps::Point& gasTransport,
        const TransportProps::Point& wallTransport,
        bool turbulenceEnabled=false,double turbulentViscosity=0.0,
        double turbulentK=0.0,double turbulentOmega=0.0)
    {
        validateLayout(d,mesh,i,k);
        if(!gasTransport.valid() || !wallTransport.valid())
            throw std::runtime_error(
                "sector wall transport cell -1 status=INVALID_INPUT "
                "detail=INVALID_CACHED_TRANSPORT_POINT");
        SectorWallHeatFluxResult r;
        const int cellJ=d.Ny()-2,ghostJ=d.Ny()-1;
        const FaceAreaVector areaVector=mesh.radialAreaVector(i,ghostJ,k,+1);
        r.area=magnitude(areaVector);
        if(r.area<=1e-30) return r;
        const FaceAreaVector cell=mesh.cellCenter(i,cellJ,k);
        const FaceAreaVector ghost=mesh.cellCenter(i,ghostJ,k);
        r.distance=0.5*std::max(((ghost.x-cell.x)*areaVector.x
            +(ghost.y-cell.y)*areaVector.y+(ghost.z-cell.z)*areaVector.z)/r.area,1e-12);
        const PrimitiveVars& gas=d.W(i,cellJ,k);
        const PrimitiveVars& wall=d.W(i,ghostJ,k);
        r.laminarConductivity=0.5*(
            gasTransport.kappa+wallTransport.kappa);
        if(!std::isfinite(r.laminarConductivity)
           ||r.laminarConductivity<=0.0)
            throw std::runtime_error(
                "invalid cached sector wall conductivity");
        const double temperatureGradient=(gas.T-wall.T)/r.distance;
        r.laminarFlux=r.laminarConductivity*temperatureGradient;
        r.totalFlux=r.laminarFlux;
        if(!turbulenceEnabled) return r;

        const double rho=d.Q(i,cellJ,k).rho;
        const double cp=gasTransport.Cp;
        const double nu=gasTransport.nu;
        if(!std::isfinite(turbulentViscosity)||turbulentViscosity<0.0
           ||!std::isfinite(turbulentK)||turbulentK<0.0
           ||!std::isfinite(turbulentOmega)||turbulentOmega<=0.0)
            throw std::runtime_error(
                "invalid cached sector wall turbulence state");
        const double normalVelocity=(gas.U[0]*areaVector.x
            +gas.U[1]*areaVector.y+gas.U[2]*areaVector.z)/r.area;
        const double speedSquared=gas.U[0]*gas.U[0]+gas.U[1]*gas.U[1]
                                 +gas.U[2]*gas.U[2];
        const double tangentialSpeed=std::sqrt(std::max(
            speedSquared-normalVelocity*normalVelocity,0.0));
        r.frictionVelocity=WallFunction::frictionVelocity(
            r.distance,tangentialSpeed,nu);
        r.yPlus=r.distance*r.frictionVelocity/std::max(nu,1e-30);
        r.Prandtl=rho*nu*cp/std::max(r.laminarConductivity,1e-30);
        r.totalFlux=WallFunction::thermalHeatFlux(
            r.distance,tangentialSpeed,nu,rho,cp,r.laminarConductivity,
            gas.T,wall.T);
        r.turbulentFlux=r.totalFlux-r.laminarFlux;
        r.turbulentViscosity=turbulentViscosity;
        r.turbulentConductivity=rho*cp*turbulentViscosity/0.9;
        r.turbulentK=turbulentK;
        r.turbulentOmega=turbulentOmega;
        return r;
    }

    static SectorWallHeatFluxResult evaluate(
        const Domain& d,const SectorMeshMetrics& mesh,
        const SpeciesRegistry& registry,
        const SectorTurbulenceOperator* turbulence,int i,int k,
        bool diagnosticMode=false)
    {
        validateLayout(d,mesh,i,k);
        if(turbulence&&!turbulence->meshMatches(d))
            throw std::runtime_error(
                "sector wall transport cell -1 status=INVALID_INPUT "
                "detail=TURBULENCE_DOMAIN_MESH_MISMATCH");
        const int cellJ=d.Ny()-2,ghostJ=d.Ny()-1;
        const auto& names=d.speciesNames();
        const TransportProps::Result gas=diagnosticMode
            ?TransportProps::evaluateDiagnosticGas(
                names,d.W(i,cellJ,k),d.Q(i,cellJ,k),
                d.thermoModel(),registry).result
            :TransportProps::evaluateCommittedLegacyGas(
                names,d.W(i,cellJ,k),d.Q(i,cellJ,k),
                d.thermoModel(),registry);
        if(!gas.success())
        {
            const int cell=k*mesh.Nx()*mesh.Nr()+cellJ*mesh.Nx()+i;
            throw std::runtime_error(
                "sector wall transport cell "+std::to_string(cell)
                +" status="+TransportProps::statusName(gas.status));
        }
        const TransportProps::Result wall=diagnosticMode
            ?TransportProps::evaluateDiagnosticGas(
                names,d.W(i,ghostJ,k),d.Q(i,ghostJ,k),
                d.thermoModel(),registry).result
            :TransportProps::evaluateCommittedLegacyGas(
                names,d.W(i,ghostJ,k),d.Q(i,ghostJ,k),
                d.thermoModel(),registry);
        if(!wall.success())
        {
            const int cell=k*mesh.Nx()*mesh.Nr()+ghostJ*mesh.Nx()+i;
            throw std::runtime_error(
                "sector wall transport cell "+std::to_string(cell)
                +" status="+TransportProps::statusName(wall.status));
        }

        double nut=0.0,kValue=0.0,omegaValue=0.0;
        if(turbulence)
        {
            nut=turbulence->turbulentViscosity(
                i,cellJ,k,d,gas.point.nu);
            kValue=turbulence->turbulentKValue(i,cellJ,k);
            omegaValue=turbulence->turbulentOmegaValue(i,cellJ,k);
        }
        return evaluateCached(
            d,mesh,i,k,gas.point,wall.point,turbulence!=nullptr,
            nut,kValue,omegaValue);
    }
};

}

#endif
