#ifndef SectorOutputWriter_H
#define SectorOutputWriter_H

#include "Domain.H"
#include "SectorMeshMetrics.H"
#include "SectorWallHeatFlux.H"
#include "SectorTurbulenceOperator.H"
#include "SpeciesRegistry.H"
#include "TransportProps.H"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace Foam
{

class SectorOutputWriter
{
public:
    static void writeFieldCsv(const std::string& path,const Domain& d,
                              const SectorMeshMetrics& mesh)
    {
        std::ofstream out(path);
        if(!out) throw std::runtime_error("cannot write sector field: "+path);
        out<<"i,j,k,x_m,y_m,z_m,rho_kg_m3,p_Pa,T_K,Ux_m_s,Uy_m_s,Uz_m_s";
        for(const auto& species:d.speciesNames()) out<<",Y_"<<species;
        out<<"\n"<<std::setprecision(12);
        for(int i=1;i<d.Nx()-1;++i) for(int j=1;j<d.Ny()-1;++j)
        for(int k=1;k<d.Nz()-1;++k)
        {
            const auto c=mesh.cellCenter(i,j,k);const auto& w=d.W(i,j,k);
            out<<i<<','<<j<<','<<k<<','<<c.x<<','<<c.y<<','<<c.z<<','
               <<d.U(i,j,k).rho<<','<<w.p<<','<<w.T<<','<<w.U[0]<<','<<w.U[1]
               <<','<<w.U[2]; for(double value:w.Y) out<<','<<value; out<<'\n';
        }
    }

    static void writeWallCsv(const std::string& path,const Domain& d,
                             const SectorMeshMetrics& mesh,
                             const SpeciesRegistry& registry,
                             const SectorTurbulenceOperator* turbulence=nullptr)
    {
        std::ofstream out(path);
        if(!out) throw std::runtime_error("cannot write wall field: "+path);
        out<<"i,k,x_m,theta_rad,radius_m,Tgas_K,Twall_K,qwall_W_m2,"
              "area_m2,wall_power_W,"
              "distance_m,k_lam_W_mK,nu_t_m2_s,k_turb_W_mK,"
              "q_lam_W_m2,q_turb_W_m2,turb_k_m2_s2,turb_omega_s-1,"
              "u_tau_m_s,y_plus,Pr\n"
           <<std::setprecision(12);
        const int cellJ=d.Ny()-2;
        for(int i=1;i<d.Nx()-1;++i) for(int k=1;k<d.Nz()-1;++k)
        {
            const SectorWallHeatFluxResult q=SectorWallHeatFlux::evaluate(
                d,mesh,registry,turbulence,i,k,
                turbulence?turbulence->diagnosticMode():false);
            if(q.area<=1e-30) continue;
            const auto cell=mesh.cellCenter(i,cellJ,k);
            const auto& gas=d.W(i,cellJ,k);const auto& wall=d.W(i,d.Ny()-1,k);
            out<<i<<','<<k<<','<<cell.x<<','<<mesh.thetaCenter(k)<<','
               <<mesh.radiusCenter(i,cellJ)<<','<<gas.T<<','<<wall.T<<','
               <<q.totalFlux<<','<<q.area<<','<<q.totalFlux*q.area<<','
               <<q.distance<<','<<q.laminarConductivity<<','
               <<q.turbulentViscosity<<','<<q.turbulentConductivity<<','
               <<q.laminarFlux<<','<<q.turbulentFlux<<','<<q.turbulentK<<','
               <<q.turbulentOmega<<','<<q.frictionVelocity<<','<<q.yPlus<<','
               <<q.Prandtl<<'\n';
        }
    }
};

} // End namespace Foam
#endif
