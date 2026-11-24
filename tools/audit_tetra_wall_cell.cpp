// Read-only audit of a tetra cell's four face fluxes from a rounded field dump.
// No boundary model is inferred for inlet/outlet faces; those cells are rejected.
#include "TetraMesh.H"
#include "CoolPropThermoModel.H"
#include "AUSMPlusUpFlux.H"
#include "TableAcoustic.H"
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Foam;
namespace
{
struct FieldCell { PrimitiveVars w{2}; float rho=0, sound=0, h=0; };
FieldCell readCell(std::ifstream& file, std::uint64_t n, std::size_t c)
{
    std::array<float,8> f{};
    std::array<float,2> y{};
    file.seekg(25+static_cast<std::streamoff>(12*n+32*c));
    file.read(reinterpret_cast<char*>(f.data()),sizeof(f));
    file.seekg(25+static_cast<std::streamoff>(44*n+8*c));
    file.read(reinterpret_cast<char*>(y.data()),sizeof(y));
    if(!file)throw std::runtime_error("truncated field dump");
    FieldCell s;
    s.w.p=f[0];s.w.T=f[1];s.w.U={f[2],f[3],f[4]};
    const double sum=double(y[0])+double(y[1]);
    if(!(sum>0.))throw std::runtime_error("invalid field composition");
    s.w.Y={y[0]/sum,y[1]/sum};
    s.rho=f[5];s.sound=f[6];s.h=f[7];
    return s;
}
}
int main(int argc,char** argv)
{
    if(argc!=4){std::cerr<<"usage: audit_tetra_wall_cell MESH FIELD CELL\n";return 2;}
    try
    {
        TetraMesh mesh;mesh.read(argv[1]);
        const std::size_t cell=std::stoull(argv[3]);
        if(cell>=mesh.cellCount())throw std::runtime_error("cell out of range");
        std::ifstream file(argv[2],std::ios::binary);
        char magic[9]{};std::uint64_t count=0,ns=0;
        file.read(magic,9);file.read(reinterpret_cast<char*>(&count),8);
        file.read(reinterpret_cast<char*>(&ns),8);
        if(!file||std::string(magic,9)!="RFAMTETF1"||count!=mesh.cellCount()||ns!=2)
            throw std::runtime_error("field header/mesh/species mismatch: magic="
                +std::string(magic,9)+" cells="+std::to_string(count)
                +" mesh_cells="+std::to_string(mesh.cellCount())+" species="+std::to_string(ns));
        auto model=thermo::makeCoolPropThermoModel();
        AUSMPlusUpFlux flux({"O2","C3H8"});
        if(const auto* table=thermo::tableAcoustic())
            flux.setConsistentDensityHook([table](double y,double p,double h,double& rho)
                {return table->stateFromPHY(y,p,h,rho);});
        std::array<double,3> netMomentum{};
        double netMass=0.,netEnergy=0.,pressureForce=0.;
        std::cout<<std::setprecision(15);
        std::cout<<"field="<<argv[2]<<" cell="<<cell<<" volume_m3="<<mesh.volume(cell)
                 <<" rounded_field=1\n";
        std::cout<<"face,zone,neighbour,area_m2,normal_x,normal_y,normal_z,"
                    "u_inside_normal_m_s,p_inside_Pa,p_neighbour_Pa,"
                    "wall_pressure_Pa,mass_flux_kg_s,energy_flux_W,"
                    "momentum_x_N,momentum_y_N,momentum_z_N\n";
        for(auto f:mesh.cellFaces(cell))
        {
            const auto& face=mesh.face(f);
            const bool owner=face.owner==static_cast<std::int32_t>(cell);
            const int other=owner?face.neighbour:face.owner;
            if(other<0&&face.boundaryZone!=0)
                throw std::runtime_error("inlet/outlet needs the actual prescribed boundary flux");
            const auto& leftField=readCell(file,count,static_cast<std::size_t>(face.owner));
            const auto& left=leftField.w;
            auto lq=model->initializePTY(left.p,left.T,left.Y);
            if(!lq.success())throw std::runtime_error("left flash: "+lq.diagnostics.backendStatus);
            auto areaVector=mesh.faceAreaVector(f);
            const double area=std::sqrt(areaVector.x*areaVector.x+areaVector.y*areaVector.y+areaVector.z*areaVector.z);
            std::vector<double> n{areaVector.x/area,areaVector.y/area,areaVector.z/area};
            PrimitiveVars right(2);thermo::EquilibriumState rq;
            if(face.neighbour<0)
            {
                right=left;rq=*lq.state;
                double un=0.;for(int k=0;k<3;++k)un+=right.U[k]*n[k];
                for(int k=0;k<3;++k)right.U[k]-=2*un*n[k];
            }
            else
            {
                const auto otherField=readCell(file,count,static_cast<std::size_t>(face.neighbour));
                right=otherField.w;
                auto q=model->initializePTY(right.p,right.T,right.Y);
                if(!q.success())throw std::runtime_error("right flash: "+q.diagnostics.backendStatus);
                rq=*q.state;
            }
            ConservedVars result(2);
            flux.computeFluxCached(left,right,n,result,
                lq.state->rho,lq.state->acoustic.soundSpeed(),lq.state->h,
                rq.rho,rq.acoustic.soundSpeed(),rq.h);
            const double sign=owner?1.:-1.;
            const auto& inside=owner?left:right;
            const double un=sign*(inside.U[0]*n[0]+inside.U[1]*n[1]+inside.U[2]*n[2]);
            const double pWall=other<0?sign*(result.rhoU[0]*n[0]+result.rhoU[1]*n[1]+result.rhoU[2]*n[2]):0.;
            pressureForce+=area*inside.p;
            netMass+=sign*area*result.rho;netEnergy+=sign*area*result.rhoE;
            for(int k=0;k<3;++k)netMomentum[k]+=sign*area*result.rhoU[k];
            std::cout<<f<<','<<face.boundaryZone<<','<<other<<','<<area;
            for(auto x:n)std::cout<<','<<sign*x;
            std::cout<<','<<un<<','<<inside.p<<','<<(other<0?0.:(owner?right.p:left.p))
                     <<','<<pWall<<','<<sign*area*result.rho<<','<<sign*area*result.rhoE;
            for(int k=0;k<3;++k)std::cout<<','<<sign*area*result.rhoU[k];
            std::cout<<'\n';
        }
        const double nm=std::sqrt(netMomentum[0]*netMomentum[0]
            +netMomentum[1]*netMomentum[1]+netMomentum[2]*netMomentum[2]);
        std::cout<<"net_mass_kg_s="<<netMass<<" net_energy_W="<<netEnergy
                 <<" net_momentum_N="<<nm<<" net_momentum_over_pA="<<nm/pressureForce
                 <<" net_momentum_per_volume_N_m3="<<nm/mesh.volume(cell)
                 <<" table_fallbacks="<<flux.consistentDensityFallbacks()<<'\n';
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
