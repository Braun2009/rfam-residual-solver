// Read-only finite-volume balance on a cell-centroid-defined core control
// volume. This is not an exact cylindrical surface or a steady-state verdict.
#include "TetraMesh.H"
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
struct Field {PrimitiveVars w{2};double rho=0.,a=0.,h=0.;};
std::array<double,3> velocity(const PrimitiveVars& w)
{return {w.U[0],w.U[1],w.U[2]};}
std::array<double,3> add(std::array<double,3> a,const std::array<double,3>& b,double s)
{for(int k=0;k<3;++k)a[k]+=s*b[k];return a;}
double magnitude(const std::array<double,3>& a)
{return std::hypot(a[0],a[1],a[2]);}
struct Budget
{
    std::size_t faces=0;
    double mass=0.,oxygen=0.,fuel=0.,energy=0.,pressureScale=0.,inboundMass=0.;
    double inboundOxygen=0.,inboundPureOxygen=0.,inboundFuelRich=0.;
    std::array<double,3> total{},convective{},pressure{};
    double inboundConvectiveNorm=0.;
    void include(const ConservedVars& flux,const std::array<double,3>& n,
                 const PrimitiveVars& upwind,double area,double sign,
                 double pressureHalf)
    {
        ++faces;
        const double m=sign*area*flux.rho;
        mass+=m;oxygen+=sign*area*flux.rhoY[0];fuel+=sign*area*flux.rhoY[1];
        energy+=sign*area*flux.rhoE;
        const auto u=velocity(upwind);
        for(int k=0;k<3;++k)
        {
            total[k]+=sign*area*flux.rhoU[k];
            convective[k]+=m*u[k];
            pressure[k]+=sign*area*pressureHalf*n[k];
        }
        pressureScale+=area*std::abs(pressureHalf);
        if(m<0.)
        {
            inboundMass-=m;
            inboundOxygen-=m*upwind.Y[0];
            if(upwind.Y[0]>.95)inboundPureOxygen-=m;
            if(upwind.Y[0]<.2)inboundFuelRich-=m;
            inboundConvectiveNorm+=-m*magnitude(u);
        }
    }
};
}
int main(int argc,char**argv)
{
    if(argc!=5){std::cerr<<"usage: audit_tetra_core_balance MESH FIELD BOARD_Z RADIUS_M\n";return 2;}
    try
    {
        TetraMesh mesh;mesh.read(argv[1]);
        const double board=std::stod(argv[3]),radius=std::stod(argv[4]);
        if(!(radius>0.))throw std::runtime_error("invalid core radius");
        std::ifstream file(argv[2],std::ios::binary);
        char magic[9]{};std::uint64_t n=0,ns=0;
        file.read(magic,9);file.read(reinterpret_cast<char*>(&n),8);
        file.read(reinterpret_cast<char*>(&ns),8);
        if(!file||std::string(magic,9)!="RFAMTETF1"||n!=mesh.cellCount()||ns!=2)
            throw std::runtime_error("field/mesh header mismatch");
        std::vector<float> centres(3*n),states(8*n),y(2*n);
        file.read(reinterpret_cast<char*>(centres.data()),centres.size()*sizeof(float));
        file.read(reinterpret_cast<char*>(states.data()),states.size()*sizeof(float));
        file.read(reinterpret_cast<char*>(y.data()),y.size()*sizeof(float));
        if(!file)throw std::runtime_error("truncated field");
        const auto state=[&](std::size_t cell)
        {
            Field f;const auto k=8*cell,j=2*cell;
            f.w.p=states[k];f.w.T=states[k+1];
            f.w.U={states[k+2],states[k+3],states[k+4]};
            const double total=double(y[j])+y[j+1];
            if(!(total>0.))throw std::runtime_error("invalid species fractions");
            f.w.Y={y[j]/total,y[j+1]/total};
            f.rho=states[k+5];f.a=states[k+6];f.h=states[k+7];return f;
        };
        const double cx=.04921534822760748,cy=.04835004867381059;
        std::vector<char> in(n,0);
        double volume=0.;std::size_t count=0;
        for(std::size_t i=0;i<n;++i)
        {
            const auto& centre=mesh.cellCentroid(i);
            const double r=std::hypot(centre.x-cx,centre.y-cy);
            if(r<radius&&centre.z>=board&&centre.z<board+.02)
            {in[i]=1;volume+=mesh.volume(i);++count;}
        }
        if(!count)throw std::runtime_error("empty core control volume");
        AUSMPlusUpFlux ausm({"O2","C3H8"});
        const auto* table=thermo::tableAcoustic();
        if(!table)throw std::runtime_error("audit needs RFAM_TETRA_TABLE");
        ausm.setConsistentDensityHook([table](double Y,double p,double h,double& rho)
            {return table->stateFromPHY(Y,p,h,rho);});
        Budget side,lower,upper,wall,all;
        for(std::size_t fi=0;fi<mesh.faceCount();++fi)
        {
            const auto& f=mesh.face(fi);
            const bool a=in[f.owner],b=f.neighbour>=0&&in[f.neighbour];
            if(a==b)continue;
            if(f.neighbour<0&&f.boundaryZone!=0&&f.boundaryZone!=30)
                throw std::runtime_error("non-wall physical boundary inside control volume");
            const auto left=state(f.owner);
            Field right;
            if(f.neighbour>=0)right=state(f.neighbour);
            else
            {
                right=left;
                const auto av=mesh.faceAreaVector(fi);
                const double area=std::hypot(av.x,av.y,av.z);
                const std::array<double,3> normal{av.x/area,av.y/area,av.z/area};
                double uNormal=0.;for(int k=0;k<3;++k)uNormal+=normal[k]*right.w.U[k];
                for(int k=0;k<3;++k)right.w.U[k]-=2.*normal[k]*uNormal;
            }
            const auto av=mesh.faceAreaVector(fi);
            const double area=std::hypot(av.x,av.y,av.z);
            const std::vector<double> normal{av.x/area,av.y/area,av.z/area};
            ConservedVars flux(2);
            ausm.computeFluxCached(left.w,right.w,normal,flux,
                left.rho,left.a,left.h,right.rho,right.a,right.h);
            const double sign=a?1.:-1.;
            // AUSM upwind velocity has the sign of the mass flux (rho>0).
            const auto& upwind=flux.rho>=0.?left.w:right.w;
            const double pHalf=(flux.rhoU[0]-flux.rho*upwind.U[0])*normal[0]
                +(flux.rhoU[1]-flux.rho*upwind.U[1])*normal[1]
                +(flux.rhoU[2]-flux.rho*upwind.U[2])*normal[2];
            Budget* category=&wall;
            if(f.neighbour>=0)
            {
                const auto& inside=mesh.cellCentroid(a?f.owner:f.neighbour);
                const auto& outside=mesh.cellCentroid(a?f.neighbour:f.owner);
                const double dr=std::hypot(outside.x-cx,outside.y-cy)
                    -std::hypot(inside.x-cx,inside.y-cy);
                if(dr>0.&&outside.z>=board&&outside.z<board+.02)category=&side;
                else if(outside.z<board)category=&lower;
                else if(outside.z>=board+.02)category=&upper;
                else category=&side;
            }
            const std::array<double,3> norm{normal[0],normal[1],normal[2]};
            category->include(flux,norm,upwind,area,sign,pHalf);
            all.include(flux,norm,upwind,area,sign,pHalf);
        }
        const auto print=[](const char* name,const Budget& s)
        {
            std::cout<<std::setprecision(12)<<name<<" faces="<<s.faces
                     <<" mass_out_kg_s="<<s.mass<<" inlet_kg_s="<<s.inboundMass
                     <<" O2_out_kg_s="<<s.oxygen<<" fuel_out_kg_s="<<s.fuel
                     <<" energy_out_W="<<s.energy
                     <<" inlet_mean_YO2="<<(s.inboundMass>0.?s.inboundOxygen/s.inboundMass:0.)
                     <<" inlet_pureO2_kg_s="<<s.inboundPureOxygen
                     <<" inlet_fuelrich_kg_s="<<s.inboundFuelRich
                     <<" force_vector_N=("<<s.total[0]<<','<<s.total[1]<<','<<s.total[2]<<')'
                     <<" convective_xy_N="<<std::hypot(s.convective[0],s.convective[1])
                     <<" inbound_convective_norm_N="<<s.inboundConvectiveNorm
                     <<" pressure_scale_N="<<s.pressureScale<<'\n';
        };
        std::cout<<"CORE_CONTROL_DIAGNOSTIC rounded_field=1 converged_unverified=1"
                 <<" cell_centroid_surface=1 cells="<<count<<" volume_m3="<<volume<<'\n';
        print("side",side);print("lower",lower);print("upper",upper);print("wall",wall);print("all",all);
        std::cout<<"table_fallbacks="<<ausm.consistentDensityFallbacks()<<'\n';
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
