// Diagnostic surface flux through the injector/chamber centroid partition.
// This is NOT an injector exit plane: crossing faces are selected by the
// positions of their neighbouring cell centroids. Field values are float32.
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
struct Cell
{
    PrimitiveVars w{2};
    double rho=0.,a=0.,h=0.;
};
}
int main(int argc,char** argv)
{
    if(argc!=5){std::cerr<<"usage: audit_tetra_crossing_flux MESH FIELD BOARD_Z OUTPUT_CSV\n";return 2;}
    try
    {
        TetraMesh mesh;mesh.read(argv[1]);
        const double board=std::stod(argv[3]);
        std::ifstream input(argv[2],std::ios::binary);
        char magic[9]{};std::uint64_t n=0,ns=0;
        input.read(magic,9);input.read(reinterpret_cast<char*>(&n),8);
        input.read(reinterpret_cast<char*>(&ns),8);
        if(!input||std::string(magic,9)!="RFAMTETF1"||n!=mesh.cellCount()||ns!=2)
            throw std::runtime_error("field header/mesh mismatch");
        std::vector<float> centres(3*n),field(8*n),ys(2*n);
        input.read(reinterpret_cast<char*>(centres.data()),centres.size()*sizeof(float));
        input.read(reinterpret_cast<char*>(field.data()),field.size()*sizeof(float));
        input.read(reinterpret_cast<char*>(ys.data()),ys.size()*sizeof(float));
        if(!input)throw std::runtime_error("truncated field");
        AUSMPlusUpFlux ausm({"O2","C3H8"});
        const auto* table=thermo::tableAcoustic();
        if(!table)throw std::runtime_error("audit requires RFAM_TETRA_TABLE");
        ausm.setConsistentDensityHook([table](double y,double p,double h,double& rho)
            {return table->stateFromPHY(y,p,h,rho);});
        std::ofstream out(argv[4]);
        if(!out)throw std::runtime_error("cannot open CSV");
        out<<std::setprecision(17)
           <<"face,area_m2,r_mm,z_mm,mass_to_chamber_kg_s,O2_to_chamber_kg_s,"
             "fuel_to_chamber_kg_s,energy_to_chamber_W\n";
        const auto state=[&](std::size_t c)
        {
            Cell s;const std::size_t i=8*c,j=2*c;
            s.w.p=field[i];s.w.T=field[i+1];
            s.w.U={field[i+2],field[i+3],field[i+4]};
            const double sum=double(ys[j])+double(ys[j+1]);
            if(!(sum>0.))throw std::runtime_error("invalid mass fractions");
            s.w.Y={ys[j]/sum,ys[j+1]/sum};
            s.rho=field[i+5];s.a=field[i+6];s.h=field[i+7];
            return s;
        };
        double areaSum=0.,mass=0.,oxygen=0.,fuel=0.,energy=0.;
        std::size_t count=0;
        for(std::size_t f=0;f<mesh.faceCount();++f)
        {
            const auto& face=mesh.face(f);
            if(face.neighbour<0)continue;
            const std::size_t P=face.owner,N=face.neighbour;
            const bool belowP=mesh.cellCentroid(P).z<board;
            if(belowP==(mesh.cellCentroid(N).z<board))continue;
            const auto lp=state(P),rn=state(N);
            const auto av=mesh.faceAreaVector(f);
            const double area=std::sqrt(av.x*av.x+av.y*av.y+av.z*av.z);
            const std::vector<double> normal{av.x/area,av.y/area,av.z/area};
            ConservedVars flux(2);
            ausm.computeFluxCached(lp.w,rn.w,normal,flux,
                lp.rho,lp.a,lp.h,rn.rho,rn.a,rn.h);
            const double sign=belowP?1.:-1.;
            const auto fc=mesh.faceCentroid(f);
            const double r=std::hypot(fc.x-0.04921534822760748,fc.y-0.04835004867381059);
            const double m=sign*area*flux.rho,o=sign*area*flux.rhoY[0];
            const double b=sign*area*flux.rhoY[1],e=sign*area*flux.rhoE;
            out<<f<<','<<area<<','<<r*1e3<<','<<(fc.z-board)*1e3<<','
               <<m<<','<<o<<','<<b<<','<<e<<'\n';
            ++count;areaSum+=area;mass+=m;oxygen+=o;fuel+=b;energy+=e;
        }
        std::cout<<std::setprecision(12)<<"CROSSING_FLUX_DIAGNOSTIC faces="<<count
                 <<" area_mm2="<<areaSum*1e6<<" mass_kg_s="<<mass
                 <<" O2_kg_s="<<oxygen<<" fuel_kg_s="<<fuel
                 <<" energy_W="<<energy<<" table_fallbacks="
                 <<ausm.consistentDensityFallbacks()
                 <<" diagnostic_only=1 solver_status_unchecked=1"<<std::endl;
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
