#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include "src/framework/Domain.H"
#include "src/framework/SectorOutputWriter.H"
#include "src/mesh/SectorMeshMetrics.H"
#include "src/thermo/SpeciesRegistry.H"

using namespace Foam;

int main()
{
    const std::string field="test_sector_field.csv",wall="test_sector_wall.csv";
    SpeciesRegistry registry;SectorMeshMetrics mesh(vortexProfilePath(),8,6,6,120.0);
    Domain d(8,6,6,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
    d.setSpecies(registry.names());std::vector<double> y(registry.size(),0.0);
    y[registry.index("CO2")]=0.6468;y[registry.index("H2O")]=0.3532;
    for(int i=0;i<8;++i)for(int j=0;j<6;++j)for(int k=0;k<6;++k)
    {auto& w=d.stagedW(i,j,k);w.p=4e6;w.T=j==5?700:1000;w.U={0,0,0};w.Y=y;}
    d.primitiveToConserved();SectorOutputWriter::writeFieldCsv(field,d,mesh);
    SectorOutputWriter::writeWallCsv(wall,d,mesh,registry);
    std::ifstream f(field),q(wall);std::string fh,qh,fd,qd;
    std::getline(f,fh);std::getline(f,fd);std::getline(q,qh);std::getline(q,qd);
    std::remove(field.c_str());std::remove(wall.c_str());
    if(fh.find("T_K")==std::string::npos||fh.find("Y_CO2")==std::string::npos||fd.empty())return 1;
    if(qh.find("qwall_W_m2")==std::string::npos
       ||qh.find("area_m2")==std::string::npos
       ||qh.find("wall_power_W")==std::string::npos
       ||qh.find("q_lam_W_m2")==std::string::npos
       ||qh.find("q_turb_W_m2")==std::string::npos
       ||qh.find("turb_k_m2_s2")==std::string::npos
       ||qh.find("turb_omega_s-1")==std::string::npos||qd.empty())return 1;
    std::cout<<"SECTOR_OUTPUT_PASS\n";return 0;
}
