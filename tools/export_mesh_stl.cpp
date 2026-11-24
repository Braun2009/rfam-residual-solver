// Structured chamber boundary preview only. STL has no volume-cell topology.
// CAD nozzle fluid reconstruction: tools/rebuild_injector_fluid.py
// Tetrahedral cells and connectivity: tools/mesh_injector_tetra.py (.vtu/.npz)
#include "SectorMeshMetrics.H"
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using Point=std::array<float,3>;
struct Triangle { Point a,b,c; };

void triangle(std::vector<Triangle>& tris,Point a,Point b,Point c)
{
    const double ux=b[0]-a[0],uy=b[1]-a[1],uz=b[2]-a[2];
    const double vx=c[0]-a[0],vy=c[1]-a[1],vz=c[2]-a[2];
    const double nx=uy*vz-uz*vy,ny=uz*vx-ux*vz,nz=ux*vy-uy*vx;
    if(nx*nx+ny*ny+nz*nz>0.0) tris.push_back({a,b,c});
}
void quad(std::vector<Triangle>& tris,Point a,Point b,Point c,Point d)
{
    triangle(tris,a,b,c);triangle(tris,a,c,d);
}

void writeStl(const std::string& path,const std::vector<Triangle>& tris)
{
    std::ofstream f(path,std::ios::binary);
    if(!f) throw std::runtime_error("cannot open STL output: "+path);
    const char header[80]="RFAM STRUCTURED CHAMBER BOUNDARY ONLY; no CAD nozzle or volume cells; units m";
    f.write(header,80);
    const auto n=static_cast<std::uint32_t>(tris.size());
    f.write(reinterpret_cast<const char*>(&n),4);
    for(const auto& t:tris)
    {
        const double ux=t.b[0]-t.a[0],uy=t.b[1]-t.a[1],uz=t.b[2]-t.a[2];
        const double vx=t.c[0]-t.a[0],vy=t.c[1]-t.a[1],vz=t.c[2]-t.a[2];
        const double nx=uy*vz-uz*vy,ny=uz*vx-ux*vz,nz=ux*vy-uy*vx;
        const double mag=std::sqrt(nx*nx+ny*ny+nz*nz);
        const float normal[3]={float(nx/mag),float(ny/mag),float(nz/mag)};
        f.write(reinterpret_cast<const char*>(normal),12);
        for(const auto* p:{&t.a,&t.b,&t.c})
            f.write(reinterpret_cast<const char*>(p->data()),12);
        const std::uint16_t attr=0;
        f.write(reinterpret_cast<const char*>(&attr),2);
    }
    if(!f) throw std::runtime_error("failed writing STL: "+path);
    std::cout<<"BOUNDARY_ONLY "<<path<<" triangles="<<n
             <<" units=m; no volume cells, no injector CAD\n";
}
}

int main(int argc,char** argv)
{
    try
    {
        if(argc!=7 && !(argc==8 && std::string(argv[7])=="--surface"))
        {
            std::cerr<<"usage: export_mesh_stl <profile> <Nx> <Nr> <Ntheta> <degrees> <out.stl> [--surface]\n"
                     <<"Removed --full/--slices/--injector-cad: these combined unconnected surfaces, not a fluid volume mesh.\n"
                     <<"Use rebuild_injector_fluid.py and mesh_injector_tetra.py for CAD fluid boundaries and tetrahedra.\n";
            return 2;
        }
        const int nx=std::stoi(argv[2]),nr=std::stoi(argv[3]),nt=std::stoi(argv[4]);
        const double degrees=std::stod(argv[5]);
        if(!(degrees>0.0 && degrees<=360.0)) throw std::invalid_argument("degrees must be in (0,360]");
        Foam::SectorMeshMetrics mesh(argv[1],nx,nr,nt,degrees);
        const bool full=std::abs(degrees-360.0)<1e-12;
        auto point=[&](int i,int j,int k)->Point
        {
            // Weld the 360-degree seam exactly in exported float coordinates.
            const double t=(full && k==nt-1)?0.0:mesh.thetaFace(k);
            const double r=mesh.wallFace(i)*mesh.etaFace(j);
            return {float(mesh.xFace(i)),float(r*std::cos(t)),float(r*std::sin(t))};
        };
        std::vector<Triangle> tris;
        // Physical boundary only. Ghost cells are numerical storage, not metal.
        for(int i=1;i<nx-1;++i) for(int k=1;k<nt-1;++k)
            quad(tris,point(i,nr-1,k),point(i,nr-1,k+1),point(i+1,nr-1,k+1),point(i+1,nr-1,k));
        for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
        {
            quad(tris,point(1,j,k),point(1,j,k+1),point(1,j+1,k+1),point(1,j+1,k));
            quad(tris,point(nx-1,j,k),point(nx-1,j+1,k),point(nx-1,j+1,k+1),point(nx-1,j,k+1));
        }
        if(!full)
            for(int i=1;i<nx-1;++i) for(int j=1;j<nr-1;++j)
            {
                quad(tris,point(i,j,1),point(i,j+1,1),point(i+1,j+1,1),point(i+1,j,1));
                quad(tris,point(i,j,nt-1),point(i+1,j,nt-1),point(i+1,j+1,nt-1),point(i,j+1,nt-1));
            }
        writeStl(argv[6],tris);
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr<<e.what()<<"\n";
        return 1;
    }
}
