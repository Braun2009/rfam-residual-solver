#include "TetraBlueprintWall.H"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace Foam;
int main(int argc,char**argv)
{
    if(argc!=2)return 2;
    try
    {
        TetraMesh mesh;mesh.read(argv[1]);
        const auto normals=tetraPlateWallNormals(mesh,0.2709527893066406);
        PrimitiveVars ambient(2);
        ambient.p=4.11e6;ambient.T=400.;ambient.Y={.784,.216};
        ambient.U={0.,0.,156.71218858116077};
        auto w=ambient;
        if(!removeBlueprintPlateImpact(w,ambient.Y,ambient.p,ambient.T,normals.at(398377)))
            throw std::runtime_error("plate wall cell did not receive projection");
        if(std::abs(w.U[2])>1e-8||w.p!=ambient.p||w.T!=ambient.T||w.Y!=ambient.Y)
            throw std::runtime_error("projection changed thermodynamics or left wall penetration");
        w=ambient;
        if(removeBlueprintPlateImpact(w,ambient.Y,ambient.p,ambient.T,normals.at(1741581)))
            throw std::runtime_error("LOX inlet cell was mistaken for ambient plate");
        w=ambient;w.Y={1.,0.};
        if(removeBlueprintPlateImpact(w,ambient.Y,ambient.p,ambient.T,normals.at(398377)))
            throw std::runtime_error("propellant region was changed");
        w=ambient;w.U[2]=-10.;
        if(removeBlueprintPlateImpact(w,ambient.Y,ambient.p,ambient.T,normals.at(398377)))
            throw std::runtime_error("wall-outgoing velocity was changed");
        std::cout<<"TETRA_BLUEPRINT_WALL_PASS\n";
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
