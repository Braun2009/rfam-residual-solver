#ifndef RFAM_TETRA_BLUEPRINT_WALL_H
#define RFAM_TETRA_BLUEPRINT_WALL_H

#include "TetraMesh.H"
#include "ConservedVars.H"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace Foam
{
// The axial continuity blueprint is not a wall-compatible velocity field.
// Restrict this diagnostic to the nearly horizontal injector plate in the
// ambient gas; inlet-zone priming and all other faces retain their states.
inline std::vector<std::array<double,3>> tetraPlateWallNormals(
    const TetraMesh& mesh,double plateZ,double axialMin=0.8)
{
    std::vector<std::array<double,3>> normals(mesh.cellCount(),{0.,0.,0.});
    for(std::size_t f=0;f<mesh.faceCount();++f)
    {
        const auto& face=mesh.face(f);
        if(face.neighbour>=0||face.boundaryZone!=0
           ||std::abs(mesh.faceCentroid(f).z-plateZ)>5e-6)continue;
        const auto a=mesh.faceAreaVector(f);
        const double mag=std::sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
        if(!(mag>0.)||a.z/mag<axialMin)continue;
        auto& n=normals[face.owner];
        n[0]+=a.x;n[1]+=a.y;n[2]+=a.z;
    }
    for(auto& n:normals)
    {
        const double mag=std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        if(mag>0.)for(auto& component:n)component/=mag;
    }
    return normals;
}

inline bool removeBlueprintPlateImpact(PrimitiveVars& w,
    const std::vector<double>& ambientY,double ambientP,double ambientT,
    const std::array<double,3>& n)
{
    if(w.Y.size()!=ambientY.size()||n[2]<0.8
       ||std::abs(w.p-ambientP)>1e-7*ambientP
       ||std::abs(w.T-ambientT)>1e-7*ambientT)return false;
    for(std::size_t i=0;i<w.Y.size();++i)
        if(std::abs(w.Y[i]-ambientY[i])>1e-7)return false;
    const double un=w.U[0]*n[0]+w.U[1]*n[1]+w.U[2]*n[2];
    if(un<=0.)return false;
    for(int d=0;d<3;++d)w.U[d]-=un*n[d];
    return true;
}
}
#endif
