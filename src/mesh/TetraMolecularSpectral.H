#ifndef RFAM_TETRA_MOLECULAR_SPECTRAL_H
#define RFAM_TETRA_MOLECULAR_SPECTRAL_H

#include "TetraMesh.H"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Foam
{
// Conservative positive diffusion spectral scale with units [m3/s], the
// same units as sum_f (|u.n|+a)*A in the local Euler CFL and LU-SGS diagonal.
inline std::vector<double> tetraMolecularSpectral(
    const TetraMesh& mesh,const std::vector<double>& rho,
    const std::vector<double>& cv,const std::vector<double>& mu,
    const std::vector<double>& lambda)
{
    const auto cells=mesh.cellCount();
    if(rho.size()!=cells||cv.size()!=cells||mu.size()!=cells||lambda.size()!=cells)
        throw std::invalid_argument("tetra molecular spectral size mismatch");
    std::vector<double> diffusivity(cells),spectral(cells,0.);
    for(std::size_t c=0;c<cells;++c)
    {
        if(!(rho[c]>0.)||!(cv[c]>0.)||!(mu[c]>0.)||!(lambda[c]>0.)
           ||!std::isfinite(rho[c])||!std::isfinite(cv[c])
           ||!std::isfinite(mu[c])||!std::isfinite(lambda[c]))
            throw std::runtime_error("invalid tetra molecular spectral state");
        diffusivity[c]=std::max((4./3.)*mu[c]/rho[c],lambda[c]/(rho[c]*cv[c]));
    }
    for(std::size_t f=0;f<mesh.faceCount();++f)
    {
        const auto& face=mesh.face(f);
        if(face.neighbour<0&&face.boundaryZone!=0)continue; // zero-gradient inlet/outlet
        const auto av=mesh.faceAreaVector(f);
        const double area=std::sqrt(av.x*av.x+av.y*av.y+av.z*av.z);
        if(!(area>0.)||!std::isfinite(area))throw std::runtime_error("invalid tetra molecular face area");
        const auto p=static_cast<std::size_t>(face.owner);
        const auto xp=mesh.cellCentroid(p);
        const auto xn=face.neighbour<0?mesh.faceCentroid(f):mesh.cellCentroid(face.neighbour);
        const double distance=std::abs((xn.x-xp.x)*av.x+(xn.y-xp.y)*av.y+(xn.z-xp.z)*av.z)/area;
        if(!(distance>0.)||!std::isfinite(distance))
            throw std::runtime_error("invalid tetra molecular normal distance");
        const double d=face.neighbour<0?diffusivity[p]
            :std::max(diffusivity[p],diffusivity[face.neighbour]);
        const double contribution=2.*d*area/distance;
        spectral[p]+=contribution;
        if(face.neighbour>=0)spectral[face.neighbour]+=contribution;
    }
    return spectral;
}
}
#endif
