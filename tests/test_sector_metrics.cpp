#include <cmath>
#include <iostream>

#include "src/mesh/SectorMeshMetrics.H"

using namespace Foam;

int main()
{
    const std::string profile = vortexProfilePath();
    SectorMeshMetrics mesh(profile, 120, 30, 20, 120.0);

    if (std::abs(mesh.length() - 0.4065) > 1e-9) return 1;
    if (std::abs(mesh.xFace(1)) > 1e-14) return 1;
    if (std::abs(mesh.xFace(mesh.Nx()-1)-mesh.length()) > 1e-14) return 1;
    if (std::abs(mesh.wallFace(1) - 0.024) > 1e-9) return 1;
    if (std::abs(mesh.wallFace(mesh.Nx()-1) - 0.060) > 1e-9) return 1;
    if (std::abs(mesh.etaFace(1)) > 1e-14) return 1;
    if (std::abs(mesh.etaFace(mesh.Nr()-1)-1.0) > 1e-14) return 1;
    if (std::abs(mesh.thetaFace(1)) > 1e-14) return 1;
    if (std::abs(mesh.thetaFace(mesh.Ntheta()-1)-mesh.sectorAngle()) > 1e-14) return 1;

    // Boundary operators need a real halo outside the wall and an axis
    // mirror halo. Clamping these faces to 0/1 silently removes both layers.
    if (!(mesh.etaFace(0)<0.0 && mesh.etaFace(mesh.Nr())>1.0))
    {
        std::cerr << "radial halo collapsed onto physical boundary\n";
        return 1;
    }
    for (int i=0;i<mesh.Nx();++i)
        if (!(mesh.xFace(i+1)>mesh.xFace(i))) return 1;
    for (int j=0;j<mesh.Nr();++j)
        if (!(mesh.etaFace(j+1)>mesh.etaFace(j))) return 1;
    for (int k=0;k<mesh.Ntheta();++k)
        if (!(mesh.thetaFace(k+1)>mesh.thetaFace(k))) return 1;
    const double nearWall=mesh.radiusCenter(10,mesh.Nr()-2);
    const double ghost=mesh.radiusCenter(10,mesh.Nr()-1);
    if (!(nearWall<mesh.wallCenter(10) && ghost>mesh.wallCenter(10))) return 1;

    double axialIntegral = 0.0;
    for (int i = 1; i < mesh.Nx()-1; ++i)
    {
        const double dx = mesh.xFace(i+1) - mesh.xFace(i);
        const double rw = mesh.wallCenter(i);
        axialIntegral += 0.5 * mesh.sectorAngle() * rw * rw * dx;
    }
    const double total = mesh.totalVolume();
    const double rel = std::abs(total - axialIntegral) / axialIntegral;
    if (rel > 1e-12)
    {
        std::cerr << "sector volume closure failed: " << rel << "\n";
        return 1;
    }

    for (int i = 1; i < mesh.Nx()-1; ++i)
    for (int j = 1; j < mesh.Nr()-1; ++j)
    for (int k = 1; k < mesh.Ntheta()-1; ++k)
    {
        if (mesh.volume(i,j,k) <= 0.0 ||
            mesh.axialArea(i,j,k) <= 0.0 ||
            mesh.azimuthalArea(i,j,k) <= 0.0)
            return 1;
    }

    SectorMeshMetrics full(profile,120,30,56,360.0);
    if (std::abs(full.totalVolume()/total-3.0)>1e-11) return 1;
    // Full-annulus seam must identify exactly the same physical locations.
    const double r=full.wallFace(10)*full.etaFace(10);
    const double a=full.thetaFace(1),b=full.thetaFace(full.Ntheta()-1);
    if (std::hypot(r*(std::cos(a)-std::cos(b)),
                   r*(std::sin(a)-std::sin(b)))>1e-14) return 1;

    std::cout << "SECTOR_METRICS_PASS volume=" << total << " m3\n";
    return 0;
}
