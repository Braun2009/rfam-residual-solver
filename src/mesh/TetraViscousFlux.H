#ifndef RFAM_TETRA_VISCOUS_FLUX_H
#define RFAM_TETRA_VISCOUS_FLUX_H
// Viscous and conductive flux assembly for the tetra finite-volume solver
// (Stage 2, first piece). Deliberately standalone so it can be validated on a
// manufactured field without touching the convective residual path.
//
//   tau    = mu (grad u + grad u^T) - (2/3) mu (div u) I
//   F_mom  = tau . n A
//   F_energy = (tau . u + lambda grad T) . n A
//
// The solver stores the total residual as R = sum(F_conv - F_visc).
#include "TetraMesh.H"
#include "TetraGradients.H"
#include "ConservedVars.H"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Foam
{
struct ViscousFaceFlux
{
    TetraVec3 momentum{0., 0., 0.};   // tau . n A          [N]
    double energy = 0.;               // (tau.u + lambda grad T) . n A  [W]
};

// One face of a structured viscous flux, given the averaged velocity and the
// face gradients of the velocity components and temperature.
inline ViscousFaceFlux viscousFaceFlux(
    const TetraVec3& velocity,
    const std::array<TetraVec3, 3>& gradVelocity,   // grad u_i = (du_i/dx, ...)
    const TetraVec3& gradTemperature,
    double mu, double lambda,
    const std::array<double, 3>& normal, double area)
{
    const double divU = gradVelocity[0][0] + gradVelocity[1][1] + gradVelocity[2][2];
    // tau_ij = mu (du_i/dx_j + du_j/dx_i) - (2/3) mu divU delta_ij
    double tau[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            tau[i][j] = mu*(gradVelocity[i][j] + gradVelocity[j][i])
                      - (i == j ? (2.0/3.0)*mu*divU : 0.0);
    ViscousFaceFlux flux;
    for (int i = 0; i < 3; ++i)
    {
        double traction = 0.0;
        for (int j = 0; j < 3; ++j) traction += tau[i][j]*normal[j];
        flux.momentum[i] = traction*area;
    }
    double work = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        double traction = 0.0;
        for (int j = 0; j < 3; ++j) traction += tau[i][j]*normal[j];
        work += velocity[i]*traction;
    }
    const double conduction = lambda*(gradTemperature[0]*normal[0]
                                     + gradTemperature[1]*normal[1]
                                     + gradTemperature[2]*normal[2]);
    flux.energy = (work + conduction)*area;
    return flux;
}

// Assemble -F_visc into the solver's outward convective residual: the state
// advances with dU/dt = -R/V. mu/lambda are supplied per cell and averaged
// to the face. Boundary wall fluxes are opt-in so the existing manufactured
// interior-flux checks retain their original behavior.
inline void addViscousResidual(
    const TetraMesh& mesh,
    const std::vector<TetraVec3>& velocity,
    const std::vector<double>& temperature,
    const std::vector<double>& mu,
    const std::vector<double>& lambda,
    std::vector<ConservedVars>& residual,
    bool noSlipAdiabaticWalls=false)
{
    const std::size_t cells = mesh.cellCount();
    if (velocity.size() != cells || temperature.size() != cells
        || mu.size() != cells || lambda.size() != cells || residual.size() != cells)
        throw std::invalid_argument("viscous residual field size mismatch");
    TetraGradients gradients(mesh);
    std::array<std::vector<double>, 3> velocityComponent;
    for (int i = 0; i < 3; ++i)
    {
        velocityComponent[i].resize(cells);
        for (std::size_t c = 0; c < cells; ++c) velocityComponent[i][c] = velocity[c][i];
    }
    std::array<std::vector<TetraVec3>, 3> gradVelocity;
    std::array<std::vector<double>, 3> boundaryVelocity;
    if(noSlipAdiabaticWalls)
    {
        for(auto& component:boundaryVelocity)component.resize(mesh.faceCount());
        for(std::size_t f=0;f<mesh.faceCount();++f)
        {
            const auto& face=mesh.face(f);
            if(face.neighbour>=0)continue;
            for(int i=0;i<3;++i)
                boundaryVelocity[i][f]=face.boundaryZone==0?0.:velocity[face.owner][i];
        }
    }
    for (int i = 0; i < 3; ++i)
        gradVelocity[i] = gradients.cellGradient(velocityComponent[i],nullptr,
            noSlipAdiabaticWalls?&boundaryVelocity[i]:nullptr);
    const auto gradTemperature = gradients.cellGradient(temperature);
    for (std::size_t f = 0; f < mesh.faceCount(); ++f)
    {
        const auto& face = mesh.face(f);
        if(face.neighbour<0&&!noSlipAdiabaticWalls)continue;
        if(face.neighbour<0&&face.boundaryZone!=0)continue; // zero-gradient inlet/outlet
        const std::size_t P = static_cast<std::size_t>(face.owner);
        const std::size_t N = face.neighbour<0?P:static_cast<std::size_t>(face.neighbour);
        const auto area = mesh.faceAreaVector(f);
        const double a = std::sqrt(area.x*area.x + area.y*area.y + area.z*area.z);
        if (!(a > 0.)) continue;
        const std::array<double, 3> normal{area.x/a, area.y/a, area.z/a};
        TetraVec3 uFace{0.5*(velocity[P][0] + velocity[N][0]),
                        0.5*(velocity[P][1] + velocity[N][1]),
                        0.5*(velocity[P][2] + velocity[N][2])};
        std::array<TetraVec3, 3> gU;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                gU[i][j] = 0.5*(gradVelocity[i][P][j] + gradVelocity[i][N][j]);
        TetraVec3 gT{0.5*(gradTemperature[P][0] + gradTemperature[N][0]),
                     0.5*(gradTemperature[P][1] + gradTemperature[N][1]),
                     0.5*(gradTemperature[P][2] + gradTemperature[N][2])};
        if(face.neighbour<0)
        {
            const auto xc=mesh.cellCentroid(P),xf=mesh.faceCentroid(f);
            const double d[3]={xf.x-xc.x,xf.y-xc.y,xf.z-xc.z};
            const double dn=d[0]*normal[0]+d[1]*normal[1]+d[2]*normal[2];
            if(!(dn>0.))throw std::runtime_error("invalid wall-normal cell distance");
            uFace={0.,0.,0.};
            for(int i=0;i<3;++i)
            {
                const double current=gU[i][0]*normal[0]+gU[i][1]*normal[1]+gU[i][2]*normal[2];
                const double needed=-velocity[P][i]/dn;
                for(int j=0;j<3;++j)gU[i][j]+=(needed-current)*normal[j];
            }
            gT={0.,0.,0.}; // adiabatic wall; u_wall=0 also makes viscous work zero
        }
        const double muFace = 0.5*(mu[P] + mu[N]), lambdaFace = 0.5*(lambda[P] + lambda[N]);
        if(!(muFace>0.)||!(lambdaFace>0.)||!std::isfinite(muFace)||!std::isfinite(lambdaFace))
            throw std::runtime_error("invalid tetra molecular transport");
        const auto flux = viscousFaceFlux(uFace, gU, gT, muFace, lambdaFace, normal, a);
        // Internal faces cancel exactly; walls exchange momentum with the solid.
        for (int i = 0; i < 3; ++i)
        {
            residual[P].rhoU[i] -= flux.momentum[i];
            if(face.neighbour>=0)residual[N].rhoU[i] += flux.momentum[i];
        }
        residual[P].rhoE -= flux.energy;
        if(face.neighbour>=0)residual[N].rhoE += flux.energy;
    }
}
}
#endif
