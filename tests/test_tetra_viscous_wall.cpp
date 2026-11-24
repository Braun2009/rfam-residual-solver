#include "TetraViscousFlux.H"
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
        const auto n=mesh.cellCount();
        const TetraVec3 u{15.,-8.,11.};
        std::vector<TetraVec3> velocity(n,u);
        std::vector<double> temperature(n,400.),mu(n,2e-5),lambda(n,0.05);
        std::vector<ConservedVars> residual(n,ConservedVars(2));
        addViscousResidual(mesh,velocity,temperature,mu,lambda,residual,true);
        double force[3]{},power=0.,totalEnergy=0.,wallArea=0.;
        for(std::size_t c=0;c<n;++c)
        {
            for(int d=0;d<3;++d)force[d]+=residual[c].rhoU[d];
            power+=u[0]*residual[c].rhoU[0]+u[1]*residual[c].rhoU[1]
                  +u[2]*residual[c].rhoU[2];
            totalEnergy+=residual[c].rhoE;
        }
        if(std::abs(totalEnergy)>1e-9)
            throw std::runtime_error("adiabatic domain has net viscous energy flux");
        const std::array<TetraVec3,3> gU{{{0.,0.,-u[0]},{0.,0.,-u[1]},{0.,0.,-u[2]}}};
        const auto wall=viscousFaceFlux({0.,0.,0.},gU,{0.,0.,0.},mu[0],lambda[0],{0.,0.,1.},1.);
        if(wall.energy!=0.)throw std::runtime_error("no-slip adiabatic wall transports energy");
        for(std::size_t f=0;f<mesh.faceCount();++f)
        {
            const auto& face=mesh.face(f);
            if(face.neighbour>=0||face.boundaryZone!=0)continue;
            const auto a=mesh.faceAreaVector(f);
            wallArea+=std::sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
        }
        if(!(wallArea>0.)||!(power>0.)||!std::isfinite(power))
            throw std::runtime_error("wall traction does not oppose moving fluid");
        double before=0.,after=0.;
        double minTime=1e9;
        for(std::size_t c=0;c<n;++c)
        {
            const double rho=50.;
            const double scale=std::max({std::abs(residual[c].rhoU[0]),
                std::abs(residual[c].rhoU[1]),std::abs(residual[c].rhoU[2])});
            if(scale>0.)minTime=std::min(minTime,0.01*rho*mesh.volume(c)/scale);
        }
        const double dt=std::min(minTime,1e-8);

        for(std::size_t c=0;c<n;++c)
        {
            const double rho=50.;
            const auto v=mesh.volume(c);
            before+=.5*rho*(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])*v;
            double newU2=0.;
            for(int d=0;d<3;++d)
            {
                const double next=u[d]-dt*residual[c].rhoU[d]/(rho*v);
                newU2+=next*next;
            }
            after+=.5*rho*newU2*v;
        }
        if(!(after<before))throw std::runtime_error("wall update increased kinetic energy");
        std::cout<<"TETRA_VISCOUS_WALL_PASS power="<<power
                 <<" kinetic_drop_J="<<before-after<<'\n';
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
