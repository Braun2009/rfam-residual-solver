#include "AUSMPlusUpFlux.H"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace Foam;
int main()
{
    try
    {
        AUSMPlusUpFlux ausm({"O2","C3H8"});
        const std::vector<double> n{0.,0.,1.};
        for(double un:{-156.71218858116077,-22.,0.,22.,156.71218858116077})
        {
            PrimitiveVars inside(2);inside.p=4.11e6;inside.T=400.;
            inside.Y={.784,.216};inside.U={8.,9.,un};
            PrimitiveVars ghost=inside;ghost.U[2]=-un;
            ConservedVars flux(2);
            ausm.computeFluxCached(inside,ghost,n,flux,42.4823,353.019,1e5,
                                   42.4823,353.019,1e5);
            if(std::abs(flux.rho)>1e-10 || std::abs(flux.rhoE)>1e-5
                || std::abs(flux.rhoY[0])>1e-10 || std::abs(flux.rhoY[1])>1e-10
                || std::abs(flux.rhoU[0])>1e-8 || std::abs(flux.rhoU[1])>1e-8
                || !std::isfinite(flux.rhoU[2]) || flux.rhoU[2]<=0.)
                throw std::runtime_error("moving slip-wall flux contract failed");
        }
        std::cout<<"TETRA_MOVING_SLIP_WALL_PASS\n";
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
