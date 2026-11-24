#include <cmath>
#include <iostream>
#include "src/thermo/RealFluidState.H"
#include "src/thermo/SpeciesRegistry.H"

using namespace Foam;

int main()
{
    SpeciesRegistry registry;
    std::vector<double> y(registry.size(),0.0);
    y[registry.index("CO2")]=0.6468; y[registry.index("H2O")]=0.3532;
    RealFluidState state; state.setSpecies(registry.names());
    const auto status=state.updatePTY(4e6,700.0,y);
    std::cout<<"BACKGROUND status="<<static_cast<int>(status)<<" rho="<<state.rho()
             <<" e="<<state.e()<<" Cv="<<state.Cv()<<" p="<<state.p()
             <<" T="<<state.T()<<"\n";
    if(status!=RealFluidState::Status::SUCCESS||!std::isfinite(state.rho())
       ||!std::isfinite(state.e())||!std::isfinite(state.Cv())||state.rho()<=0.0)
        return 1;
    return 0;
}
