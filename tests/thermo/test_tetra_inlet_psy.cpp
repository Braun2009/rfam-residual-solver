#include "CoolPropThermoModel.H"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main()
{
    using namespace Foam;
    if(!thermo::coolPropThermoAvailable())return 77;
    try
    {
        auto model=thermo::makeCoolPropThermoModel();
        const auto supply=model->initializePTY(6e6,130.,{1.,0.});
        if(!supply.success())throw std::runtime_error("oxygen supply flash failed");
        thermo::EquilibriumSeed seed;
        seed.temperature=130.;
        const auto face=model->initializePSY(4995827.5,supply.state->s,{1.,0.},&seed);
        if(!face.success())throw std::runtime_error("oxygen face PS flash failed: "+face.diagnostics.backendStatus);
        if(std::abs(face.state->p-4995827.5)>0.1
           ||std::abs(face.state->s-supply.state->s)>
                1e-8*std::max(std::abs(supply.state->s),1.))
            throw std::runtime_error("oxygen face entropy/pressure mismatch");
        std::cout<<"TETRA_INLET_PSY_PASS\n";
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
