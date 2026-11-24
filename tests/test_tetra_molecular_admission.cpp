#include "TetraMolecularAdmission.H"
#include <iostream>
#include <stdexcept>

int main()
{
    using namespace Foam;
    using Region=thermo::EquilibriumRegion;
    try
    {
        for(auto region:{Region::SINGLE_LIQUID,Region::SINGLE_VAPOR,
                         Region::SUPERCRITICAL,Region::UNCLASSIFIED_SINGLE_PHASE})
        {
            if(!tetraAdmitsMolecularTransport(region,-1.))
                throw std::runtime_error("single phase was rejected");
            if(tetraAdmitsMolecularTransport(region,0.5))
                throw std::runtime_error("two-phase flash was accepted");
        }
        for(auto region:{Region::TWO_PHASE,Region::UNKNOWN,Region::CONDENSED_SOLID})
            if(tetraAdmitsMolecularTransport(region,-1.))
                throw std::runtime_error("unvalidated transport region was admitted");
        std::cout<<"TETRA_MOLECULAR_ADMISSION_PASS\n";
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
