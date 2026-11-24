#include "TetraInletAllocation.H"
#include <cmath>
#include <iostream>
#include <stdexcept>

int main()
{
    using Foam::allocateTetraInletFlow;
    const auto require=[](bool ok,const char* why){if(!ok)throw std::runtime_error(why);};
    try
    {
        const std::vector<double> area{1.,2.,3.};
        const auto uniform=allocateTetraInletFlow(area,{10.,10.,10.},6.);
        require(std::abs(uniform[0]-1.)<1e-12&&std::abs(uniform[1]-2.)<1e-12
            &&std::abs(uniform[2]-3.)<1e-12,"uniform mass flux changed");
        const auto blocked=allocateTetraInletFlow(area,{0.,10.,10.},6.);
        require(blocked[0]==0.&&std::abs(blocked[1]-2.4)<1e-12
            &&std::abs(blocked[2]-3.6)<1e-12,"blocked face was not redistributed");
        const auto capped=allocateTetraInletFlow(area,{0.,1.,10.},6.);
        require(capped[0]==0.&&std::abs(capped[1]-1.)<1e-12
            &&std::abs(capped[2]-5.)<1e-12,"face capacity not enforced");
        const auto shortfall=allocateTetraInletFlow(area,{0.,1.,2.},6.);
        require(shortfall==std::vector<double>({0.,1.,2.}),"shortfall invented supply");
        require(allocateTetraInletFlow(area,{0.,0.,0.},6.)==std::vector<double>({0.,0.,0.}),
                "fully blocked inlet produced flow");
        bool rejected=false;
        try{allocateTetraInletFlow(area,{0.,-1.,2.},1.);}
        catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"negative capacity accepted");
        std::cout<<"TETRA_INLET_ALLOCATION_PASS\n";
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
