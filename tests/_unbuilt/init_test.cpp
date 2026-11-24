#include <iostream>
#include "src/variables/ConservedVars.H"
#include "src/framework/Domain.H"
#include "src/thermo/SpeciesRegistry.H"
using namespace Foam;
int main(){
    const int N=480; const double L=0.012;
    SpeciesRegistry reg;
    auto species=reg.names();
    Domain dom(N,L,species.size());
    dom.setSpecies(species);
    std::cout << "domain ok\n" << std::flush;
    const int iC3H8=reg.index("C3H8"),iO2=reg.index("O2"),iN2=reg.index("N2"),
              iCO2=reg.index("CO2"),iH2O=reg.index("H2O");
    for(int i=0;i<N;++i){
        PrimitiveVars w(species.size());
        w.p=101325.0; w.T=(dom.x(i)<0.002)?2500.0:350.0;
        if(dom.x(i)<0.002){ w.Y[iCO2]=0.25; w.Y[iH2O]=0.1; w.Y[iN2]=0.65; }
        else { w.Y[iC3H8]=0.1; w.Y[iO2]=0.3; w.Y[iN2]=0.6; }
        w.U[0]=w.U[1]=w.U[2]=0.0;
        dom.stagedW(i)=w;
    }
    std::cout << "init ok\n" << std::flush;
    dom.primitiveToConserved();
    std::cout << "p2c ok, Q(240).rho=" << dom.Q(240).rho << "\n" << std::flush;
    return 0;
}
