#include "TetraMolecularSpectral.H"
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
        std::vector<double> rho(n,42.48),cv(n,700.),mu(n,1.8e-5),lambda(n,0.03);
        const auto base=tetraMolecularSpectral(mesh,rho,cv,mu,lambda);
        for(auto& v:mu)v*=10.;
        for(auto& v:lambda)v*=10.;
        const auto strong=tetraMolecularSpectral(mesh,rho,cv,mu,lambda);
        for(std::size_t c=0;c<n;++c)
            if(!(base[c]>0.)||!std::isfinite(base[c])
                ||std::abs(strong[c]/base[c]-10.)>1e-12)
                throw std::runtime_error("molecular spectral positivity/scaling failed");
        cv[0]=0.;bool rejected=false;
        try{tetraMolecularSpectral(mesh,rho,cv,mu,lambda);}
        catch(const std::runtime_error&){rejected=true;}
        if(!rejected)throw std::runtime_error("invalid molecular heat capacity admitted");
        std::cout<<"TETRA_MOLECULAR_SPECTRAL_PASS cells="<<n<<'\n';
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
