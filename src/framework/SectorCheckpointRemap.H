#ifndef SectorCheckpointRemap_H
#define SectorCheckpointRemap_H

#include "Domain.H"
#include "SectorTurbulenceOperator.H"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Foam
{

class SectorCheckpointRemap
{
    struct Coordinate
    {
        int lower,upper;
        double fraction;
    };

    static Coordinate coordinate(int targetIndex,int targetSize,int sourceSize)
    {
        const double logical=(targetIndex-0.5)/(targetSize-2.0);
        const double source=logical*(sourceSize-2.0)+0.5;
        const int lower=std::clamp(static_cast<int>(std::floor(source)),0,sourceSize-1);
        const int upper=std::clamp(lower+1,0,sourceSize-1);
        return {lower,upper,std::clamp(source-lower,0.0,1.0)};
    }

    template<class Getter>
    static double interpolate(int i,int j,int k,const Domain& source,
                              const Domain& target,Getter getter)
    {
        const Coordinate x=coordinate(i,target.Nx(),source.Nx());
        const Coordinate r=coordinate(j,target.Ny(),source.Ny());
        const Coordinate t=coordinate(k,target.Nz(),source.Nz());
        double value=0.0;
        for(int a=0;a<2;++a) for(int b=0;b<2;++b) for(int c=0;c<2;++c)
        {
            const int si=a?x.upper:x.lower,sj=b?r.upper:r.lower,sk=c?t.upper:t.lower;
            const double weight=(a?x.fraction:1.0-x.fraction)
                *(b?r.fraction:1.0-r.fraction)*(c?t.fraction:1.0-t.fraction);
            value+=weight*getter(si,sj,sk);
        }
        return value;
    }

public:
    static void apply(const Domain& source,const SectorTurbulenceOperator* sourceTurbulence,
                      Domain& target,SectorTurbulenceOperator* targetTurbulence)
    {
        if(!source.stateReady()||source.tripped())
            throw std::runtime_error("checkpoint remap source state is not ready");
        if(source.stepTransactionActive()||target.stepTransactionActive())
            throw std::runtime_error("checkpoint remap cannot run during a step transaction");
        if(source.nSpecies()!=target.nSpecies()
           ||source.speciesNames()!=target.speciesNames())
            throw std::runtime_error("checkpoint remap species mismatch");
        if(source.thermoModelId()!=target.thermoModelId())
            throw std::runtime_error("checkpoint remap thermo model mismatch");
        if(static_cast<bool>(sourceTurbulence)!=static_cast<bool>(targetTurbulence))
            throw std::runtime_error("checkpoint remap turbulence mismatch");
        const std::vector<double>* sourceK=nullptr;
        const std::vector<double>* sourceOmega=nullptr;
        if(sourceTurbulence)
        {
            if(sourceTurbulence->stepTransactionActive()
               ||targetTurbulence->stepTransactionActive())
                throw std::runtime_error(
                    "checkpoint remap turbulence transaction is active");
            if(!sourceTurbulence->meshMatches(source)
               ||!targetTurbulence->meshMatches(target))
                throw std::runtime_error("checkpoint remap turbulence mesh mismatch");
            sourceK=&sourceTurbulence->kField();
            sourceOmega=&sourceTurbulence->omegaField();
            if(!sourceTurbulence->fieldsAreRestorable(*sourceK,*sourceOmega))
                throw std::runtime_error("invalid checkpoint remap turbulence fields");
        }

        std::vector<ConservedVars> mappedField(
            static_cast<std::size_t>(target.N()),
            ConservedVars(target.nSpecies()));
        std::vector<double> temperatureSeeds(
            static_cast<std::size_t>(target.N()),300.0);
        const auto targetId=[&](int i,int j,int k)
        { return static_cast<std::size_t>(k*target.Nx()*target.Ny()+j*target.Nx()+i); };
        for(int i=0;i<target.Nx();++i) for(int j=0;j<target.Ny();++j)
        for(int k=0;k<target.Nz();++k)
        {
            const std::size_t id=targetId(i,j,k);
            ConservedVars& mapped=mappedField[id];
            mapped.rho=interpolate(i,j,k,source,target,
                [&](int a,int b,int c){return source.U(a,b,c).rho;});
            mapped.rhoE=interpolate(i,j,k,source,target,
                [&](int a,int b,int c){return source.U(a,b,c).rhoE;});
            for(int m=0;m<3;++m) mapped.rhoU[m]=interpolate(i,j,k,source,target,
                [&](int a,int b,int c){return source.U(a,b,c).rhoU[m];});
            for(int s=0;s<target.nSpecies();++s) mapped.rhoY[s]=interpolate(
                i,j,k,source,target,[&](int a,int b,int c){return source.U(a,b,c).rhoY[s];});
            mapped.rhoZ=interpolate(i,j,k,source,target,
                [&](int a,int b,int c){return source.U(a,b,c).rhoZ;});
            mapped.rhoC=interpolate(i,j,k,source,target,
                [&](int a,int b,int c){return source.U(a,b,c).rhoC;});
            temperatureSeeds[id]=interpolate(i,j,k,source,target,
                [&](int a,int b,int c){return source.W(a,b,c).T;});
        }
        std::vector<double> mappedK,mappedOmega;
        if(sourceTurbulence)
        {
            mappedK.resize(static_cast<std::size_t>(target.N()));
            mappedOmega.resize(static_cast<std::size_t>(target.N()));
            const auto sourceId=[&](int i,int j,int k)
            { return static_cast<std::size_t>(
                k*source.Nx()*source.Ny()+j*source.Nx()+i); };
            for(int i=0;i<target.Nx();++i) for(int j=0;j<target.Ny();++j)
            for(int k=0;k<target.Nz();++k)
            {
                const std::size_t id=targetId(i,j,k);
                mappedK[id]=interpolate(i,j,k,source,target,[&](int a,int b,int c)
                    {return (*sourceK)[sourceId(a,b,c)];});
                mappedOmega[id]=interpolate(i,j,k,source,target,[&](int a,int b,int c)
                    {return (*sourceOmega)[sourceId(a,b,c)];});
            }
        }

        if(targetTurbulence
           &&!targetTurbulence->fieldsAreRestorable(mappedK,mappedOmega))
            throw std::runtime_error("invalid remapped turbulence fields");
        if(!target.replaceConservedField(mappedField,temperatureSeeds))
        {
            const Domain::FailureInfo& failure=target.failure();
            throw std::runtime_error(
                "checkpoint remap state recovery failed at cell "
                +std::to_string(failure.cell)+": "+failure.detail);
        }
        if(targetTurbulence)
            targetTurbulence->restoreFields(
                std::move(mappedK),std::move(mappedOmega));
    }
};

} // End namespace Foam
#endif
