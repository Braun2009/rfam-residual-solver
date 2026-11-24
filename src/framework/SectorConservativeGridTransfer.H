#ifndef SectorConservativeGridTransfer_H
#define SectorConservativeGridTransfer_H

/*---------------------------------------------------------------------------*\
    SectorConservativeGridTransfer - aligned conservative sector-grid transfer.

    This is infrastructure for coarse-grid continuation and later geometric
    multigrid. It is deliberately separate from the legacy trilinear restart
    remapper:

      fine -> coarse : volume-integrated conservative restriction;
      coarse -> fine : piecewise-constant conservative prolongation.

    Active grids must cover the same mapped sector and differ by an integer
    ratio in every direction. Halo states are copied from the nearest mapped
    active state; the solver's boundary operators remain authoritative before
    the next residual evaluation. No physical/restart clock is changed.
\*---------------------------------------------------------------------------*/

#include "Domain.H"
#include "SectorMeshMetrics.H"
#include "SectorTurbulenceOperator.H"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Foam
{

class SectorConservativeGridTransfer
{
public:
    enum class Direction { RESTRICT, PROLONG };

    struct Plan
    {
        Direction direction=Direction::RESTRICT;
        int ratioI=0,ratioJ=0,ratioK=0;
    };

private:
    static int activeSize(int total)
    { return total-2; }

    static bool sameGeometry(
        const SectorMeshMetrics& source,
        const SectorMeshMetrics& target)
    {
        const auto close=[](double a,double b)
        {
            return std::abs(a-b)
                <=1.0e-12*std::max({std::abs(a),std::abs(b),1.0});
        };
        return close(source.length(),target.length())
            &&close(source.sectorAngle(),target.sectorAngle());
    }

    static int linearIndex(const Domain& domain,int i,int j,int k)
    { return k*domain.Nx()*domain.Ny()+j*domain.Nx()+i; }

    static void addScaled(
        ConservedVars& target,const ConservedVars& source,double scale)
    {
        target.rho+=scale*source.rho;
        target.rhoE+=scale*source.rhoE;
        for(std::size_t component=0;component<target.rhoU.size();++component)
            target.rhoU[component]+=scale*source.rhoU[component];
        for(std::size_t component=0;component<target.rhoY.size();++component)
            target.rhoY[component]+=scale*source.rhoY[component];
        target.rhoZ+=scale*source.rhoZ;
        target.rhoC+=scale*source.rhoC;
    }

    static void divide(ConservedVars& value,double denominator)
    {
        value.rho/=denominator;
        value.rhoE/=denominator;
        for(double& component:value.rhoU) component/=denominator;
        for(double& component:value.rhoY) component/=denominator;
        value.rhoZ/=denominator;
        value.rhoC/=denominator;
    }

    static void mixPair(
        ConservedVars& left,ConservedVars& right,
        double leftVolume,double rightVolume,double fraction)
    {
        const ConservedVars oldLeft=left;
        const ConservedVars oldRight=right;
        const double totalVolume=leftVolume+rightVolume;
        const double leftWeight=fraction*rightVolume/totalVolume;
        const double rightWeight=fraction*leftVolume/totalVolume;
        auto mix=[&](double a,double b,double& newA,double& newB)
        {
            newA=a+leftWeight*(b-a);
            newB=b+rightWeight*(a-b);
        };
        mix(oldLeft.rho,oldRight.rho,left.rho,right.rho);
        mix(oldLeft.rhoE,oldRight.rhoE,left.rhoE,right.rhoE);
        for(std::size_t component=0;component<left.rhoU.size();++component)
            mix(oldLeft.rhoU[component],oldRight.rhoU[component],
                left.rhoU[component],right.rhoU[component]);
        for(std::size_t component=0;component<left.rhoY.size();++component)
            mix(oldLeft.rhoY[component],oldRight.rhoY[component],
                left.rhoY[component],right.rhoY[component]);
        mix(oldLeft.rhoZ,oldRight.rhoZ,left.rhoZ,right.rhoZ);
        mix(oldLeft.rhoC,oldRight.rhoC,left.rhoC,right.rhoC);
    }

public:
    static Plan plan(const Domain& source,const Domain& target)
    {
        const int sourceI=activeSize(source.Nx());
        const int sourceJ=activeSize(source.Ny());
        const int sourceK=activeSize(source.Nz());
        const int targetI=activeSize(target.Nx());
        const int targetJ=activeSize(target.Ny());
        const int targetK=activeSize(target.Nz());
        if(sourceI<=0||sourceJ<=0||sourceK<=0
            ||targetI<=0||targetJ<=0||targetK<=0)
            throw std::invalid_argument(
                "sector grid transfer requires active cells in every direction");

        const bool restrict=sourceI>=targetI&&sourceJ>=targetJ&&sourceK>=targetK
            &&sourceI%targetI==0&&sourceJ%targetJ==0&&sourceK%targetK==0;
        const bool prolong=targetI>=sourceI&&targetJ>=sourceJ&&targetK>=sourceK
            &&targetI%sourceI==0&&targetJ%sourceJ==0&&targetK%sourceK==0;
        if(restrict&&!(sourceI==targetI&&sourceJ==targetJ&&sourceK==targetK))
            return {Direction::RESTRICT,
                sourceI/targetI,sourceJ/targetJ,sourceK/targetK};
        if(prolong&&!(sourceI==targetI&&sourceJ==targetJ&&sourceK==targetK))
            return {Direction::PROLONG,
                targetI/sourceI,targetJ/sourceJ,targetK/sourceK};
        throw std::invalid_argument(
            "sector grid transfer requires distinct integer-aligned active grids");
    }

    static const char* directionName(Direction direction)
    {
        return direction==Direction::RESTRICT?"restriction":"prolongation";
    }

    static void apply(
        const Domain& source,
        const SectorMeshMetrics& sourceMesh,
        const SectorTurbulenceOperator* sourceTurbulence,
        Domain& target,
        const SectorMeshMetrics& targetMesh,
        SectorTurbulenceOperator* targetTurbulence,
        int prolongationSmoothingSweeps=0)
    {
        if(!source.stateReady()||source.tripped())
            throw std::runtime_error(
                "conservative grid-transfer source state is not ready");
        if(source.stepTransactionActive()||target.stepTransactionActive())
            throw std::runtime_error(
                "conservative grid transfer cannot run during a transaction");
        if(sourceMesh.Nx()!=source.Nx()
            ||sourceMesh.Nr()!=source.Ny()
            ||sourceMesh.Ntheta()!=source.Nz()
            ||targetMesh.Nx()!=target.Nx()
            ||targetMesh.Nr()!=target.Ny()
            ||targetMesh.Ntheta()!=target.Nz())
            throw std::runtime_error(
                "conservative grid-transfer mesh/domain mismatch");
        if(!sameGeometry(sourceMesh,targetMesh))
            throw std::runtime_error(
                "conservative grid-transfer geometry mismatch");
        if(source.nSpecies()!=target.nSpecies()
            ||source.speciesNames()!=target.speciesNames())
            throw std::runtime_error(
                "conservative grid-transfer species mismatch");
        if(source.thermoModelId()!=target.thermoModelId())
            throw std::runtime_error(
                "conservative grid-transfer thermo model mismatch");
        if(static_cast<bool>(sourceTurbulence)
            !=static_cast<bool>(targetTurbulence))
            throw std::runtime_error(
                "conservative grid-transfer turbulence mismatch");

        const Plan transfer=plan(source,target);
        const std::vector<double>* sourceK=nullptr;
        const std::vector<double>* sourceOmega=nullptr;
        if(sourceTurbulence)
        {
            if(sourceTurbulence->stepTransactionActive()
                ||targetTurbulence->stepTransactionActive()
                ||!sourceTurbulence->meshMatches(source)
                ||!targetTurbulence->meshMatches(target))
                throw std::runtime_error(
                    "conservative grid-transfer turbulence contract failed");
            sourceK=&sourceTurbulence->kField();
            sourceOmega=&sourceTurbulence->omegaField();
            if(!sourceTurbulence->fieldsAreRestorable(*sourceK,*sourceOmega))
                throw std::runtime_error(
                    "conservative grid-transfer source turbulence is invalid");
        }

        std::vector<ConservedVars> mapped(
            static_cast<std::size_t>(target.N()),
            ConservedVars(target.nSpecies()));
        std::vector<double> temperature(
            static_cast<std::size_t>(target.N()),300.0);
        std::vector<double> mappedK,mappedOmega;
        if(sourceTurbulence)
        {
            mappedK.assign(static_cast<std::size_t>(target.N()),0.0);
            mappedOmega.assign(static_cast<std::size_t>(target.N()),0.0);
        }

        auto targetId=[&](int i,int j,int k)
        { return static_cast<std::size_t>(linearIndex(target,i,j,k)); };
        auto sourceId=[&](int i,int j,int k)
        { return static_cast<std::size_t>(linearIndex(source,i,j,k)); };

        for(int i=1;i<target.Nx()-1;++i)
        for(int j=1;j<target.Ny()-1;++j)
        for(int k=1;k<target.Nz()-1;++k)
        {
            const std::size_t id=targetId(i,j,k);
            if(transfer.direction==Direction::RESTRICT)
            {
                ConservedVars integrated(target.nSpecies());
                double integratedTemperature=0.0;
                double integratedK=0.0,integratedOmega=0.0;
                double sourceVolumeSum=0.0;
                const int beginI=1+(i-1)*transfer.ratioI;
                const int beginJ=1+(j-1)*transfer.ratioJ;
                const int beginK=1+(k-1)*transfer.ratioK;
                for(int si=beginI;si<beginI+transfer.ratioI;++si)
                for(int sj=beginJ;sj<beginJ+transfer.ratioJ;++sj)
                for(int sk=beginK;sk<beginK+transfer.ratioK;++sk)
                {
                    const double volume=sourceMesh.volume(si,sj,sk);
                    addScaled(integrated,source.U(si,sj,sk),volume);
                    integratedTemperature+=volume*source.W(si,sj,sk).T;
                    sourceVolumeSum+=volume;
                    if(sourceTurbulence)
                    {
                        integratedK+=volume*(*sourceK)[sourceId(si,sj,sk)];
                        integratedOmega+=volume
                            *(*sourceOmega)[sourceId(si,sj,sk)];
                    }
                }
                const double targetVolume=targetMesh.volume(i,j,k);
                // Divide the integrated source inventory by the actual target
                // volume. This preserves every conserved inventory even when
                // mapped wall curvature makes nested geometric volumes differ
                // slightly from an algebraic sum of fine-cell volumes.
                divide(integrated,targetVolume);
                mapped[id]=std::move(integrated);
                temperature[id]=integratedTemperature/sourceVolumeSum;
                if(sourceTurbulence)
                {
                    mappedK[id]=integratedK/sourceVolumeSum;
                    mappedOmega[id]=integratedOmega/sourceVolumeSum;
                }
            }
            else
            {
                const int si=1+(i-1)/transfer.ratioI;
                const int sj=1+(j-1)/transfer.ratioJ;
                const int sk=1+(k-1)/transfer.ratioK;
                double childVolumeSum=0.0;
                const int beginI=1+(si-1)*transfer.ratioI;
                const int beginJ=1+(sj-1)*transfer.ratioJ;
                const int beginK=1+(sk-1)*transfer.ratioK;
                for(int ti=beginI;ti<beginI+transfer.ratioI;++ti)
                for(int tj=beginJ;tj<beginJ+transfer.ratioJ;++tj)
                for(int tk=beginK;tk<beginK+transfer.ratioK;++tk)
                    childVolumeSum+=targetMesh.volume(ti,tj,tk);
                const double inventoryScale=sourceMesh.volume(si,sj,sk)
                    /childVolumeSum;
                mapped[id]=source.U(si,sj,sk);
                mapped[id].rho*=inventoryScale;
                mapped[id].rhoE*=inventoryScale;
                for(double& component:mapped[id].rhoU)
                    component*=inventoryScale;
                for(double& component:mapped[id].rhoY)
                    component*=inventoryScale;
                mapped[id].rhoZ*=inventoryScale;
                mapped[id].rhoC*=inventoryScale;
                temperature[id]=source.W(si,sj,sk).T;
                if(sourceTurbulence)
                {
                    mappedK[id]=(*sourceK)[sourceId(si,sj,sk)];
                    mappedOmega[id]=(*sourceOmega)[sourceId(si,sj,sk)];
                }
            }
        }

        if(transfer.direction==Direction::PROLONG
            &&prolongationSmoothingSweeps>0)
        {
            // Conservative pairwise smoothing removes piecewise-constant
            // parent interfaces before a fine-grid post-smoother sees them.
            // Every pair update is a convex combination and exactly preserves
            // V_left*U_left + V_right*U_right for every conserved component.
            const double fraction=0.25;
            for(int sweep=0;sweep<prolongationSmoothingSweeps;++sweep)
            for(int direction=0;direction<3;++direction)
            for(int parity=0;parity<2;++parity)
            {
                for(int i=1;i<target.Nx()-1;++i)
                for(int j=1;j<target.Ny()-1;++j)
                for(int k=1;k<target.Nz()-1;++k)
                {
                    const int coordinate=direction==0?i:(direction==1?j:k);
                    if((coordinate-1)%2!=parity) continue;
                    int ni=i,nj=j,nk=k;
                    if(direction==0) ++ni;
                    else if(direction==1) ++nj;
                    else ++nk;
                    if(ni>=target.Nx()-1
                        ||nj>=target.Ny()-1
                        ||nk>=target.Nz()-1) continue;
                    const std::size_t left=targetId(i,j,k);
                    const std::size_t right=targetId(ni,nj,nk);
                    const double leftVolume=targetMesh.volume(i,j,k);
                    const double rightVolume=targetMesh.volume(ni,nj,nk);
                    mixPair(mapped[left],mapped[right],
                        leftVolume,rightVolume,fraction);
                    const double totalVolume=leftVolume+rightVolume;
                    const double leftWeight=fraction*rightVolume/totalVolume;
                    const double rightWeight=fraction*leftVolume/totalVolume;
                    const double oldLeftTemperature=temperature[left];
                    const double oldRightTemperature=temperature[right];
                    temperature[left]=oldLeftTemperature+leftWeight
                        *(oldRightTemperature-oldLeftTemperature);
                    temperature[right]=oldRightTemperature+rightWeight
                        *(oldLeftTemperature-oldRightTemperature);
                    if(sourceTurbulence)
                    {
                        const double oldLeftK=mappedK[left];
                        const double oldRightK=mappedK[right];
                        const double oldLeftOmega=mappedOmega[left];
                        const double oldRightOmega=mappedOmega[right];
                        mappedK[left]=oldLeftK+leftWeight*(oldRightK-oldLeftK);
                        mappedK[right]=oldRightK+rightWeight*(oldLeftK-oldRightK);
                        mappedOmega[left]=oldLeftOmega+leftWeight
                            *(oldRightOmega-oldLeftOmega);
                        mappedOmega[right]=oldRightOmega+rightWeight
                            *(oldLeftOmega-oldRightOmega);
                    }
                }
            }
        }

        // Halos are restart seeds only. Use the nearest active transferred
        // state; physical boundary operators overwrite them before residual use.
        for(int i=0;i<target.Nx();++i)
        for(int j=0;j<target.Ny();++j)
        for(int k=0;k<target.Nz();++k)
        {
            if(i>0&&i<target.Nx()-1
                &&j>0&&j<target.Ny()-1
                &&k>0&&k<target.Nz()-1) continue;
            const int ai=std::clamp(i,1,target.Nx()-2);
            const int aj=std::clamp(j,1,target.Ny()-2);
            const int ak=std::clamp(k,1,target.Nz()-2);
            const std::size_t id=targetId(i,j,k);
            const std::size_t activeId=targetId(ai,aj,ak);
            mapped[id]=mapped[activeId];
            temperature[id]=temperature[activeId];
            if(sourceTurbulence)
            {
                mappedK[id]=mappedK[activeId];
                mappedOmega[id]=mappedOmega[activeId];
            }
        }

        if(sourceTurbulence
            &&!targetTurbulence->fieldsAreRestorable(mappedK,mappedOmega))
            throw std::runtime_error(
                "conservative grid-transfer target turbulence is invalid");
        if(!target.replaceConservedField(mapped,temperature))
        {
            const Domain::FailureInfo& failure=target.failure();
            throw std::runtime_error(
                "conservative grid-transfer state recovery failed at cell "
                +std::to_string(failure.cell)+": "+failure.detail);
        }
        if(sourceTurbulence)
            targetTurbulence->restoreFields(
                std::move(mappedK),std::move(mappedOmega));
    }
};

} // End namespace Foam

#endif
