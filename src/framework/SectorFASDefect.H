#ifndef SectorFASDefect_H
#define SectorFASDefect_H

/*---------------------------------------------------------------------------*\
    SectorFASDefect - aligned conservative FAS defect construction.

    With the steady rate convention used by SteadyResidualEvaluator,

        F(U) = (U_candidate(U)-U)/dt,

    the coarse FAS equation is

        F_H(U_H) + S_H = 0,
        S_H = I_h^H F_h(U_h) - F_H(I_h^H U_h).

    Therefore the forced coarse residual at the restricted fine state is
    exactly the volume-restricted fine residual. Residual rates are restricted
    directly as signed fields; they are not passed through thermodynamic state
    recovery and may contain arbitrary signs.
\*---------------------------------------------------------------------------*/

#include "SectorConservativeGridTransfer.H"
#include "SteadyResidualEvaluator.H"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Foam
{

class SectorFASDefect
{
public:
    struct Result
    {
        std::vector<ConservedVars> restrictedFineResidual;
        std::vector<ConservedVars> tau;
        std::vector<ConservedVars> coarseSource;
        double maximumInitialConsistencyError=0.0;
        double maximumRestrictionInventoryError=0.0;
    };

private:
    static int index(const Domain& domain,int i,int j,int k)
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

    static ConservedVars difference(
        const ConservedVars& left,const ConservedVars& right)
    {
        ConservedVars result(left.rhoY.size());
        result.rho=left.rho-right.rho;
        result.rhoE=left.rhoE-right.rhoE;
        for(std::size_t component=0;component<result.rhoU.size();++component)
            result.rhoU[component]=left.rhoU[component]-right.rhoU[component];
        for(std::size_t component=0;component<result.rhoY.size();++component)
            result.rhoY[component]=left.rhoY[component]-right.rhoY[component];
        result.rhoZ=left.rhoZ-right.rhoZ;
        result.rhoC=left.rhoC-right.rhoC;
        return result;
    }

    static double maximumRelativeDifference(
        const ConservedVars& left,const ConservedVars& right)
    {
        double maximum=0.0;
        auto compare=[&](double a,double b)
        {
            maximum=std::max(maximum,
                std::abs(a-b)/std::max({std::abs(a),std::abs(b),1.0}));
        };
        compare(left.rho,right.rho);
        compare(left.rhoE,right.rhoE);
        for(std::size_t component=0;component<left.rhoU.size();++component)
            compare(left.rhoU[component],right.rhoU[component]);
        for(std::size_t component=0;component<left.rhoY.size();++component)
            compare(left.rhoY[component],right.rhoY[component]);
        compare(left.rhoZ,right.rhoZ);
        compare(left.rhoC,right.rhoC);
        return maximum;
    }

    static ConservedVars activeInventory(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<ConservedVars>& field)
    {
        ConservedVars total(domain.nSpecies());
        for(int i=1;i<domain.Nx()-1;++i)
        for(int j=1;j<domain.Ny()-1;++j)
        for(int k=1;k<domain.Nz()-1;++k)
            addScaled(total,
                field[static_cast<std::size_t>(index(domain,i,j,k))],
                mesh.volume(i,j,k));
        return total;
    }

public:
    static SteadyResidualEvaluator::Result restrictResidual(
        const Domain& fine,
        const SectorMeshMetrics& fineMesh,
        const SteadyResidualEvaluator::Result& fineResidual,
        const Domain& coarse,
        const SectorMeshMetrics& coarseMesh)
    {
        if(fine.nSpecies()!=coarse.nSpecies()
            ||fine.speciesNames()!=coarse.speciesNames()
            ||fineResidual.field.size()!=static_cast<std::size_t>(fine.N())
            ||fineResidual.active.size()!=static_cast<std::size_t>(fine.N()))
            throw std::runtime_error("FAS residual restriction contract mismatch");
        const SectorConservativeGridTransfer::Plan transfer=
            SectorConservativeGridTransfer::plan(fine,coarse);
        if(transfer.direction
            !=SectorConservativeGridTransfer::Direction::RESTRICT)
            throw std::runtime_error(
                "FAS residual restriction requires an aligned fine-to-coarse grid pair");

        SteadyResidualEvaluator::Result result;
        result.referenceStep=fineResidual.referenceStep;
        result.field.assign(
            static_cast<std::size_t>(coarse.N()),
            ConservedVars(coarse.nSpecies()));
        result.active.assign(static_cast<std::size_t>(coarse.N()),0U);
        for(int ci=1;ci<coarse.Nx()-1;++ci)
        for(int cj=1;cj<coarse.Ny()-1;++cj)
        for(int ck=1;ck<coarse.Nz()-1;++ck)
        {
            const int coarseCell=index(coarse,ci,cj,ck);
            ConservedVars restricted(coarse.nSpecies());
            const int beginI=1+(ci-1)*transfer.ratioI;
            const int beginJ=1+(cj-1)*transfer.ratioJ;
            const int beginK=1+(ck-1)*transfer.ratioK;
            for(int fi=beginI;fi<beginI+transfer.ratioI;++fi)
            for(int fj=beginJ;fj<beginJ+transfer.ratioJ;++fj)
            for(int fk=beginK;fk<beginK+transfer.ratioK;++fk)
            {
                const int fineCell=index(fine,fi,fj,fk);
                if(fineResidual.active[static_cast<std::size_t>(fineCell)]==0U)
                    throw std::runtime_error(
                        "FAS fine residual excludes a restricted child cell");
                addScaled(restricted,
                    fineResidual.field[static_cast<std::size_t>(fineCell)],
                    fineMesh.volume(fi,fj,fk));
            }
            divide(restricted,coarseMesh.volume(ci,cj,ck));
            result.field[static_cast<std::size_t>(coarseCell)]=
                std::move(restricted);
            result.active[static_cast<std::size_t>(coarseCell)]=1U;
        }
        return result;
    }

    static Result build(
        const Domain& fine,
        const SectorMeshMetrics& fineMesh,
        const SteadyResidualEvaluator::Result& fineResidual,
        const Domain& coarse,
        const SectorMeshMetrics& coarseMesh,
        const SteadyResidualEvaluator::Result& coarseResidual)
    {
        if(coarseResidual.field.size()!=static_cast<std::size_t>(coarse.N())
            ||coarseResidual.active.size()!=static_cast<std::size_t>(coarse.N()))
            throw std::runtime_error("FAS coarse residual field contract mismatch");
        const SteadyResidualEvaluator::Result restrictedFine=
            restrictResidual(
                fine,fineMesh,fineResidual,coarse,coarseMesh);

        Result result;
        result.restrictedFineResidual=restrictedFine.field;
        result.tau.assign(
            static_cast<std::size_t>(coarse.N()),
            ConservedVars(coarse.nSpecies()));
        result.coarseSource.assign(
            static_cast<std::size_t>(coarse.N()),
            ConservedVars(coarse.nSpecies()));

        for(int ci=1;ci<coarse.Nx()-1;++ci)
        for(int cj=1;cj<coarse.Ny()-1;++cj)
        for(int ck=1;ck<coarse.Nz()-1;++ck)
        {
            const int coarseCell=index(coarse,ci,cj,ck);
            if(coarseResidual.active[static_cast<std::size_t>(coarseCell)]==0U)
                throw std::runtime_error(
                    "FAS coarse residual excludes an active coarse cell");
            const ConservedVars& restricted=
                result.restrictedFineResidual[
                    static_cast<std::size_t>(coarseCell)];
            result.tau[static_cast<std::size_t>(coarseCell)]=difference(
                coarseResidual.field[static_cast<std::size_t>(coarseCell)],
                restricted);
            result.coarseSource[static_cast<std::size_t>(coarseCell)]=difference(
                restricted,
                coarseResidual.field[static_cast<std::size_t>(coarseCell)]);

            ConservedVars forced=
                coarseResidual.field[static_cast<std::size_t>(coarseCell)];
            addScaled(forced,
                result.coarseSource[static_cast<std::size_t>(coarseCell)],1.0);
            result.maximumInitialConsistencyError=std::max(
                result.maximumInitialConsistencyError,
                maximumRelativeDifference(forced,restricted));
        }

        const ConservedVars fineInventory=
            activeInventory(fine,fineMesh,fineResidual.field);
        const ConservedVars restrictedInventory=
            activeInventory(coarse,coarseMesh,result.restrictedFineResidual);
        result.maximumRestrictionInventoryError=
            maximumRelativeDifference(fineInventory,restrictedInventory);
        return result;
    }
};

} // End namespace Foam

#endif
