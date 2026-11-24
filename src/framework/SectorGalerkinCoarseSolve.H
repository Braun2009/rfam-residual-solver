#ifndef SectorGalerkinCoarseSolve_H
#define SectorGalerkinCoarseSolve_H

/*---------------------------------------------------------------------------*\
    SectorGalerkinCoarseSolve - stationary coarse solve preconditioned by the
    same forward/backward approximate factorisation as steady LU-SGS.

    For the positive scalar RAP operator

        A x = D x - sum_f C_f x_f,

    each sweep evaluates r=b-Ax and applies

        (D+L) D^-1 (D+U) delta = r,
        x <- x + delta.

    Every conserved component shares the same positive stencil. Halo
    corrections are fixed at zero, matching the steady active-unknown contract.
\*---------------------------------------------------------------------------*/

#include "Operator.H"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Foam
{

class SectorGalerkinCoarseSolve
{
public:
    struct Result
    {
        std::vector<ConservedVars> correction;
        int sweeps=0;
        double initialResidualNorm=0.0;
        double finalResidualNorm=0.0;
        double residualRatio=1.0;
    };

private:
    static bool active(const Domain& domain,int cell)
    {
        if(cell<0||cell>=domain.N()) return false;
        const int nx=domain.Nx();
        const int ny=domain.Ny();
        const int plane=nx*ny;
        const int i=cell%nx;
        const int j=(cell/nx)%ny;
        const int k=cell/plane;
        return i>0&&i<nx-1&&j>0&&j<ny-1&&k>0&&k<domain.Nz()-1;
    }

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

    static double squaredNorm(const ConservedVars& value)
    {
        double result=value.rho*value.rho+value.rhoE*value.rhoE
            +value.rhoZ*value.rhoZ+value.rhoC*value.rhoC;
        for(double component:value.rhoU) result+=component*component;
        for(double component:value.rhoY) result+=component*component;
        return result;
    }

    static double fieldNorm(
        const Domain& domain,const std::vector<ConservedVars>& field)
    {
        double sum=0.0;
        std::size_t count=0;
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,cell)) continue;
            sum+=squaredNorm(field[static_cast<std::size_t>(cell)]);
            count+=static_cast<std::size_t>(7+domain.nSpecies());
        }
        return std::sqrt(sum/std::max<std::size_t>(count,1U));
    }

    static std::vector<ConservedVars> residual(
        const Domain& domain,
        const StructuredImplicitStencil& stencil,
        double pseudoFactor,
        const std::vector<ConservedVars>& rightHandSide,
        const std::vector<ConservedVars>& correction)
    {
        const int nx=domain.Nx();
        const int plane=nx*domain.Ny();
        std::vector<ConservedVars> result=rightHandSide;
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,cell))
            {
                result[static_cast<std::size_t>(cell)]=
                    ConservedVars(domain.nSpecies());
                continue;
            }
            addScaled(result[static_cast<std::size_t>(cell)],
                correction[static_cast<std::size_t>(cell)],
                -pseudoFactor*stencil.diagonal[
                    static_cast<std::size_t>(cell)]);
            const auto addNeighbor=[&](int neighbor,int face)
            {
                if(!active(domain,neighbor)) return;
                addScaled(result[static_cast<std::size_t>(cell)],
                    correction[static_cast<std::size_t>(neighbor)],
                    stencil.face[face][static_cast<std::size_t>(cell)]);
            };
            addNeighbor(cell-1,StructuredImplicitStencil::I_MINUS);
            addNeighbor(cell+1,StructuredImplicitStencil::I_PLUS);
            addNeighbor(cell-nx,StructuredImplicitStencil::J_MINUS);
            addNeighbor(cell+nx,StructuredImplicitStencil::J_PLUS);
            addNeighbor(cell-plane,StructuredImplicitStencil::K_MINUS);
            addNeighbor(cell+plane,StructuredImplicitStencil::K_PLUS);
        }
        return result;
    }

public:
    static Result solve(
        const Domain& domain,
        const StructuredImplicitStencil& stencil,
        double pseudoCfl,
        const std::vector<ConservedVars>& rightHandSide,
        int maximumSweeps=20,
        double relativeTolerance=1.0e-8)
    {
        if(!std::isfinite(pseudoCfl)||pseudoCfl<=0.0
            ||maximumSweeps<=0
            ||!std::isfinite(relativeTolerance)
            ||relativeTolerance<=0.0||relativeTolerance>=1.0)
            throw std::invalid_argument(
                "invalid Galerkin coarse solve controls");
        if(stencil.diagonal.size()!=static_cast<std::size_t>(domain.N())
            ||rightHandSide.size()!=static_cast<std::size_t>(domain.N()))
            throw std::runtime_error(
                "Galerkin coarse solve dimensions do not match Domain");
        for(const auto& face:stencil.face)
            if(face.size()!=static_cast<std::size_t>(domain.N()))
                throw std::runtime_error(
                    "Galerkin coarse solve face dimensions do not match Domain");
        for(const ConservedVars& value:rightHandSide)
            if(value.rhoY.size()!=static_cast<std::size_t>(domain.nSpecies()))
                throw std::runtime_error(
                    "Galerkin coarse right-hand-side species mismatch");

        const int nx=domain.Nx();
        const int plane=nx*domain.Ny();
        const double pseudoFactor=1.0+1.0/pseudoCfl;
        std::vector<double> diagonal(
            static_cast<std::size_t>(domain.N()),0.0);
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,cell)) continue;
            const double value=pseudoFactor*stencil.diagonal[
                static_cast<std::size_t>(cell)];
            if(!std::isfinite(value)||value<=0.0)
                throw std::runtime_error(
                    "Galerkin coarse solve diagonal is not positive");
            diagonal[static_cast<std::size_t>(cell)]=value;
        }

        Result result;
        result.correction.assign(
            static_cast<std::size_t>(domain.N()),
            ConservedVars(domain.nSpecies()));
        std::vector<ConservedVars> currentResidual=residual(
            domain,stencil,pseudoFactor,rightHandSide,result.correction);
        result.initialResidualNorm=fieldNorm(domain,currentResidual);
        if(!std::isfinite(result.initialResidualNorm))
            throw std::runtime_error(
                "Galerkin coarse initial residual norm is not finite");
        if(result.initialResidualNorm==0.0)
        {
            result.finalResidualNorm=0.0;
            result.residualRatio=0.0;
            return result;
        }

        std::vector<ConservedVars> forward(
            static_cast<std::size_t>(domain.N()),
            ConservedVars(domain.nSpecies()));
        std::vector<ConservedVars> increment(
            static_cast<std::size_t>(domain.N()),
            ConservedVars(domain.nSpecies()));
        for(int sweep=0;sweep<maximumSweeps;++sweep)
        {
            for(int cell=0;cell<domain.N();++cell)
            {
                forward[static_cast<std::size_t>(cell)]=
                    ConservedVars(domain.nSpecies());
                increment[static_cast<std::size_t>(cell)]=
                    ConservedVars(domain.nSpecies());
            }
            for(int cell=0;cell<domain.N();++cell)
            {
                if(!active(domain,cell)) continue;
                ConservedVars value=
                    currentResidual[static_cast<std::size_t>(cell)];
                const auto lower=[&](int neighbor,int face)
                {
                    if(!active(domain,neighbor)) return;
                    addScaled(value,
                        forward[static_cast<std::size_t>(neighbor)],
                        stencil.face[face][static_cast<std::size_t>(cell)]);
                };
                lower(cell-1,StructuredImplicitStencil::I_MINUS);
                lower(cell-nx,StructuredImplicitStencil::J_MINUS);
                lower(cell-plane,StructuredImplicitStencil::K_MINUS);
                divide(value,diagonal[static_cast<std::size_t>(cell)]);
                forward[static_cast<std::size_t>(cell)]=std::move(value);
            }
            for(int cell=domain.N()-1;cell>=0;--cell)
            {
                if(!active(domain,cell)) continue;
                ConservedVars value=
                    forward[static_cast<std::size_t>(cell)];
                const auto upper=[&](int neighbor,int face)
                {
                    if(!active(domain,neighbor)) return;
                    addScaled(value,
                        increment[static_cast<std::size_t>(neighbor)],
                        stencil.face[face][static_cast<std::size_t>(cell)]
                            /diagonal[static_cast<std::size_t>(cell)]);
                };
                upper(cell+1,StructuredImplicitStencil::I_PLUS);
                upper(cell+nx,StructuredImplicitStencil::J_PLUS);
                upper(cell+plane,StructuredImplicitStencil::K_PLUS);
                increment[static_cast<std::size_t>(cell)]=std::move(value);
            }
            for(int cell=0;cell<domain.N();++cell)
                if(active(domain,cell))
                    addScaled(result.correction[
                        static_cast<std::size_t>(cell)],
                        increment[static_cast<std::size_t>(cell)],1.0);

            currentResidual=residual(
                domain,stencil,pseudoFactor,
                rightHandSide,result.correction);
            result.finalResidualNorm=fieldNorm(domain,currentResidual);
            result.sweeps=sweep+1;
            result.residualRatio=result.finalResidualNorm
                /result.initialResidualNorm;
            if(!std::isfinite(result.finalResidualNorm))
                throw std::runtime_error(
                    "Galerkin coarse residual norm is not finite");
            if(result.residualRatio<=relativeTolerance) break;
        }
        return result;
    }
};

} // End namespace Foam

#endif
