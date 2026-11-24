#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../src/framework/SectorGalerkinCoarseSolve.H"
#include "../src/framework/SectorGalerkinCoarseStencil.H"
#include "../src/mesh/SectorMeshMetrics.H"

using namespace Foam;

namespace
{

int index(const Domain& domain,int i,int j,int k)
{ return k*domain.Nx()*domain.Ny()+j*domain.Nx()+i; }

bool active(const Domain& domain,int i,int j,int k)
{
    return i>0&&i<domain.Nx()-1
        &&j>0&&j<domain.Ny()-1
        &&k>0&&k<domain.Nz()-1;
}

} // namespace

int main()
{
    const std::string profile=vortexProfilePath();
    SectorMeshMetrics fineMesh(profile,10,8,6,120.0);
    SectorMeshMetrics coarseMesh(profile,6,5,4,120.0);
    Domain fine(10,8,6,fineMesh.length(),0.065,
        fineMesh.sectorAngle(),2);
    Domain coarse(6,5,4,coarseMesh.length(),0.065,
        coarseMesh.sectorAngle(),2);

    StructuredImplicitStencil fineStencil;
    fineStencil.ensureSize(static_cast<std::size_t>(fine.N()));
    for(int i=1;i<fine.Nx()-1;++i)
    for(int j=1;j<fine.Ny()-1;++j)
    for(int k=1;k<fine.Nz()-1;++k)
    {
        const int cell=index(fine,i,j,k);
        const double spatial=1.0+0.03*i+0.02*j+0.01*k;
        fineStencil.diagonal[static_cast<std::size_t>(cell)]=20.0*spatial;
        fineStencil.face[StructuredImplicitStencil::I_MINUS]
            [static_cast<std::size_t>(cell)]=0.8+0.01*j;
        fineStencil.face[StructuredImplicitStencil::I_PLUS]
            [static_cast<std::size_t>(cell)]=0.9+0.01*k;
        fineStencil.face[StructuredImplicitStencil::J_MINUS]
            [static_cast<std::size_t>(cell)]=0.7+0.01*i;
        fineStencil.face[StructuredImplicitStencil::J_PLUS]
            [static_cast<std::size_t>(cell)]=0.6+0.02*k;
        fineStencil.face[StructuredImplicitStencil::K_MINUS]
            [static_cast<std::size_t>(cell)]=0.5+0.01*i;
        fineStencil.face[StructuredImplicitStencil::K_PLUS]
            [static_cast<std::size_t>(cell)]=0.4+0.02*j;
    }

    const double pseudoCfl=2.5;
    const double pseudoFactor=1.0+1.0/pseudoCfl;
    const SectorGalerkinCoarseStencil::Result aggregated=
        SectorGalerkinCoarseStencil::aggregate(
            fine,fineMesh,fineStencil,coarse,coarseMesh,pseudoCfl);
    if(aggregated.activeCoarseCells!=24
        ||!(aggregated.minimumActualDiagonal>0.0)
        ||!(aggregated.minimumDiagonalDominance>0.0)
        ||!(aggregated.maximumFaceCoupling>0.0))
    {
        std::cerr << "invalid Galerkin aggregate diagnostics\n";
        return 1;
    }

    std::vector<double> coarseCorrection(
        static_cast<std::size_t>(coarse.N()),0.0);
    for(int i=1;i<coarse.Nx()-1;++i)
    for(int j=1;j<coarse.Ny()-1;++j)
    for(int k=1;k<coarse.Nz()-1;++k)
        coarseCorrection[static_cast<std::size_t>(index(coarse,i,j,k))]=
            std::sin(0.4*i)+0.2*j-0.15*k;

    std::vector<double> fineCorrection(
        static_cast<std::size_t>(fine.N()),0.0);
    for(int i=1;i<fine.Nx()-1;++i)
    for(int j=1;j<fine.Ny()-1;++j)
    for(int k=1;k<fine.Nz()-1;++k)
    {
        const int ci=1+(i-1)/2;
        const int cj=1+(j-1)/2;
        const int ck=1+(k-1)/2;
        fineCorrection[static_cast<std::size_t>(index(fine,i,j,k))]=
            coarseCorrection[static_cast<std::size_t>(
                index(coarse,ci,cj,ck))];
    }

    std::vector<double> fineAction(
        static_cast<std::size_t>(fine.N()),0.0);
    for(int i=1;i<fine.Nx()-1;++i)
    for(int j=1;j<fine.Ny()-1;++j)
    for(int k=1;k<fine.Nz()-1;++k)
    {
        const int cell=index(fine,i,j,k);
        double value=pseudoFactor*fineStencil.diagonal[
            static_cast<std::size_t>(cell)]
            *fineCorrection[static_cast<std::size_t>(cell)];
        const auto subtract=[&](int face,int ni,int nj,int nk)
        {
            value-=fineStencil.face[face][static_cast<std::size_t>(cell)]
                *fineCorrection[static_cast<std::size_t>(
                    index(fine,ni,nj,nk))];
        };
        subtract(StructuredImplicitStencil::I_MINUS,i-1,j,k);
        subtract(StructuredImplicitStencil::I_PLUS,i+1,j,k);
        subtract(StructuredImplicitStencil::J_MINUS,i,j-1,k);
        subtract(StructuredImplicitStencil::J_PLUS,i,j+1,k);
        subtract(StructuredImplicitStencil::K_MINUS,i,j,k-1);
        subtract(StructuredImplicitStencil::K_PLUS,i,j,k+1);
        fineAction[static_cast<std::size_t>(cell)]=value;
    }

    double maximumRelativeError=0.0;
    for(int ci=1;ci<coarse.Nx()-1;++ci)
    for(int cj=1;cj<coarse.Ny()-1;++cj)
    for(int ck=1;ck<coarse.Nz()-1;++ck)
    {
        double restricted=0.0;
        for(int fi=1+(ci-1)*2;fi<1+ci*2;++fi)
        for(int fj=1+(cj-1)*2;fj<1+cj*2;++fj)
        for(int fk=1+(ck-1)*2;fk<1+ck*2;++fk)
            restricted+=fineMesh.volume(fi,fj,fk)
                *fineAction[static_cast<std::size_t>(index(fine,fi,fj,fk))];
        restricted/=coarseMesh.volume(ci,cj,ck);

        const int cell=index(coarse,ci,cj,ck);
        double coarseAction=pseudoFactor*aggregated.stencil.diagonal[
            static_cast<std::size_t>(cell)]
            *coarseCorrection[static_cast<std::size_t>(cell)];
        const auto subtract=[&](int face,int ni,int nj,int nk)
        {
            coarseAction-=aggregated.stencil.face[face]
                [static_cast<std::size_t>(cell)]
                *coarseCorrection[static_cast<std::size_t>(
                    index(coarse,ni,nj,nk))];
        };
        subtract(StructuredImplicitStencil::I_MINUS,ci-1,cj,ck);
        subtract(StructuredImplicitStencil::I_PLUS,ci+1,cj,ck);
        subtract(StructuredImplicitStencil::J_MINUS,ci,cj-1,ck);
        subtract(StructuredImplicitStencil::J_PLUS,ci,cj+1,ck);
        subtract(StructuredImplicitStencil::K_MINUS,ci,cj,ck-1);
        subtract(StructuredImplicitStencil::K_PLUS,ci,cj,ck+1);
        maximumRelativeError=std::max(maximumRelativeError,
            std::abs(coarseAction-restricted)
            /std::max({std::abs(coarseAction),std::abs(restricted),1.0}));
    }

    if(maximumRelativeError>3.0e-14)
    {
        std::cerr << "Galerkin RAP identity error="
                  << maximumRelativeError << "\n";
        return 1;
    }

    std::vector<ConservedVars> rightHandSide(
        static_cast<std::size_t>(coarse.N()),ConservedVars(2));
    for(int ci=1;ci<coarse.Nx()-1;++ci)
    for(int cj=1;cj<coarse.Ny()-1;++cj)
    for(int ck=1;ck<coarse.Nz()-1;++ck)
    {
        const int cell=index(coarse,ci,cj,ck);
        ConservedVars& value=rightHandSide[static_cast<std::size_t>(cell)];
        const double scalar=coarseCorrection[static_cast<std::size_t>(cell)];
        value.rho=pseudoFactor*aggregated.stencil.diagonal[
            static_cast<std::size_t>(cell)]*scalar;
        const auto subtract=[&](int face,int ni,int nj,int nk)
        {
            value.rho-=aggregated.stencil.face[face]
                [static_cast<std::size_t>(cell)]
                *coarseCorrection[static_cast<std::size_t>(
                    index(coarse,ni,nj,nk))];
        };
        subtract(StructuredImplicitStencil::I_MINUS,ci-1,cj,ck);
        subtract(StructuredImplicitStencil::I_PLUS,ci+1,cj,ck);
        subtract(StructuredImplicitStencil::J_MINUS,ci,cj-1,ck);
        subtract(StructuredImplicitStencil::J_PLUS,ci,cj+1,ck);
        subtract(StructuredImplicitStencil::K_MINUS,ci,cj,ck-1);
        subtract(StructuredImplicitStencil::K_PLUS,ci,cj,ck+1);
    }
    const SectorGalerkinCoarseSolve::Result solved=
        SectorGalerkinCoarseSolve::solve(
            coarse,aggregated.stencil,pseudoCfl,rightHandSide,40,1.0e-12);
    double maximumSolveError=0.0;
    for(int ci=1;ci<coarse.Nx()-1;++ci)
    for(int cj=1;cj<coarse.Ny()-1;++cj)
    for(int ck=1;ck<coarse.Nz()-1;++ck)
    {
        const int cell=index(coarse,ci,cj,ck);
        maximumSolveError=std::max(maximumSolveError,
            std::abs(solved.correction[static_cast<std::size_t>(cell)].rho
                -coarseCorrection[static_cast<std::size_t>(cell)]));
    }
    if(solved.residualRatio>1.0e-11||maximumSolveError>1.0e-10)
    {
        std::cerr << "Galerkin coarse solve ratio/error="
                  << solved.residualRatio << "/"
                  << maximumSolveError << "\n";
        return 1;
    }

    std::cout << "rap_error=" << maximumRelativeError
              << " diagonal=" << aggregated.minimumActualDiagonal
              << "/" << aggregated.maximumActualDiagonal
              << " dominance=" << aggregated.minimumDiagonalDominance
              << "\ncoarse_solve_sweeps/ratio/error="
              << solved.sweeps << "/" << solved.residualRatio
              << "/" << maximumSolveError
              << "\nSECTOR_GALERKIN_COARSE_STENCIL_PASS\n";
    return 0;
}
