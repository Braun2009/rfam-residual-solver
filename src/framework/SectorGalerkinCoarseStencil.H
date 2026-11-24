#ifndef SectorGalerkinCoarseStencil_H
#define SectorGalerkinCoarseStencil_H

/*---------------------------------------------------------------------------*\
    SectorGalerkinCoarseStencil - volume-weighted RAP aggregation of the
    positive scalar structured LU-SGS approximation.

    Fine corrections are prolonged piecewise constantly and residual rates are
    restricted by cell volume. The returned operator stencil is scaled so that
    SteadyLUSGSIntegrator, with the same pseudo CFL, reconstructs

        A_H = R A_h P

    including fine couplings internal to one coarse aggregate. Halo correction
    values are zero, matching the steady integrator's active-unknown contract.
\*---------------------------------------------------------------------------*/

#include "Operator.H"
#include "SectorConservativeGridTransfer.H"
#include "SectorMeshMetrics.H"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Foam
{

class SectorGalerkinCoarseStencil
{
public:
    struct Result
    {
        StructuredImplicitStencil stencil;
        double minimumActualDiagonal=std::numeric_limits<double>::infinity();
        double maximumActualDiagonal=0.0;
        double maximumFaceCoupling=0.0;
        double minimumDiagonalDominance=std::numeric_limits<double>::infinity();
        int activeCoarseCells=0;
    };

private:
    static int index(const Domain& domain,int i,int j,int k)
    { return k*domain.Nx()*domain.Ny()+j*domain.Nx()+i; }

    static bool active(const Domain& domain,int i,int j,int k)
    {
        return i>0&&i<domain.Nx()-1
            &&j>0&&j<domain.Ny()-1
            &&k>0&&k<domain.Nz()-1;
    }

public:
    static Result aggregate(
        const Domain& fine,
        const SectorMeshMetrics& fineMesh,
        const StructuredImplicitStencil& fineOperatorStencil,
        const Domain& coarse,
        const SectorMeshMetrics& coarseMesh,
        double pseudoCfl)
    {
        if(!std::isfinite(pseudoCfl)||pseudoCfl<=0.0)
            throw std::invalid_argument(
                "Galerkin coarse stencil pseudo CFL must be positive");
        if(fineOperatorStencil.diagonal.size()
                !=static_cast<std::size_t>(fine.N()))
            throw std::runtime_error(
                "fine operator stencil dimensions do not match Domain");
        for(const auto& face:fineOperatorStencil.face)
            if(face.size()!=static_cast<std::size_t>(fine.N()))
                throw std::runtime_error(
                    "fine operator face stencil dimensions do not match Domain");
        const SectorConservativeGridTransfer::Plan transfer=
            SectorConservativeGridTransfer::plan(fine,coarse);
        if(transfer.direction
            !=SectorConservativeGridTransfer::Direction::RESTRICT)
            throw std::runtime_error(
                "Galerkin aggregation requires an aligned fine-to-coarse grid pair");

        const double pseudoFactor=1.0+1.0/pseudoCfl;
        Result result;
        result.stencil.ensureSize(static_cast<std::size_t>(coarse.N()));
        for(int ci=1;ci<coarse.Nx()-1;++ci)
        for(int cj=1;cj<coarse.Ny()-1;++cj)
        for(int ck=1;ck<coarse.Nz()-1;++ck)
        {
            const int coarseCell=index(coarse,ci,cj,ck);
            const double coarseVolume=coarseMesh.volume(ci,cj,ck);
            if(!std::isfinite(coarseVolume)||coarseVolume<=0.0)
                throw std::runtime_error(
                    "Galerkin coarse cell volume is invalid");
            double actualDiagonal=0.0;
            double faceCoupling[StructuredImplicitStencil::FACE_COUNT]={};
            const int beginI=1+(ci-1)*transfer.ratioI;
            const int beginJ=1+(cj-1)*transfer.ratioJ;
            const int beginK=1+(ck-1)*transfer.ratioK;
            for(int fi=beginI;fi<beginI+transfer.ratioI;++fi)
            for(int fj=beginJ;fj<beginJ+transfer.ratioJ;++fj)
            for(int fk=beginK;fk<beginK+transfer.ratioK;++fk)
            {
                if(!active(fine,fi,fj,fk))
                    throw std::runtime_error(
                        "Galerkin aggregate contains a non-active fine child");
                const int fineCell=index(fine,fi,fj,fk);
                const double weight=fineMesh.volume(fi,fj,fk)/coarseVolume;
                const double rawDiagonal=fineOperatorStencil.diagonal[
                    static_cast<std::size_t>(fineCell)];
                if(!std::isfinite(rawDiagonal)||rawDiagonal<0.0
                    ||!std::isfinite(weight)||weight<=0.0)
                    throw std::runtime_error(
                        "Galerkin fine diagonal or volume weight is invalid");
                actualDiagonal+=weight*pseudoFactor*rawDiagonal;

                auto aggregateFace=[&](
                    int face,int ni,int nj,int nk,
                    int neighborCi,int neighborCj,int neighborCk)
                {
                    const double coupling=fineOperatorStencil.face[face][
                        static_cast<std::size_t>(fineCell)];
                    if(!std::isfinite(coupling)||coupling<0.0)
                        throw std::runtime_error(
                            "Galerkin fine face coupling is invalid");
                    if(coupling==0.0) return;
                    if(neighborCi==ci&&neighborCj==cj&&neighborCk==ck)
                        actualDiagonal-=weight*coupling;
                    else
                    {
                        // A structured aggregate may couple only to the same
                        // coarse cell or the corresponding adjacent coarse
                        // face, including a zero-correction halo aggregate.
                        const int distance=std::abs(neighborCi-ci)
                            +std::abs(neighborCj-cj)
                            +std::abs(neighborCk-ck);
                        if(distance!=1||ni<0||ni>=fine.Nx()
                            ||nj<0||nj>=fine.Ny()
                            ||nk<0||nk>=fine.Nz())
                            throw std::runtime_error(
                                "Galerkin fine face does not map to an adjacent coarse aggregate");
                        faceCoupling[face]+=weight*coupling;
                    }
                };
                const auto coarseI=[&](int i)
                { return i<=0?0:(i>=fine.Nx()-1?coarse.Nx()-1
                    :1+(i-1)/transfer.ratioI); };
                const auto coarseJ=[&](int j)
                { return j<=0?0:(j>=fine.Ny()-1?coarse.Ny()-1
                    :1+(j-1)/transfer.ratioJ); };
                const auto coarseK=[&](int k)
                { return k<=0?0:(k>=fine.Nz()-1?coarse.Nz()-1
                    :1+(k-1)/transfer.ratioK); };
                aggregateFace(StructuredImplicitStencil::I_MINUS,
                    fi-1,fj,fk,coarseI(fi-1),cj,ck);
                aggregateFace(StructuredImplicitStencil::I_PLUS,
                    fi+1,fj,fk,coarseI(fi+1),cj,ck);
                aggregateFace(StructuredImplicitStencil::J_MINUS,
                    fi,fj-1,fk,ci,coarseJ(fj-1),ck);
                aggregateFace(StructuredImplicitStencil::J_PLUS,
                    fi,fj+1,fk,ci,coarseJ(fj+1),ck);
                aggregateFace(StructuredImplicitStencil::K_MINUS,
                    fi,fj,fk-1,ci,cj,coarseK(fk-1));
                aggregateFace(StructuredImplicitStencil::K_PLUS,
                    fi,fj,fk+1,ci,cj,coarseK(fk+1));
            }

            double couplingSum=0.0;
            for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
            {
                const double coupling=faceCoupling[face];
                result.stencil.face[face][
                    static_cast<std::size_t>(coarseCell)]=coupling;
                couplingSum+=coupling;
                result.maximumFaceCoupling=std::max(
                    result.maximumFaceCoupling,coupling);
            }
            if(!std::isfinite(actualDiagonal)||actualDiagonal<=0.0)
                throw std::runtime_error(
                    "Galerkin coarse actual diagonal is not positive");
            result.stencil.diagonal[
                static_cast<std::size_t>(coarseCell)]=
                actualDiagonal/pseudoFactor;
            result.minimumActualDiagonal=std::min(
                result.minimumActualDiagonal,actualDiagonal);
            result.maximumActualDiagonal=std::max(
                result.maximumActualDiagonal,actualDiagonal);
            result.minimumDiagonalDominance=std::min(
                result.minimumDiagonalDominance,
                actualDiagonal-couplingSum);
            ++result.activeCoarseCells;
        }
        if(result.activeCoarseCells==0)
            throw std::runtime_error(
                "Galerkin aggregation produced no active coarse cells");
        return result;
    }
};

} // End namespace Foam

#endif
