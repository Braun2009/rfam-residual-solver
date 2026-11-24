#ifndef SectorStructuredMomentumInverse_H
#define SectorStructuredMomentumInverse_H

/*---------------------------------------------------------------------------*\
    SectorStructuredMomentumInverse - scalar-stencil inverse for Cartesian
    momentum vectors on the mapped rotational sector.

    The positive operator is the same scalar approximation used by steady
    LU-SGS,

        M x = (1 + 1/CFL) D x - sum_f C_f R_f x_f,

    applied independently to the three momentum components. R_f is identity
    on ordinary faces and the authoritative +/- sector-angle Cartesian
    rotation on wrapped azimuthal faces. Physical-boundary corrections are
    zero. The solve is a residual-correction symmetric Gauss-Seidel iteration;
    it publishes no Domain state and fails closed on invalid stencils.
\*---------------------------------------------------------------------------*/

#include "Operator.H"
#include "SectorMeshMetrics.H"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Foam
{

class SectorStructuredMomentumInverse
{
public:
    using Vector=std::array<double,3>;

    struct Result
    {
        std::vector<Vector> correction;
        int sweeps=0;
        double initialResidualNorm=0.0;
        double finalResidualNorm=0.0;
        double residualRatio=1.0;
        double minimumDiagonal=0.0;
        double maximumDiagonal=0.0;
        double maximumFaceCoupling=0.0;
        double minimumDiagonalDominance=0.0;
    };

private:
    static int index(const Domain& domain,int i,int j,int k)
    { return k*domain.Nx()*domain.Ny()+j*domain.Nx()+i; }

    static bool active(
        const Domain& domain,const std::vector<unsigned char>& activeMask,
        int cell)
    {
        return cell>=0&&cell<domain.N()
            &&activeMask[static_cast<std::size_t>(cell)]!=0U;
    }

    static int wrappedNeighbor(
        const Domain& domain,int i,int j,int k,int face)
    {
        switch(face)
        {
            case StructuredImplicitStencil::I_MINUS:
                return index(domain,i-1,j,k);
            case StructuredImplicitStencil::I_PLUS:
                return index(domain,i+1,j,k);
            case StructuredImplicitStencil::J_MINUS:
                return index(domain,i,j-1,k);
            case StructuredImplicitStencil::J_PLUS:
                return index(domain,i,j+1,k);
            case StructuredImplicitStencil::K_MINUS:
                return index(domain,i,j,k==1?domain.Nz()-2:k-1);
            case StructuredImplicitStencil::K_PLUS:
                return index(domain,i,j,k==domain.Nz()-2?1:k+1);
            default:return -1;
        }
    }

    static Vector vectorInCellFrame(
        const Domain& domain,const SectorMeshMetrics& mesh,
        int k,int face,const Vector& value)
    {
        double angle=0.0;
        if(face==StructuredImplicitStencil::K_MINUS&&k==1)
            angle=-mesh.sectorAngle();
        else if(face==StructuredImplicitStencil::K_PLUS
            &&k==domain.Nz()-2)
            angle=mesh.sectorAngle();
        if(angle==0.0) return value;
        const double cosine=std::cos(angle);
        const double sine=std::sin(angle);
        return {value[0],
            cosine*value[1]-sine*value[2],
            sine*value[1]+cosine*value[2]};
    }

    static void addScaled(Vector& target,const Vector& source,double scale)
    {
        for(int component=0;component<3;++component)
            target[static_cast<std::size_t>(component)]+=scale
                *source[static_cast<std::size_t>(component)];
    }

    static std::vector<Vector> apply(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const StructuredImplicitStencil& stencil,
        const std::vector<double>& diagonal,
        const std::vector<Vector>& field)
    {
        std::vector<Vector> result(
            static_cast<std::size_t>(domain.N()),Vector{0.0,0.0,0.0});
        const int nx=domain.Nx(),ny=domain.Ny();
        const int plane=nx*ny;
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,activeMask,cell)) continue;
            const int i=cell%nx;
            const int j=(cell/nx)%ny;
            const int k=cell/plane;
            Vector value=field[static_cast<std::size_t>(cell)];
            for(double& component:value)
                component*=diagonal[static_cast<std::size_t>(cell)];
            for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
            {
                const int neighbor=wrappedNeighbor(domain,i,j,k,face);
                if(!active(domain,activeMask,neighbor)) continue;
                const Vector neighborValue=vectorInCellFrame(
                    domain,mesh,k,face,
                    field[static_cast<std::size_t>(neighbor)]);
                addScaled(value,neighborValue,
                    -stencil.face[face][static_cast<std::size_t>(cell)]);
            }
            result[static_cast<std::size_t>(cell)]=value;
        }
        return result;
    }

    static double norm(
        const Domain& domain,const std::vector<unsigned char>& activeMask,
        const std::vector<Vector>& field)
    {
        double sum=0.0;
        std::size_t count=0;
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,activeMask,cell)) continue;
            for(double component:field[static_cast<std::size_t>(cell)])
            {
                sum+=component*component;
                ++count;
            }
        }
        return std::sqrt(sum/std::max<std::size_t>(count,1U));
    }

public:
    static std::vector<Vector> applyOperator(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const StructuredImplicitStencil& stencil,double pseudoCfl,
        const std::vector<Vector>& field)
    {
        if(activeMask.size()!=static_cast<std::size_t>(domain.N())
            ||field.size()!=static_cast<std::size_t>(domain.N())
            ||stencil.diagonal.size()!=static_cast<std::size_t>(domain.N())
            ||!std::isfinite(pseudoCfl)||pseudoCfl<=0.0)
            throw std::runtime_error(
                "structured momentum inverse inputs do not match Domain");
        for(const auto& face:stencil.face)
            if(face.size()!=static_cast<std::size_t>(domain.N()))
                throw std::runtime_error(
                    "structured momentum inverse face dimensions do not match");
        const double pseudoFactor=1.0+1.0/pseudoCfl;
        std::vector<double> diagonal(
            static_cast<std::size_t>(domain.N()),0.0);
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,activeMask,cell)) continue;
            const double value=pseudoFactor*stencil.diagonal[
                static_cast<std::size_t>(cell)];
            if(!std::isfinite(value)||value<=0.0)
                throw std::runtime_error(
                    "structured momentum inverse diagonal is not positive");
            diagonal[static_cast<std::size_t>(cell)]=value;
        }
        return apply(domain,mesh,activeMask,stencil,diagonal,field);
    }

    static Result solve(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const StructuredImplicitStencil& stencil,double pseudoCfl,
        const std::vector<Vector>& rightHandSide,
        int maximumSweeps=80,double relativeTolerance=1.0e-8,
        int minimumSweeps=0)
    {
        if(activeMask.size()!=static_cast<std::size_t>(domain.N())
            ||rightHandSide.size()!=static_cast<std::size_t>(domain.N())
            ||stencil.diagonal.size()!=static_cast<std::size_t>(domain.N())
            ||!std::isfinite(pseudoCfl)||pseudoCfl<=0.0
            ||maximumSweeps<=0||minimumSweeps<0
            ||minimumSweeps>maximumSweeps
            ||!std::isfinite(relativeTolerance)
            ||relativeTolerance<=0.0||relativeTolerance>=1.0)
            throw std::runtime_error(
                "invalid structured momentum inverse controls");
        for(const auto& face:stencil.face)
            if(face.size()!=static_cast<std::size_t>(domain.N()))
                throw std::runtime_error(
                    "structured momentum inverse face dimensions do not match");

        const double pseudoFactor=1.0+1.0/pseudoCfl;
        std::vector<double> diagonal(
            static_cast<std::size_t>(domain.N()),0.0);
        Result result;
        result.minimumDiagonal=std::numeric_limits<double>::infinity();
        result.minimumDiagonalDominance=
            std::numeric_limits<double>::infinity();
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,activeMask,cell)) continue;
            const double value=pseudoFactor*stencil.diagonal[
                static_cast<std::size_t>(cell)];
            if(!std::isfinite(value)||value<=0.0)
                throw std::runtime_error(
                    "structured momentum inverse diagonal is not positive");
            diagonal[static_cast<std::size_t>(cell)]=value;
            result.minimumDiagonal=std::min(result.minimumDiagonal,value);
            result.maximumDiagonal=std::max(result.maximumDiagonal,value);
            double couplingSum=0.0;
            for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
            {
                const double coupling=stencil.face[face][
                    static_cast<std::size_t>(cell)];
                if(!std::isfinite(coupling)||coupling<0.0)
                    throw std::runtime_error(
                        "structured momentum inverse coupling is invalid");
                couplingSum+=coupling;
                result.maximumFaceCoupling=std::max(
                    result.maximumFaceCoupling,coupling);
            }
            result.minimumDiagonalDominance=std::min(
                result.minimumDiagonalDominance,value-couplingSum);
        }
        if(!std::isfinite(result.minimumDiagonal)
            ||!std::isfinite(result.minimumDiagonalDominance)
            ||result.minimumDiagonalDominance<=0.0)
            throw std::runtime_error(
                "structured momentum inverse is not strictly dominant");

        result.correction.assign(
            static_cast<std::size_t>(domain.N()),Vector{0.0,0.0,0.0});
        std::vector<Vector> currentResidual=rightHandSide;
        result.initialResidualNorm=norm(domain,activeMask,currentResidual);
        if(!std::isfinite(result.initialResidualNorm))
            throw std::runtime_error(
                "structured momentum inverse initial residual is not finite");
        if(result.initialResidualNorm==0.0)
        {
            result.finalResidualNorm=0.0;
            result.residualRatio=0.0;
            return result;
        }

        const int nx=domain.Nx(),ny=domain.Ny();
        const int plane=nx*ny;
        std::vector<Vector> increment(
            static_cast<std::size_t>(domain.N()),Vector{0.0,0.0,0.0});
        for(int sweep=0;sweep<maximumSweeps;++sweep)
        {
            // Symmetric Gauss-Seidel solve for M delta = residual.
            for(int pass=0;pass<2;++pass)
            {
                const int begin=pass==0?0:domain.N()-1;
                const int end=pass==0?domain.N():-1;
                const int step=pass==0?1:-1;
                for(int cell=begin;cell!=end;cell+=step)
                {
                    if(!active(domain,activeMask,cell)) continue;
                    const int i=cell%nx;
                    const int j=(cell/nx)%ny;
                    const int k=cell/plane;
                    Vector value=currentResidual[
                        static_cast<std::size_t>(cell)];
                    for(int face=0;
                        face<StructuredImplicitStencil::FACE_COUNT;++face)
                    {
                        const int neighbor=wrappedNeighbor(
                            domain,i,j,k,face);
                        if(!active(domain,activeMask,neighbor)) continue;
                        const Vector neighborIncrement=vectorInCellFrame(
                            domain,mesh,k,face,
                            increment[static_cast<std::size_t>(neighbor)]);
                        addScaled(value,neighborIncrement,
                            stencil.face[face][
                                static_cast<std::size_t>(cell)]);
                    }
                    for(double& component:value)
                        component/=diagonal[
                            static_cast<std::size_t>(cell)];
                    increment[static_cast<std::size_t>(cell)]=value;
                }
            }
            for(int cell=0;cell<domain.N();++cell)
                if(active(domain,activeMask,cell))
                    addScaled(result.correction[
                        static_cast<std::size_t>(cell)],
                        increment[static_cast<std::size_t>(cell)],1.0);

            const std::vector<Vector> action=apply(
                domain,mesh,activeMask,stencil,diagonal,
                result.correction);
            currentResidual=rightHandSide;
            for(int cell=0;cell<domain.N();++cell)
                if(active(domain,activeMask,cell))
                    addScaled(currentResidual[
                        static_cast<std::size_t>(cell)],
                        action[static_cast<std::size_t>(cell)],-1.0);
            result.finalResidualNorm=norm(
                domain,activeMask,currentResidual);
            result.residualRatio=result.finalResidualNorm
                /result.initialResidualNorm;
            result.sweeps=sweep+1;
            if(!std::isfinite(result.finalResidualNorm))
                throw std::runtime_error(
                    "structured momentum inverse residual is not finite");
            if(result.sweeps>=minimumSweeps
                &&result.residualRatio<=relativeTolerance) break;
            std::fill(increment.begin(),increment.end(),
                Vector{0.0,0.0,0.0});
        }
        return result;
    }
};

} // End namespace Foam

#endif
