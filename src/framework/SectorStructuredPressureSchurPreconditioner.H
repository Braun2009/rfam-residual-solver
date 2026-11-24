#ifndef SectorStructuredPressureSchurPreconditioner_H
#define SectorStructuredPressureSchurPreconditioner_H

/*---------------------------------------------------------------------------*\
    SectorStructuredPressureSchurPreconditioner - matrix-free pressure Schur
    solve using the same scalar structured operator as steady LU-SGS.

      M_rho delta_rho + B delta_m = F_rho
      M_m   delta_m   + G delta_p = F_m
      delta_p = a^2 delta_rho

    gives

      [M_rho/a^2 - B M_m^-1 G + D_AUSM] delta_p
          = F_rho - B M_m^-1 F_m.

    D_AUSM is the exact frozen-base pressure-diffusion derivative of the
    production AUSM+-up mass flux, preassembled once per nonlinear state.
    M_rho and M_m use the positive D-C_face stencil and pseudo factor
    (1+1/CFL). M_m^-1 is the fixed-sweep rotationally periodic structured
    momentum inverse. G is the mapped least-squares pressure gradient; B is
    the conservative face divergence with controlled inlet mass flux, solid
    axis/wall flux, extrapolating outlet flux, and rotational periodicity.

    The pressure system is solved by matrix-free GMRES. No nonlinear state is
    published here. Failure to reach the complete linear residual tolerance is
    fail-close and produces no correction direction.
\*---------------------------------------------------------------------------*/

#include "SectorPressureSchurPreconditioner.H"
#include "SectorStructuredMomentumInverse.H"
#include "AUSMPlusUpFlux.H"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Foam
{

class SectorStructuredPressureSchurPreconditioner
{
public:
    using Vector=SectorStructuredMomentumInverse::Vector;

    struct Audit
    {
        int momentumSolves=0;
        int momentumSweeps=0;
        double maximumMomentumResidualRatio=0.0;
        int acousticPreconditionerSolves=0;
        int acousticPreconditionerSweeps=0;
        double maximumAcousticPreconditionerResidualRatio=0.0;
        int leastSquaresApplications=0;
        int leastSquaresFallbacks=0;
        int pressureDiffusionFaces=0;
        double minimumPressureDiffusionMobility=
            std::numeric_limits<double>::infinity();
        double maximumPressureDiffusionMobility=0.0;
        double maximumPressureDiffusionCoefficient=0.0;
        double maximumPressureDiffusionReciprocityError=0.0;
    };

    struct Result
    {
        std::vector<double> pressure;
        std::vector<Vector> velocity;
        int krylovIterations=0;
        double initialResidualNorm=0.0;
        double finalResidualNorm=0.0;
        double residualRatio=1.0;
        double maximumRelativePressure=0.0;
        double maximumVelocityMach=0.0;
        Audit audit;
    };

private:
    struct PressureDiffusionStencil
    {
        std::vector<double> face[StructuredImplicitStencil::FACE_COUNT];
        Audit audit;
    };

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

    static FaceAreaVector faceArea(
        const SectorMeshMetrics& mesh,int i,int j,int k,int face)
    {
        switch(face)
        {
            case StructuredImplicitStencil::I_MINUS:
                return mesh.axialAreaVector(i,j,k,-1);
            case StructuredImplicitStencil::I_PLUS:
                return mesh.axialAreaVector(i+1,j,k,+1);
            case StructuredImplicitStencil::J_MINUS:
                return mesh.radialAreaVector(i,j,k,-1);
            case StructuredImplicitStencil::J_PLUS:
                return mesh.radialAreaVector(i,j+1,k,+1);
            case StructuredImplicitStencil::K_MINUS:
                return mesh.azimuthalAreaVector(i,j,k,-1);
            case StructuredImplicitStencil::K_PLUS:
                return mesh.azimuthalAreaVector(i,j,k+1,+1);
            default:return FaceAreaVector{0.0,0.0,0.0};
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

    static double dot(const Vector& lhs,const FaceAreaVector& rhs)
    { return lhs[0]*rhs.x+lhs[1]*rhs.y+lhs[2]*rhs.z; }

    static std::vector<double> divergence(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const std::vector<Vector>& field)
    {
        std::vector<double> result(
            static_cast<std::size_t>(domain.N()),0.0);
        const int nx=domain.Nx(),ny=domain.Ny();
        const int plane=nx*ny;
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,activeMask,cell)) continue;
            const int i=cell%nx;
            const int j=(cell/nx)%ny;
            const int k=cell/plane;
            double flux=0.0;
            for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
            {
                const int neighbor=wrappedNeighbor(domain,i,j,k,face);
                const bool neighborActive=active(
                    domain,activeMask,neighbor);
                // The controlled injection, axis and wall have zero correction
                // mass flux. The steady outlet extrapolates the cell value.
                if(!neighborActive
                    &&face!=StructuredImplicitStencil::I_PLUS) continue;
                Vector faceValue=field[static_cast<std::size_t>(cell)];
                if(neighborActive)
                {
                    const Vector neighborValue=vectorInCellFrame(
                        domain,mesh,k,face,
                        field[static_cast<std::size_t>(neighbor)]);
                    for(int component=0;component<3;++component)
                        faceValue[static_cast<std::size_t>(component)]=0.5*(
                            faceValue[static_cast<std::size_t>(component)]
                            +neighborValue[
                                static_cast<std::size_t>(component)]);
                }
                flux+=dot(faceValue,faceArea(mesh,i,j,k,face));
            }
            result[static_cast<std::size_t>(cell)]=
                flux/mesh.volume(i,j,k);
        }
        return result;
    }

    static PrimitiveVars neighborPrimitiveInCellFrame(
        const Domain& domain,const SectorMeshMetrics& mesh,
        int k,int face,const PrimitiveVars& primitive)
    {
        PrimitiveVars result=primitive;
        double angle=0.0;
        if(face==StructuredImplicitStencil::K_MINUS&&k==1)
            angle=-mesh.sectorAngle();
        else if(face==StructuredImplicitStencil::K_PLUS
            &&k==domain.Nz()-2)
            angle=mesh.sectorAngle();
        if(angle!=0.0)
        {
            const double cosine=std::cos(angle);
            const double sine=std::sin(angle);
            const double radial=primitive.U[1];
            const double tangential=primitive.U[2];
            result.U[1]=cosine*radial-sine*tangential;
            result.U[2]=sine*radial+cosine*tangential;
        }
        return result;
    }

    static std::vector<double> applyScalarOperator(
        const Domain& domain,const std::vector<unsigned char>& activeMask,
        const StructuredImplicitStencil& stencil,double pseudoCfl,
        const std::vector<double>& field)
    {
        const double pseudoFactor=1.0+1.0/pseudoCfl;
        std::vector<double> result(
            static_cast<std::size_t>(domain.N()),0.0);
        const int nx=domain.Nx(),ny=domain.Ny();
        const int plane=nx*ny;
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,activeMask,cell)) continue;
            const int i=cell%nx;
            const int j=(cell/nx)%ny;
            const int k=cell/plane;
            double value=pseudoFactor*stencil.diagonal[
                static_cast<std::size_t>(cell)]
                *field[static_cast<std::size_t>(cell)];
            for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
            {
                const int neighbor=wrappedNeighbor(domain,i,j,k,face);
                if(!active(domain,activeMask,neighbor)) continue;
                value-=stencil.face[face][static_cast<std::size_t>(cell)]
                    *field[static_cast<std::size_t>(neighbor)];
            }
            result[static_cast<std::size_t>(cell)]=value;
        }
        return result;
    }

    static PressureDiffusionStencil assembleAUSMPressureDiffusion(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask)
    {
        PressureDiffusionStencil stencil;
        for(auto& face:stencil.face)
            face.assign(static_cast<std::size_t>(domain.N()),0.0);
        const int nx=domain.Nx(),ny=domain.Ny();
        const int plane=nx*ny;
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,activeMask,cell)) continue;
            const int i=cell%nx;
            const int j=(cell/nx)%ny;
            const int k=cell/plane;
            const PrimitiveVars& left=domain.W(cell);
            const thermo::EquilibriumState& leftQ=domain.Q(cell);
            for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
            {
                const int neighbor=wrappedNeighbor(domain,i,j,k,face);
                if(!active(domain,activeMask,neighbor)) continue;
                const PrimitiveVars right=neighborPrimitiveInCellFrame(
                    domain,mesh,k,face,domain.W(neighbor));
                const thermo::EquilibriumState& rightQ=domain.Q(neighbor);
                const FaceAreaVector area=faceArea(mesh,i,j,k,face);
                const double areaMagnitude=std::sqrt(
                    area.x*area.x+area.y*area.y+area.z*area.z);
                if(!std::isfinite(areaMagnitude)||areaMagnitude<=1.0e-30)
                    throw std::runtime_error(
                        "structured pressure Schur AUSM face area is invalid");
                const std::vector<double> normal={
                    area.x/areaMagnitude,area.y/areaMagnitude,
                    area.z/areaMagnitude};
                const double mobility=ausmPressureDiffusionMassMobility(
                    left,right,normal,leftQ.rho,
                    leftQ.acoustic.soundSpeed(),rightQ.rho,
                    rightQ.acoustic.soundSpeed());
                if(!std::isfinite(mobility)||mobility<0.0)
                    throw std::runtime_error(
                        "structured pressure Schur AUSM pressure mobility is invalid");
                const double coefficient=mobility*areaMagnitude
                    /mesh.volume(i,j,k);
                if(!std::isfinite(coefficient)||coefficient<0.0)
                    throw std::runtime_error(
                        "structured pressure Schur AUSM pressure coefficient is invalid");
                stencil.face[face][static_cast<std::size_t>(cell)]=
                    coefficient;
                ++stencil.audit.pressureDiffusionFaces;
                stencil.audit.minimumPressureDiffusionMobility=std::min(
                    stencil.audit.minimumPressureDiffusionMobility,mobility);
                stencil.audit.maximumPressureDiffusionMobility=std::max(
                    stencil.audit.maximumPressureDiffusionMobility,mobility);
                stencil.audit.maximumPressureDiffusionCoefficient=std::max(
                    stencil.audit.maximumPressureDiffusionCoefficient,
                    coefficient);
            }
        }
        const int reciprocalFace[StructuredImplicitStencil::FACE_COUNT]={
            StructuredImplicitStencil::I_PLUS,
            StructuredImplicitStencil::I_MINUS,
            StructuredImplicitStencil::J_PLUS,
            StructuredImplicitStencil::J_MINUS,
            StructuredImplicitStencil::K_PLUS,
            StructuredImplicitStencil::K_MINUS};
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,activeMask,cell)) continue;
            const int i=cell%nx;
            const int j=(cell/nx)%ny;
            const int k=cell/plane;
            for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
            {
                const int neighbor=wrappedNeighbor(domain,i,j,k,face);
                if(!active(domain,activeMask,neighbor)) continue;
                const int ni=neighbor%nx;
                const int nj=(neighbor/nx)%ny;
                const int nk=neighbor/plane;
                const double left=stencil.face[face][
                    static_cast<std::size_t>(cell)]*mesh.volume(i,j,k);
                const double right=stencil.face[reciprocalFace[face]][
                    static_cast<std::size_t>(neighbor)]
                    *mesh.volume(ni,nj,nk);
                stencil.audit.maximumPressureDiffusionReciprocityError=
                    std::max(
                        stencil.audit.maximumPressureDiffusionReciprocityError,
                        std::abs(left-right)/std::max(
                            {std::abs(left),std::abs(right),1.0e-300}));
            }
        }
        if(!std::isfinite(stencil.audit.minimumPressureDiffusionMobility))
            stencil.audit.minimumPressureDiffusionMobility=0.0;
        return stencil;
    }

    static void addAUSMPressureDiffusionAction(
        const Domain& domain,const std::vector<unsigned char>& activeMask,
        const PressureDiffusionStencil& stencil,
        const std::vector<double>& pressure,
        std::vector<double>& action)
    {
        const int nx=domain.Nx(),ny=domain.Ny();
        const int plane=nx*ny;
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,activeMask,cell)) continue;
            const int i=cell%nx;
            const int j=(cell/nx)%ny;
            const int k=cell/plane;
            for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
            {
                const int neighbor=wrappedNeighbor(domain,i,j,k,face);
                if(!active(domain,activeMask,neighbor)) continue;
                const double coefficient=stencil.face[face][
                    static_cast<std::size_t>(cell)];
                action[static_cast<std::size_t>(cell)]+=coefficient
                    *(pressure[static_cast<std::size_t>(cell)]
                        -pressure[static_cast<std::size_t>(neighbor)]);
            }
        }
    }

    static std::vector<Vector> pressureGradient(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const std::vector<unsigned char>& inletPressureDirichlet,
        const std::vector<double>& pressure,Audit& audit)
    {
        std::vector<Vector> result(
            static_cast<std::size_t>(domain.N()),Vector{0.0,0.0,0.0});
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,activeMask,cell)) continue;
            bool applied=false;
            const FaceAreaVector gradient=
                SectorPressureSchurPreconditioner::reconstructPressureGradient(
                    domain,mesh,activeMask,inletPressureDirichlet,
                    pressure,cell,true,&applied);
            if(applied) ++audit.leastSquaresApplications;
            else ++audit.leastSquaresFallbacks;
            result[static_cast<std::size_t>(cell)]={
                gradient.x,gradient.y,gradient.z};
        }
        return result;
    }

    static std::vector<Vector> momentumInverse(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const StructuredImplicitStencil& stencil,double pseudoCfl,
        const std::vector<Vector>& rhs,Audit& audit)
    {
        constexpr int fixedSweeps=24;
        const auto solved=SectorStructuredMomentumInverse::solve(
            domain,mesh,activeMask,stencil,pseudoCfl,rhs,
            fixedSweeps,1.0e-14,fixedSweeps);
        ++audit.momentumSolves;
        audit.momentumSweeps+=solved.sweeps;
        audit.maximumMomentumResidualRatio=std::max(
            audit.maximumMomentumResidualRatio,solved.residualRatio);
        return solved.correction;
    }

    static std::vector<double> acousticPreconditioner(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const StructuredImplicitStencil& stencil,
        const std::vector<unsigned char>& inletPressureDirichlet,
        double pseudoCfl,const std::vector<double>& rhs,Audit& audit)
    {
        std::vector<ConservedVars> residual(
            static_cast<std::size_t>(domain.N()),
            ConservedVars(domain.nSpecies()));
        for(int cell=0;cell<domain.N();++cell)
            if(active(domain,activeMask,cell))
                residual[static_cast<std::size_t>(cell)].rho=
                    rhs[static_cast<std::size_t>(cell)];
        constexpr int fixedSweeps=24;
        // M = D*(1+1/CFL)-C, so the diagonal momentum mobility is
        // 1/[D*(1+1/CFL)]. SectorPressureSchurPreconditioner represents
        // mobility as CFL_eff/D; therefore CFL_eff=CFL/(CFL+1).
        const double effectiveCfl=pseudoCfl/(pseudoCfl+1.0);
        const auto solved=SectorPressureSchurPreconditioner::build(
            domain,mesh,residual,activeMask,stencil,
            inletPressureDirichlet,effectiveCfl,
            fixedSweeps,1.0e-14,false,0,true,fixedSweeps);
        ++audit.acousticPreconditionerSolves;
        audit.acousticPreconditionerSweeps+=solved.sweeps;
        audit.maximumAcousticPreconditionerResidualRatio=std::max(
            audit.maximumAcousticPreconditionerResidualRatio,
            solved.residualRatio);
        return solved.pressure;
    }

    static void mergeAudit(Audit& target,const Audit& source)
    {
        target.momentumSolves+=source.momentumSolves;
        target.momentumSweeps+=source.momentumSweeps;
        target.maximumMomentumResidualRatio=std::max(
            target.maximumMomentumResidualRatio,
            source.maximumMomentumResidualRatio);
        target.acousticPreconditionerSolves+=
            source.acousticPreconditionerSolves;
        target.acousticPreconditionerSweeps+=
            source.acousticPreconditionerSweeps;
        target.maximumAcousticPreconditionerResidualRatio=std::max(
            target.maximumAcousticPreconditionerResidualRatio,
            source.maximumAcousticPreconditionerResidualRatio);
        target.leastSquaresApplications+=source.leastSquaresApplications;
        target.leastSquaresFallbacks+=source.leastSquaresFallbacks;
        target.pressureDiffusionFaces+=source.pressureDiffusionFaces;
        target.minimumPressureDiffusionMobility=std::min(
            target.minimumPressureDiffusionMobility,
            source.minimumPressureDiffusionMobility);
        target.maximumPressureDiffusionMobility=std::max(
            target.maximumPressureDiffusionMobility,
            source.maximumPressureDiffusionMobility);
        target.maximumPressureDiffusionCoefficient=std::max(
            target.maximumPressureDiffusionCoefficient,
            source.maximumPressureDiffusionCoefficient);
        target.maximumPressureDiffusionReciprocityError=std::max(
            target.maximumPressureDiffusionReciprocityError,
            source.maximumPressureDiffusionReciprocityError);
    }

    static double scalarDot(
        const Domain& domain,const std::vector<unsigned char>& activeMask,
        const std::vector<double>& lhs,const std::vector<double>& rhs)
    {
        double result=0.0;
        for(int cell=0;cell<domain.N();++cell)
            if(active(domain,activeMask,cell))
                result+=lhs[static_cast<std::size_t>(cell)]
                    *rhs[static_cast<std::size_t>(cell)];
        return result;
    }

    static double scalarNorm(
        const Domain& domain,const std::vector<unsigned char>& activeMask,
        const std::vector<double>& field)
    {
        double sum=0.0;
        for(int cell=0;cell<domain.N();++cell)
            if(active(domain,activeMask,cell))
            {
                const double value=field[static_cast<std::size_t>(cell)];
                sum+=value*value;
            }
        return std::sqrt(sum);
    }

    static void addScaled(
        const Domain& domain,const std::vector<unsigned char>& activeMask,
        std::vector<double>& target,const std::vector<double>& source,
        double scale)
    {
        for(int cell=0;cell<domain.N();++cell)
            if(active(domain,activeMask,cell))
                target[static_cast<std::size_t>(cell)]+=scale
                    *source[static_cast<std::size_t>(cell)];
    }

    static void scale(
        const Domain& domain,const std::vector<unsigned char>& activeMask,
        std::vector<double>& field,double factor)
    {
        for(int cell=0;cell<domain.N();++cell)
            if(active(domain,activeMask,cell))
                field[static_cast<std::size_t>(cell)]*=factor;
    }

    static std::vector<double> applyPressureOperatorWithStencil(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const StructuredImplicitStencil& stencil,
        const PressureDiffusionStencil& pressureDiffusion,
        const std::vector<unsigned char>& inletPressureDirichlet,
        double pseudoCfl,const std::vector<double>& pressure,
        Audit& audit)
    {
        std::vector<double> densityCorrection(
            static_cast<std::size_t>(domain.N()),0.0);
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,activeMask,cell)) continue;
            const double sound=domain.Q(cell).acoustic.soundSpeed();
            if(!std::isfinite(sound)||sound<=0.0)
                throw std::runtime_error(
                    "structured pressure Schur acoustics are invalid");
            densityCorrection[static_cast<std::size_t>(cell)]=
                pressure[static_cast<std::size_t>(cell)]/(sound*sound);
        }
        std::vector<double> result=applyScalarOperator(
            domain,activeMask,stencil,pseudoCfl,densityCorrection);
        const std::vector<Vector> gradient=pressureGradient(
            domain,mesh,activeMask,inletPressureDirichlet,pressure,audit);
        const std::vector<Vector> momentumResponse=momentumInverse(
            domain,mesh,activeMask,stencil,pseudoCfl,gradient,audit);
        const std::vector<double> responseDivergence=divergence(
            domain,mesh,activeMask,momentumResponse);
        for(int cell=0;cell<domain.N();++cell)
            if(active(domain,activeMask,cell))
                result[static_cast<std::size_t>(cell)]
                    -=responseDivergence[static_cast<std::size_t>(cell)];
        addAUSMPressureDiffusionAction(
            domain,activeMask,pressureDiffusion,pressure,result);
        return result;
    }

public:
    static std::vector<double> applyPressureOperator(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const StructuredImplicitStencil& stencil,
        const std::vector<unsigned char>& inletPressureDirichlet,
        double pseudoCfl,const std::vector<double>& pressure,
        Audit& audit)
    {
        if(activeMask.size()!=static_cast<std::size_t>(domain.N())
            ||pressure.size()!=static_cast<std::size_t>(domain.N())
            ||stencil.diagonal.size()!=static_cast<std::size_t>(domain.N())
            ||inletPressureDirichlet.size()!=static_cast<std::size_t>(
                domain.Ny()*domain.Nz())
            ||!std::isfinite(pseudoCfl)||pseudoCfl<=0.0)
            throw std::runtime_error(
                "structured pressure Schur inputs do not match Domain");
        const PressureDiffusionStencil pressureDiffusion=
            assembleAUSMPressureDiffusion(domain,mesh,activeMask);
        mergeAudit(audit,pressureDiffusion.audit);
        return applyPressureOperatorWithStencil(
            domain,mesh,activeMask,stencil,pressureDiffusion,
            inletPressureDirichlet,pseudoCfl,pressure,audit);
    }

    static Result build(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<ConservedVars>& residual,
        const std::vector<unsigned char>& activeMask,
        const StructuredImplicitStencil& stencil,
        const std::vector<unsigned char>& inletPressureDirichlet,
        double pseudoCfl,int maximumKrylovIterations=24,
        double relativeTolerance=1.0e-8)
    {
        if(residual.size()!=static_cast<std::size_t>(domain.N())
            ||activeMask.size()!=static_cast<std::size_t>(domain.N())
            ||stencil.diagonal.size()!=static_cast<std::size_t>(domain.N())
            ||maximumKrylovIterations<=0
            ||maximumKrylovIterations>160
            ||!std::isfinite(relativeTolerance)
            ||relativeTolerance<=0.0||relativeTolerance>=1.0)
            throw std::runtime_error(
                "invalid structured pressure Schur controls");
        for(const auto& face:stencil.face)
            if(face.size()!=static_cast<std::size_t>(domain.N()))
                throw std::runtime_error(
                    "structured pressure Schur face dimensions do not match");

        Result result;
        result.pressure.assign(
            static_cast<std::size_t>(domain.N()),0.0);
        result.velocity.assign(
            static_cast<std::size_t>(domain.N()),Vector{0.0,0.0,0.0});
        const PressureDiffusionStencil pressureDiffusion=
            assembleAUSMPressureDiffusion(domain,mesh,activeMask);
        mergeAudit(result.audit,pressureDiffusion.audit);

        std::vector<Vector> momentumResidual(
            static_cast<std::size_t>(domain.N()),Vector{0.0,0.0,0.0});
        for(int cell=0;cell<domain.N();++cell)
            if(active(domain,activeMask,cell))
                momentumResidual[static_cast<std::size_t>(cell)]=
                    {residual[static_cast<std::size_t>(cell)].rhoU[0],
                     residual[static_cast<std::size_t>(cell)].rhoU[1],
                     residual[static_cast<std::size_t>(cell)].rhoU[2]};
        const std::vector<Vector> predictedMomentum=momentumInverse(
            domain,mesh,activeMask,stencil,pseudoCfl,
            momentumResidual,result.audit);
        const std::vector<double> predictedDivergence=divergence(
            domain,mesh,activeMask,predictedMomentum);
        std::vector<double> rightHandSide(
            static_cast<std::size_t>(domain.N()),0.0);
        for(int cell=0;cell<domain.N();++cell)
            if(active(domain,activeMask,cell))
                rightHandSide[static_cast<std::size_t>(cell)]=
                    residual[static_cast<std::size_t>(cell)].rho
                    -predictedDivergence[static_cast<std::size_t>(cell)];

        result.initialResidualNorm=scalarNorm(
            domain,activeMask,rightHandSide);
        if(!std::isfinite(result.initialResidualNorm))
            throw std::runtime_error(
                "structured pressure Schur initial residual is not finite");
        if(result.initialResidualNorm==0.0)
        {
            result.finalResidualNorm=0.0;
            result.residualRatio=0.0;
            return result;
        }

        const int restartDimension=maximumKrylovIterations;
        std::vector<double> trueResidual=rightHandSide;
        result.finalResidualNorm=result.initialResidualNorm;
        result.residualRatio=1.0;
        while(result.krylovIterations<maximumKrylovIterations
            &&result.residualRatio>relativeTolerance)
        {
            Audit initialPreconditionerAudit;
            std::vector<double> first=acousticPreconditioner(
                domain,mesh,activeMask,stencil,
                inletPressureDirichlet,pseudoCfl,
                trueResidual,initialPreconditionerAudit);
            mergeAudit(result.audit,initialPreconditionerAudit);
            const double krylovInitialNorm=scalarNorm(
                domain,activeMask,first);
            if(!std::isfinite(krylovInitialNorm)
                ||krylovInitialNorm<=0.0)
                throw std::runtime_error(
                    "structured pressure Schur preconditioned residual is invalid");
            scale(domain,activeMask,first,1.0/krylovInitialNorm);

            const int dimension=std::min(
                restartDimension,
                maximumKrylovIterations-result.krylovIterations);
            std::vector<std::vector<double>> basis;
            basis.reserve(static_cast<std::size_t>(dimension+1));
            basis.push_back(std::move(first));
            std::vector<std::vector<double>> hessenberg(
                static_cast<std::size_t>(dimension+1),
                std::vector<double>(
                    static_cast<std::size_t>(dimension),0.0));
            std::vector<double> cosine(
                static_cast<std::size_t>(dimension),0.0);
            std::vector<double> sine(
                static_cast<std::size_t>(dimension),0.0);
            std::vector<double> transformed(
                static_cast<std::size_t>(dimension+1),0.0);
            transformed[0]=krylovInitialNorm;
            int used=0;
            for(int column=0;column<dimension;++column)
            {
                Audit actionAudit;
                const std::vector<double> rawValue=
                    applyPressureOperatorWithStencil(
                        domain,mesh,activeMask,stencil,pressureDiffusion,
                        inletPressureDirichlet,pseudoCfl,
                        basis[static_cast<std::size_t>(column)],
                        actionAudit);
                Audit preconditionerAudit;
                std::vector<double> value=acousticPreconditioner(
                    domain,mesh,activeMask,stencil,
                    inletPressureDirichlet,pseudoCfl,
                    rawValue,preconditionerAudit);
                mergeAudit(result.audit,actionAudit);
                mergeAudit(result.audit,preconditionerAudit);
                for(int orthogonalization=0;orthogonalization<2;
                    ++orthogonalization)
                for(int row=0;row<=column;++row)
                {
                    const double projection=scalarDot(
                        domain,activeMask,
                        basis[static_cast<std::size_t>(row)],value);
                    hessenberg[static_cast<std::size_t>(row)]
                        [static_cast<std::size_t>(column)]+=projection;
                    addScaled(domain,activeMask,value,
                        basis[static_cast<std::size_t>(row)],-projection);
                }
                const double nextNorm=scalarNorm(
                    domain,activeMask,value);
                hessenberg[static_cast<std::size_t>(column+1)]
                    [static_cast<std::size_t>(column)]=nextNorm;
                if(nextNorm>1.0e-300)
                {
                    scale(domain,activeMask,value,1.0/nextNorm);
                    basis.push_back(std::move(value));
                }
                else basis.push_back(std::vector<double>(
                    static_cast<std::size_t>(domain.N()),0.0));

                for(int row=0;row<column;++row)
                {
                    const double upper=hessenberg[
                        static_cast<std::size_t>(row)]
                        [static_cast<std::size_t>(column)];
                    const double lower=hessenberg[
                        static_cast<std::size_t>(row+1)]
                        [static_cast<std::size_t>(column)];
                    hessenberg[static_cast<std::size_t>(row)]
                        [static_cast<std::size_t>(column)]=
                        cosine[static_cast<std::size_t>(row)]*upper
                        +sine[static_cast<std::size_t>(row)]*lower;
                    hessenberg[static_cast<std::size_t>(row+1)]
                        [static_cast<std::size_t>(column)]=
                        -sine[static_cast<std::size_t>(row)]*upper
                        +cosine[static_cast<std::size_t>(row)]*lower;
                }
                const double diagonal=hessenberg[
                    static_cast<std::size_t>(column)]
                    [static_cast<std::size_t>(column)];
                const double subdiagonal=hessenberg[
                    static_cast<std::size_t>(column+1)]
                    [static_cast<std::size_t>(column)];
                const double radius=std::hypot(
                    diagonal,subdiagonal);
                cosine[static_cast<std::size_t>(column)]=radius>0.0
                    ?diagonal/radius:1.0;
                sine[static_cast<std::size_t>(column)]=radius>0.0
                    ?subdiagonal/radius:0.0;
                hessenberg[static_cast<std::size_t>(column)]
                    [static_cast<std::size_t>(column)]=radius;
                hessenberg[static_cast<std::size_t>(column+1)]
                    [static_cast<std::size_t>(column)]=0.0;
                const double previous=transformed[
                    static_cast<std::size_t>(column)];
                transformed[static_cast<std::size_t>(column)]=
                    cosine[static_cast<std::size_t>(column)]*previous;
                transformed[static_cast<std::size_t>(column+1)]=
                    -sine[static_cast<std::size_t>(column)]*previous;
                used=column+1;
                if(std::abs(transformed[static_cast<std::size_t>(used)])
                        /krylovInitialNorm<=relativeTolerance)
                    break;
            }
            if(used<=0)
                throw std::runtime_error(
                    "structured pressure Schur Krylov basis is empty");
            std::vector<double> coefficient(
                static_cast<std::size_t>(used),0.0);
            for(int row=used-1;row>=0;--row)
            {
                double value=transformed[static_cast<std::size_t>(row)];
                for(int column=row+1;column<used;++column)
                    value-=hessenberg[static_cast<std::size_t>(row)]
                        [static_cast<std::size_t>(column)]
                        *coefficient[static_cast<std::size_t>(column)];
                const double diagonal=hessenberg[
                    static_cast<std::size_t>(row)]
                    [static_cast<std::size_t>(row)];
                if(!std::isfinite(diagonal)
                    ||std::abs(diagonal)<1.0e-30)
                    throw std::runtime_error(
                        "structured pressure Schur Krylov system is singular");
                coefficient[static_cast<std::size_t>(row)]=
                    value/diagonal;
            }
            for(int column=0;column<used;++column)
                addScaled(domain,activeMask,result.pressure,
                    basis[static_cast<std::size_t>(column)],
                    coefficient[static_cast<std::size_t>(column)]);
            result.krylovIterations+=used;

            Audit residualAudit;
            const std::vector<double> action=
                applyPressureOperatorWithStencil(
                    domain,mesh,activeMask,stencil,pressureDiffusion,
                    inletPressureDirichlet,pseudoCfl,
                    result.pressure,residualAudit);
            mergeAudit(result.audit,residualAudit);
            trueResidual=rightHandSide;
            addScaled(domain,activeMask,trueResidual,action,-1.0);
            result.finalResidualNorm=scalarNorm(
                domain,activeMask,trueResidual);
            result.residualRatio=result.finalResidualNorm
                /result.initialResidualNorm;
            if(!std::isfinite(result.residualRatio))
                throw std::runtime_error(
                    "structured pressure Schur residual is not finite");
        }
        if(result.residualRatio>relativeTolerance)
            throw std::runtime_error(
                "structured pressure Schur Krylov solve did not converge: ratio="
                +std::to_string(result.residualRatio));

        const std::vector<Vector> gradient=pressureGradient(
            domain,mesh,activeMask,inletPressureDirichlet,
            result.pressure,result.audit);
        const std::vector<Vector> pressureMomentum=momentumInverse(
            domain,mesh,activeMask,stencil,pseudoCfl,
            gradient,result.audit);
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active(domain,activeMask,cell)) continue;
            Vector momentum=predictedMomentum[
                static_cast<std::size_t>(cell)];
            for(int component=0;component<3;++component)
                momentum[static_cast<std::size_t>(component)]
                    -=pressureMomentum[static_cast<std::size_t>(cell)]
                        [static_cast<std::size_t>(component)];
            const double rho=domain.Q(cell).rho;
            if(!std::isfinite(rho)||rho<=0.0)
                throw std::runtime_error(
                    "structured pressure Schur density is invalid");
            result.velocity[static_cast<std::size_t>(cell)]={
                momentum[0]/rho,momentum[1]/rho,momentum[2]/rho};
            result.maximumRelativePressure=std::max(
                result.maximumRelativePressure,
                std::abs(result.pressure[static_cast<std::size_t>(cell)])
                    /std::max(std::abs(domain.W(cell).p),1.0));
            const double velocityMagnitude=std::sqrt(
                result.velocity[static_cast<std::size_t>(cell)][0]
                    *result.velocity[static_cast<std::size_t>(cell)][0]
                +result.velocity[static_cast<std::size_t>(cell)][1]
                    *result.velocity[static_cast<std::size_t>(cell)][1]
                +result.velocity[static_cast<std::size_t>(cell)][2]
                    *result.velocity[static_cast<std::size_t>(cell)][2]);
            result.maximumVelocityMach=std::max(
                result.maximumVelocityMach,
                velocityMagnitude/std::max(
                    domain.Q(cell).acoustic.soundSpeed(),1.0));
        }
        return result;
    }
};

} // End namespace Foam

#endif
