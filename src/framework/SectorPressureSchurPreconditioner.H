#ifndef SectorPressureSchurPreconditioner_H
#define SectorPressureSchurPreconditioner_H

/*---------------------------------------------------------------------------*\
    SectorPressureSchurPreconditioner - pressure/velocity low-frequency solve.

    A local pseudo-time elimination of density and momentum gives the positive
    pressure Helmholtz approximation

      [1/(a^2 dtau) - div(dtau_f grad)] dp
          = F_rho - div((dtau F_m)_f),
      du = dtau/rho [F_m - grad(dp)].

    The graph uses mapped face areas, normal center distances and cell volumes.
    Interior/periodic faces share one arithmetic pseudo-time mobility so the
    pressure flux is conservative when neighboring local time scales differ.
    Azimuthal faces wrap pressure as a scalar while momentum predictors rotate
    through the sector angle in the same Cartesian convention as Domain's
    cyclic halo. Physical-boundary correction values remain homogeneous; the
    nonlinear boundary hook is authoritative when a trial is evaluated.
\*---------------------------------------------------------------------------*/

#include "Operator.H"
#include "SectorMeshMetrics.H"
#include "SectorStructuredMomentumInverse.H"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Foam
{

class SectorPressureSchurPreconditioner
{
public:
    struct Result
    {
        std::vector<double> pressure;
        std::vector<std::array<double,3>> velocity;
        int sweeps=0;
        double initialResidualNorm=0.0;
        double finalResidualNorm=0.0;
        double residualRatio=1.0;
        double maximumRelativePressure=0.0;
        double maximumVelocityMach=0.0;
        double minimumFacePseudoTime=0.0;
        double maximumFacePseudoTime=0.0;
        double maximumFaceReciprocityError=0.0;
        double maximumFaceNonOrthogonalityDegrees=0.0;
        int leastSquaresGradientApplications=0;
        int leastSquaresGradientFallbacks=0;
        int nonOrthogonalCorrections=0;
        double nonOrthogonalResidualRatio=1.0;
        int structuredSchurIterations=0;
        double structuredSchurResidualRatio=1.0;
        int structuredMomentumSolves=0;
        int structuredMomentumSweeps=0;
        double structuredMomentumMaximumResidualRatio=0.0;
        int inletPressureDirichletFaces=0;
        int inletDirichletFaces=0;
        int inletNeumannFaces=0;
    };

private:
    struct System
    {
        std::vector<double> diagonal;
        std::vector<double> face[StructuredImplicitStencil::FACE_COUNT];
        std::vector<double> facePseudoTime[
            StructuredImplicitStencil::FACE_COUNT];
        std::vector<double> pseudoTime;
        double minimumFacePseudoTime=std::numeric_limits<double>::infinity();
        double maximumFacePseudoTime=0.0;
        double maximumFaceReciprocityError=0.0;
        double maximumFaceNonOrthogonalityDegrees=0.0;
        int inletPressureDirichletFaces=0;
        int inletDirichletFaces=0;
        int inletNeumannFaces=0;
    };

    struct StructuredMomentumAudit
    {
        int solves=0;
        int sweeps=0;
        double maximumResidualRatio=0.0;
    };

    static int index(const Domain& domain,int i,int j,int k)
    { return k*domain.Nx()*domain.Ny()+j*domain.Nx()+i; }

    static bool active(const Domain& domain,int cell)
    {
        if(cell<0||cell>=domain.N()) return false;
        const int nx=domain.Nx(),ny=domain.Ny();
        const int i=cell%nx;
        const int j=(cell/nx)%ny;
        const int k=cell/(nx*ny);
        return i>0&&i<nx-1&&j>0&&j<ny-1&&k>0&&k<domain.Nz()-1;
    }

    static double magnitude(const FaceAreaVector& area)
    { return std::sqrt(area.x*area.x+area.y*area.y+area.z*area.z); }

    static double normalDistance(
        const FaceAreaVector& center,
        const FaceAreaVector& neighbor,
        const FaceAreaVector& area)
    {
        const double areaMagnitude=magnitude(area);
        if(areaMagnitude<=1.0e-30) return 0.0;
        return std::abs((neighbor.x-center.x)*area.x
            +(neighbor.y-center.y)*area.y
            +(neighbor.z-center.z)*area.z)/areaMagnitude;
    }

    static double nonOrthogonalityDegrees(
        const FaceAreaVector& center,
        const FaceAreaVector& neighbor,
        const FaceAreaVector& area)
    {
        const double dx=neighbor.x-center.x;
        const double dy=neighbor.y-center.y;
        const double dz=neighbor.z-center.z;
        const double distance=std::sqrt(dx*dx+dy*dy+dz*dz);
        const double areaMagnitude=magnitude(area);
        if(distance<=1.0e-30||areaMagnitude<=1.0e-30) return 0.0;
        const double cosine=std::clamp(std::abs(
            (dx*area.x+dy*area.y+dz*area.z)
                /(distance*areaMagnitude)),0.0,1.0);
        return std::acos(cosine)*180.0/M_PI;
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

    static FaceAreaVector geometryNeighborCenter(
        const SectorMeshMetrics& mesh,int i,int j,int k,int face)
    {
        switch(face)
        {
            case StructuredImplicitStencil::I_MINUS:
                return mesh.cellCenter(i-1,j,k);
            case StructuredImplicitStencil::I_PLUS:
                return mesh.cellCenter(i+1,j,k);
            case StructuredImplicitStencil::J_MINUS:
                return mesh.cellCenter(i,j-1,k);
            case StructuredImplicitStencil::J_PLUS:
                return mesh.cellCenter(i,j+1,k);
            case StructuredImplicitStencil::K_MINUS:
                return mesh.cellCenter(i,j,k-1);
            case StructuredImplicitStencil::K_PLUS:
                return mesh.cellCenter(i,j,k+1);
            default:return mesh.cellCenter(i,j,k);
        }
    }

    static std::array<double,3> momentumInCellFrame(
        const SectorMeshMetrics& mesh,const Domain& domain,
        int k,int face,const std::vector<double>& momentum)
    {
        double angle=0.0;
        if(face==StructuredImplicitStencil::K_MINUS&&k==1)
            angle=-mesh.sectorAngle();
        else if(face==StructuredImplicitStencil::K_PLUS
            &&k==domain.Nz()-2)
            angle=mesh.sectorAngle();
        if(momentum.size()!=3U)
            throw std::runtime_error(
                "pressure Schur momentum vector is not three-dimensional");
        if(angle==0.0) return {momentum[0],momentum[1],momentum[2]};
        const double cosine=std::cos(angle);
        const double sine=std::sin(angle);
        return {momentum[0],
            cosine*momentum[1]-sine*momentum[2],
            sine*momentum[1]+cosine*momentum[2]};
    }

    static bool solveThreeByThree(
        std::array<std::array<double,3>,3> matrix,
        std::array<double,3> rhs,std::array<double,3>& solution)
    {
        double matrixMaximum=0.0;
        for(const auto& row:matrix)
            for(double value:row)
                matrixMaximum=std::max(matrixMaximum,std::abs(value));
        if(!std::isfinite(matrixMaximum)||matrixMaximum<=0.0) return false;
        for(int column=0;column<3;++column)
        {
            int pivot=column;
            for(int row=column+1;row<3;++row)
                if(std::abs(matrix[row][column])
                    >std::abs(matrix[pivot][column])) pivot=row;
            if(!std::isfinite(matrix[pivot][column])
                ||std::abs(matrix[pivot][column])<1.0e-12*matrixMaximum)
                return false;
            if(pivot!=column)
            {
                std::swap(matrix[pivot],matrix[column]);
                std::swap(rhs[pivot],rhs[column]);
            }
            const double diagonal=matrix[column][column];
            for(int entry=column;entry<3;++entry)
                matrix[column][entry]/=diagonal;
            rhs[column]/=diagonal;
            for(int row=0;row<3;++row)
            {
                if(row==column) continue;
                const double factor=matrix[row][column];
                for(int entry=column;entry<3;++entry)
                    matrix[row][entry]-=factor*matrix[column][entry];
                rhs[row]-=factor*rhs[column];
            }
        }
        solution=rhs;
        return std::isfinite(solution[0])&&std::isfinite(solution[1])
            &&std::isfinite(solution[2]);
    }

    static FaceAreaVector greenGaussPressureGradient(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const std::vector<unsigned char>& inletPressureDirichlet,
        const std::vector<double>& pressure,int cell)
    {
        const int nx=domain.Nx(),ny=domain.Ny();
        const int i=cell%nx;
        const int j=(cell/nx)%ny;
        const int k=cell/(nx*ny);
        const double volume=mesh.volume(i,j,k);
        FaceAreaVector gradient{0.0,0.0,0.0};
        for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
        {
            const int neighbor=wrappedNeighbor(domain,i,j,k,face);
            const bool neighborActive=active(domain,neighbor)
                &&activeMask[static_cast<std::size_t>(neighbor)]!=0U;
            const double cellPressure=pressure[static_cast<std::size_t>(cell)];
            const bool inletDirichlet=
                face==StructuredImplicitStencil::I_MINUS
                &&inletPressureDirichlet[static_cast<std::size_t>(
                    k*domain.Ny()+j)]!=0U;
            const double facePressure=neighborActive
                ?0.5*(cellPressure+pressure[
                    static_cast<std::size_t>(neighbor)])
                :(inletDirichlet?0.0:cellPressure);
            const FaceAreaVector area=faceArea(mesh,i,j,k,face);
            gradient.x+=facePressure*area.x/volume;
            gradient.y+=facePressure*area.y/volume;
            gradient.z+=facePressure*area.z/volume;
        }
        return gradient;
    }

    static bool leastSquaresPressureGradient(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const std::vector<unsigned char>& inletPressureDirichlet,
        const std::vector<double>& pressure,int cell,
        FaceAreaVector& gradient)
    {
        const int nx=domain.Nx(),ny=domain.Ny();
        const int i=cell%nx;
        const int j=(cell/nx)%ny;
        const int k=cell/(nx*ny);
        const FaceAreaVector center=mesh.cellCenter(i,j,k);
        const double centerPressure=pressure[static_cast<std::size_t>(cell)];
        std::array<std::array<double,3>,3> matrix{};
        std::array<double,3> rhs{};
        int observations=0;
        for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
        {
            const int neighbor=wrappedNeighbor(domain,i,j,k,face);
            const bool neighborActive=active(domain,neighbor)
                &&activeMask[static_cast<std::size_t>(neighbor)]!=0U;
            const bool inletDirichlet=!neighborActive
                &&face==StructuredImplicitStencil::I_MINUS
                &&inletPressureDirichlet[static_cast<std::size_t>(
                    k*domain.Ny()+j)]!=0U;
            if(!neighborActive&&!inletDirichlet) continue;
            FaceAreaVector sample=geometryNeighborCenter(mesh,i,j,k,face);
            double samplePressure=neighborActive
                ?pressure[static_cast<std::size_t>(neighbor)]:0.0;
            if(inletDirichlet)
            {
                sample=mesh.cellCenter(i,j,k);
                sample.x=mesh.xFace(i);
            }
            const std::array<double,3> displacement{
                sample.x-center.x,sample.y-center.y,sample.z-center.z};
            const double distanceSquared=
                displacement[0]*displacement[0]
                +displacement[1]*displacement[1]
                +displacement[2]*displacement[2];
            if(!std::isfinite(distanceSquared)||distanceSquared<=1.0e-30)
                continue;
            const double difference=samplePressure-centerPressure;
            const double weight=1.0/distanceSquared;
            for(int row=0;row<3;++row)
            {
                rhs[row]+=weight*displacement[row]*difference;
                for(int column=0;column<3;++column)
                    matrix[row][column]+=weight
                        *displacement[row]*displacement[column];
            }
            ++observations;
        }
        std::array<double,3> solution{};
        if(observations<3||!solveThreeByThree(matrix,rhs,solution))
            return false;
        gradient={solution[0],solution[1],solution[2]};
        return true;
    }

    static FaceAreaVector gradientInCellFrame(
        const SectorMeshMetrics& mesh,const Domain& domain,
        int k,int face,const FaceAreaVector& gradient)
    {
        double angle=0.0;
        if(face==StructuredImplicitStencil::K_MINUS&&k==1)
            angle=-mesh.sectorAngle();
        else if(face==StructuredImplicitStencil::K_PLUS
            &&k==domain.Nz()-2)
            angle=mesh.sectorAngle();
        if(angle==0.0) return gradient;
        const double cosine=std::cos(angle);
        const double sine=std::sin(angle);
        return {gradient.x,
            cosine*gradient.y-sine*gradient.z,
            sine*gradient.y+cosine*gradient.z};
    }

    static std::vector<FaceAreaVector> pressureGradientField(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const std::vector<unsigned char>& inletPressureDirichlet,
        const std::vector<double>& pressure,bool useLeastSquares,
        int* applications=nullptr,int* fallbacks=nullptr)
    {
        std::vector<FaceAreaVector> gradients(
            static_cast<std::size_t>(domain.N()),
            FaceAreaVector{0.0,0.0,0.0});
        int applied=0,fallback=0;
        for(int cell=0;cell<domain.N();++cell)
        {
            if(activeMask[static_cast<std::size_t>(cell)]==0U) continue;
            FaceAreaVector gradient{0.0,0.0,0.0};
            const bool leastSquaresApplied=useLeastSquares
                &&leastSquaresPressureGradient(
                    domain,mesh,activeMask,inletPressureDirichlet,
                    pressure,cell,gradient);
            if(leastSquaresApplied) ++applied;
            else
            {
                if(useLeastSquares) ++fallback;
                gradient=greenGaussPressureGradient(
                    domain,mesh,activeMask,inletPressureDirichlet,
                    pressure,cell);
            }
            gradients[static_cast<std::size_t>(cell)]=gradient;
        }
        if(applications) *applications=applied;
        if(fallbacks) *fallbacks=fallback;
        return gradients;
    }

    static System assemble(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const StructuredImplicitStencil& spectralStencil,
        const std::vector<unsigned char>& inletPressureDirichlet,
        double pseudoCfl,bool zeroInletCorrectionMassFlux)
    {
        if(activeMask.size()!=static_cast<std::size_t>(domain.N())
            ||spectralStencil.diagonal.size()
                !=static_cast<std::size_t>(domain.N())
            ||inletPressureDirichlet.size()!=static_cast<std::size_t>(
                domain.Ny()*domain.Nz())
            ||!std::isfinite(pseudoCfl)||pseudoCfl<=0.0)
            throw std::runtime_error(
                "pressure Schur inputs do not match Domain");
        System system;
        system.diagonal.assign(static_cast<std::size_t>(domain.N()),0.0);
        system.pseudoTime.assign(static_cast<std::size_t>(domain.N()),0.0);
        for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
        {
            system.face[face].assign(
                static_cast<std::size_t>(domain.N()),0.0);
            system.facePseudoTime[face].assign(
                static_cast<std::size_t>(domain.N()),0.0);
        }

        const int nx=domain.Nx(),ny=domain.Ny();
        const int plane=nx*ny;
        for(int cell=0;cell<domain.N();++cell)
        {
            if(activeMask[static_cast<std::size_t>(cell)]==0U) continue;
            const double spectral=spectralStencil.diagonal[
                static_cast<std::size_t>(cell)];
            const double dtau=std::isfinite(spectral)&&spectral>0.0
                ?pseudoCfl/spectral:0.0;
            const double sound=domain.Q(cell).acoustic.soundSpeed();
            if(!std::isfinite(dtau)||dtau<=0.0
                ||!std::isfinite(sound)||sound<=0.0)
                throw std::runtime_error(
                    "pressure Schur local time/acoustics are invalid");
            system.pseudoTime[static_cast<std::size_t>(cell)]=dtau;
            system.diagonal[static_cast<std::size_t>(cell)]=
                1.0/(sound*sound*dtau);
        }
        for(int cell=0;cell<domain.N();++cell)
        {
            if(activeMask[static_cast<std::size_t>(cell)]==0U) continue;
            const int i=cell%nx;
            const int j=(cell/nx)%ny;
            const int k=cell/plane;
            const double cellPseudoTime=system.pseudoTime[
                static_cast<std::size_t>(cell)];
            double diagonal=system.diagonal[static_cast<std::size_t>(cell)];
            const double volume=mesh.volume(i,j,k);
            const FaceAreaVector center=mesh.cellCenter(i,j,k);
            for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
            {
                const int neighbor=wrappedNeighbor(domain,i,j,k,face);
                const bool neighborActive=active(domain,neighbor)
                    &&activeMask[static_cast<std::size_t>(neighbor)]!=0U;
                const bool physicalInlet=
                    !neighborActive
                    &&face==StructuredImplicitStencil::I_MINUS;
                const bool pressureDirichlet=physicalInlet
                    &&inletPressureDirichlet[static_cast<std::size_t>(
                        k*domain.Ny()+j)]!=0U;
                if(pressureDirichlet)
                    ++system.inletPressureDirichletFaces;
                const bool continuityDirichlet=pressureDirichlet
                    &&!zeroInletCorrectionMassFlux;
                double facePseudoTime=0.0;
                if(neighborActive)
                    facePseudoTime=0.5*(cellPseudoTime
                        +system.pseudoTime[static_cast<std::size_t>(neighbor)]);
                else if(continuityDirichlet)
                    facePseudoTime=cellPseudoTime;
                else
                {
                    if(physicalInlet) ++system.inletNeumannFaces;
                    continue;
                }
                const FaceAreaVector area=faceArea(mesh,i,j,k,face);
                const FaceAreaVector neighborCenter=
                    geometryNeighborCenter(mesh,i,j,k,face);
                system.maximumFaceNonOrthogonalityDegrees=std::max(
                    system.maximumFaceNonOrthogonalityDegrees,
                    nonOrthogonalityDegrees(center,neighborCenter,area));
                double distance=normalDistance(center,neighborCenter,area);
                if(!neighborActive
                    &&face==StructuredImplicitStencil::I_MINUS)
                    distance=std::abs(center.x-mesh.xFace(i));
                const double coefficient=facePseudoTime*magnitude(area)
                    /(std::max(distance,1.0e-30)*volume);
                if(!std::isfinite(facePseudoTime)||facePseudoTime<=0.0
                    ||!std::isfinite(coefficient)||coefficient<0.0)
                    throw std::runtime_error(
                        "pressure Schur face mobility is invalid");
                diagonal+=coefficient;
                if(neighborActive)
                {
                    system.face[face][static_cast<std::size_t>(cell)]=
                        coefficient;
                    system.facePseudoTime[face][
                        static_cast<std::size_t>(cell)]=facePseudoTime;
                    system.minimumFacePseudoTime=std::min(
                        system.minimumFacePseudoTime,facePseudoTime);
                    system.maximumFacePseudoTime=std::max(
                        system.maximumFacePseudoTime,facePseudoTime);
                }
                else ++system.inletDirichletFaces;
                // Injection patches prescribe homogeneous pressure correction.
                // Reflecting inlet background, outlet, axis and wall extrapolate
                // pressure and therefore add no matrix flux.
            }
            system.diagonal[static_cast<std::size_t>(cell)]=diagonal;
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
            if(activeMask[static_cast<std::size_t>(cell)]==0U) continue;
            const int i=cell%nx;
            const int j=(cell/nx)%ny;
            const int k=cell/plane;
            for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
            {
                const int neighbor=wrappedNeighbor(domain,i,j,k,face);
                if(!active(domain,neighbor)
                    ||activeMask[static_cast<std::size_t>(neighbor)]==0U)
                    continue;
                const double left=system.face[face][
                    static_cast<std::size_t>(cell)]*mesh.volume(i,j,k);
                const int ni=neighbor%nx;
                const int nj=(neighbor/nx)%ny;
                const int nk=neighbor/plane;
                const double right=system.face[reciprocalFace[face]][
                    static_cast<std::size_t>(neighbor)]
                    *mesh.volume(ni,nj,nk);
                system.maximumFaceReciprocityError=std::max(
                    system.maximumFaceReciprocityError,
                    std::abs(left-right)/std::max(
                        {std::abs(left),std::abs(right),1.0e-300}));
            }
        }
        if(!std::isfinite(system.minimumFacePseudoTime))
            system.minimumFacePseudoTime=0.0;
        return system;
    }

    static std::vector<double> apply(
        const Domain& domain,const std::vector<unsigned char>& activeMask,
        const System& system,const std::vector<double>& pressure)
    {
        std::vector<double> result(static_cast<std::size_t>(domain.N()),0.0);
        const int nx=domain.Nx(),ny=domain.Ny();
        const int plane=nx*ny;
        for(int cell=0;cell<domain.N();++cell)
        {
            if(activeMask[static_cast<std::size_t>(cell)]==0U) continue;
            const int i=cell%nx;
            const int j=(cell/nx)%ny;
            const int k=cell/plane;
            double value=system.diagonal[static_cast<std::size_t>(cell)]
                *pressure[static_cast<std::size_t>(cell)];
            for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
            {
                const int neighbor=wrappedNeighbor(domain,i,j,k,face);
                if(!active(domain,neighbor)) continue;
                value-=system.face[face][static_cast<std::size_t>(cell)]
                    *pressure[static_cast<std::size_t>(neighbor)];
            }
            result[static_cast<std::size_t>(cell)]=value;
        }
        return result;
    }

    static double norm(
        const std::vector<unsigned char>& activeMask,
        const std::vector<double>& values)
    {
        double sum=0.0;
        std::size_t count=0;
        for(std::size_t cell=0;cell<values.size();++cell)
            if(activeMask[cell]!=0U)
            {
                sum+=values[cell]*values[cell];
                ++count;
            }
        return std::sqrt(sum/std::max<std::size_t>(count,1U));
    }

    static std::vector<double> nonOrthogonalEquationContribution(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const System& system,
        const std::vector<FaceAreaVector>& gradients)
    {
        std::vector<double> contribution(
            static_cast<std::size_t>(domain.N()),0.0);
        const int nx=domain.Nx(),ny=domain.Ny();
        const int plane=nx*ny;
        for(int cell=0;cell<domain.N();++cell)
        {
            if(activeMask[static_cast<std::size_t>(cell)]==0U) continue;
            const int i=cell%nx;
            const int j=(cell/nx)%ny;
            const int k=cell/plane;
            const FaceAreaVector center=mesh.cellCenter(i,j,k);
            const double inverseVolume=1.0/mesh.volume(i,j,k);
            for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
            {
                const int neighbor=wrappedNeighbor(domain,i,j,k,face);
                if(!active(domain,neighbor)
                    ||activeMask[static_cast<std::size_t>(neighbor)]==0U)
                    continue;
                const double facePseudoTime=system.facePseudoTime[face][
                    static_cast<std::size_t>(cell)];
                if(facePseudoTime<=0.0) continue;
                const FaceAreaVector area=faceArea(mesh,i,j,k,face);
                const FaceAreaVector neighborCenter=
                    geometryNeighborCenter(mesh,i,j,k,face);
                const std::array<double,3> displacement{
                    neighborCenter.x-center.x,
                    neighborCenter.y-center.y,
                    neighborCenter.z-center.z};
                const double projection=area.x*displacement[0]
                    +area.y*displacement[1]+area.z*displacement[2];
                const double areaSquared=area.x*area.x
                    +area.y*area.y+area.z*area.z;
                if(!std::isfinite(projection)||projection<=1.0e-30
                    ||!std::isfinite(areaSquared)||areaSquared<=1.0e-30)
                    throw std::runtime_error(
                        "pressure Schur nonorthogonal face geometry is invalid");
                const FaceAreaVector orthogonal{
                    displacement[0]*areaSquared/projection,
                    displacement[1]*areaSquared/projection,
                    displacement[2]*areaSquared/projection};
                const FaceAreaVector tangential{
                    area.x-orthogonal.x,
                    area.y-orthogonal.y,
                    area.z-orthogonal.z};
                const FaceAreaVector neighborGradient=gradientInCellFrame(
                    mesh,domain,k,face,
                    gradients[static_cast<std::size_t>(neighbor)]);
                const FaceAreaVector faceGradient{
                    0.5*(gradients[static_cast<std::size_t>(cell)].x
                        +neighborGradient.x),
                    0.5*(gradients[static_cast<std::size_t>(cell)].y
                        +neighborGradient.y),
                    0.5*(gradients[static_cast<std::size_t>(cell)].z
                        +neighborGradient.z)};
                const double nonOrthogonalFlux=facePseudoTime*inverseVolume
                    *(faceGradient.x*tangential.x
                        +faceGradient.y*tangential.y
                        +faceGradient.z*tangential.z);
                if(!std::isfinite(nonOrthogonalFlux))
                    throw std::runtime_error(
                        "pressure Schur nonorthogonal flux is not finite");
                // The full equation is A_orth(p) - div(dtau grad_t p)=rhs.
                contribution[static_cast<std::size_t>(cell)]
                    -=nonOrthogonalFlux;
            }
        }
        return contribution;
    }

public:
    static FaceAreaVector reconstructPressureGradient(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<unsigned char>& activeMask,
        const std::vector<unsigned char>& inletPressureDirichlet,
        const std::vector<double>& pressure,int cell,
        bool useLeastSquares,bool* leastSquaresApplied=nullptr)
    {
        if(activeMask.size()!=static_cast<std::size_t>(domain.N())
            ||pressure.size()!=static_cast<std::size_t>(domain.N())
            ||inletPressureDirichlet.size()!=static_cast<std::size_t>(
                domain.Ny()*domain.Nz())
            ||cell<0||cell>=domain.N()
            ||activeMask[static_cast<std::size_t>(cell)]==0U)
            throw std::runtime_error(
                "pressure Schur gradient inputs do not match Domain");
        FaceAreaVector gradient{0.0,0.0,0.0};
        const bool applied=useLeastSquares&&leastSquaresPressureGradient(
            domain,mesh,activeMask,inletPressureDirichlet,
            pressure,cell,gradient);
        if(!applied)
            gradient=greenGaussPressureGradient(
                domain,mesh,activeMask,inletPressureDirichlet,
                pressure,cell);
        if(leastSquaresApplied) *leastSquaresApplied=applied;
        return gradient;
    }

    static Result build(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<ConservedVars>& residual,
        const std::vector<unsigned char>& activeMask,
        const StructuredImplicitStencil& spectralStencil,
        const std::vector<unsigned char>& inletPressureDirichlet,
        double pseudoCfl,int maximumSweeps=120,
        double relativeTolerance=1.0e-8,
        bool useLeastSquaresGradient=false,
        int nonOrthogonalCorrections=0,
        bool zeroInletCorrectionMassFlux=false,
        int minimumSweeps=0)
    {
        if(residual.size()!=static_cast<std::size_t>(domain.N())
            ||maximumSweeps<=0||!std::isfinite(relativeTolerance)
            ||relativeTolerance<=0.0||relativeTolerance>=1.0
            ||nonOrthogonalCorrections<0||nonOrthogonalCorrections>12
            ||(nonOrthogonalCorrections>0&&!useLeastSquaresGradient)
            ||minimumSweeps<0||minimumSweeps>maximumSweeps)
            throw std::runtime_error("invalid pressure Schur controls");
        const System system=assemble(
            domain,mesh,activeMask,spectralStencil,
            inletPressureDirichlet,pseudoCfl,
            zeroInletCorrectionMassFlux);
        std::vector<double> rhs(static_cast<std::size_t>(domain.N()),0.0);
        const int nx=domain.Nx(),ny=domain.Ny();
        const int plane=nx*ny;
        for(int cell=0;cell<domain.N();++cell)
        {
            if(activeMask[static_cast<std::size_t>(cell)]==0U) continue;
            const int i=cell%nx;
            const int j=(cell/nx)%ny;
            const int k=cell/plane;
            double momentumDivergence=0.0;
            for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
            {
                const int neighbor=wrappedNeighbor(domain,i,j,k,face);
                const FaceAreaVector area=faceArea(mesh,i,j,k,face);
                const bool neighborActive=active(domain,neighbor)
                    &&activeMask[static_cast<std::size_t>(neighbor)]!=0U;
                // Correction mass flux is fixed to zero at injection, axis
                // and wall. At the steady extrapolating outlet use the cell
                // predictor; periodic/interior faces use a centered flux.
                if(!neighborActive
                    &&face!=StructuredImplicitStencil::I_PLUS) continue;
                const std::array<double,3> neighborMomentum=neighborActive
                    ?momentumInCellFrame(mesh,domain,k,face,
                        residual[static_cast<std::size_t>(neighbor)].rhoU)
                    :std::array<double,3>{0.0,0.0,0.0};
                for(int direction=0;direction<3;++direction)
                {
                    const double centerMomentum=residual[
                        static_cast<std::size_t>(cell)].rhoU[
                            static_cast<std::size_t>(direction)];
                    const double centerPredictor=system.pseudoTime[
                        static_cast<std::size_t>(cell)]*centerMomentum;
                    const double facePredictor=neighborActive
                        ?0.5*(centerPredictor+system.pseudoTime[
                            static_cast<std::size_t>(neighbor)]
                                *neighborMomentum[
                                    static_cast<std::size_t>(direction)])
                        :centerPredictor;
                    const double areaComponent=direction==0?area.x:
                        (direction==1?area.y:area.z);
                    momentumDivergence+=facePredictor*areaComponent;
                }
            }
            momentumDivergence/=mesh.volume(i,j,k);
            rhs[static_cast<std::size_t>(cell)]=
                residual[static_cast<std::size_t>(cell)].rho
                -momentumDivergence;
        }

        Result result;
        result.pressure.assign(static_cast<std::size_t>(domain.N()),0.0);
        result.velocity.assign(static_cast<std::size_t>(domain.N()),
            std::array<double,3>{0.0,0.0,0.0});
        result.minimumFacePseudoTime=system.minimumFacePseudoTime;
        result.maximumFacePseudoTime=system.maximumFacePseudoTime;
        result.maximumFaceReciprocityError=
            system.maximumFaceReciprocityError;
        result.maximumFaceNonOrthogonalityDegrees=
            system.maximumFaceNonOrthogonalityDegrees;
        result.inletPressureDirichletFaces=
            system.inletPressureDirichletFaces;
        result.inletDirichletFaces=system.inletDirichletFaces;
        result.inletNeumannFaces=system.inletNeumannFaces;
        result.initialResidualNorm=norm(activeMask,rhs);
        if(result.initialResidualNorm==0.0)
        {
            result.finalResidualNorm=0.0;
            result.residualRatio=0.0;
            return result;
        }

        auto solveOrthogonal=[&](const std::vector<double>& equationRhs)
        {
            const double equationNorm=std::max(
                norm(activeMask,equationRhs),1.0e-300);
            int sweepsThisSolve=0;
            for(int sweep=0;sweep<maximumSweeps;++sweep)
            {
                for(int pass=0;pass<2;++pass)
                {
                    const int begin=pass==0?0:domain.N()-1;
                    const int end=pass==0?domain.N():-1;
                    const int increment=pass==0?1:-1;
                    for(int cell=begin;cell!=end;cell+=increment)
                    {
                        if(activeMask[static_cast<std::size_t>(cell)]==0U)
                            continue;
                        const int i=cell%nx;
                        const int j=(cell/nx)%ny;
                        const int k=cell/plane;
                        double value=equationRhs[
                            static_cast<std::size_t>(cell)];
                        for(int face=0;
                            face<StructuredImplicitStencil::FACE_COUNT;
                            ++face)
                        {
                            const int neighbor=wrappedNeighbor(
                                domain,i,j,k,face);
                            if(!active(domain,neighbor)) continue;
                            value+=system.face[face][
                                static_cast<std::size_t>(cell)]
                                *result.pressure[
                                    static_cast<std::size_t>(neighbor)];
                        }
                        result.pressure[static_cast<std::size_t>(cell)]=
                            value/system.diagonal[
                                static_cast<std::size_t>(cell)];
                    }
                }
                ++result.sweeps;
                ++sweepsThisSolve;
                std::vector<double> equationResidual=equationRhs;
                const std::vector<double> action=apply(
                    domain,activeMask,system,result.pressure);
                for(std::size_t cell=0;
                    cell<equationResidual.size();++cell)
                    equationResidual[cell]-=action[cell];
                const double residualNorm=norm(
                    activeMask,equationResidual);
                if(!std::isfinite(residualNorm))
                    throw std::runtime_error(
                        "pressure Schur residual is not finite");
                if(sweepsThisSolve>=minimumSweeps
                    &&residualNorm/equationNorm<=relativeTolerance) break;
            }
        };

        std::vector<double> equationRhs=rhs;
        solveOrthogonal(equationRhs);
        auto fullEquationResidual=[&]()
        {
            std::vector<double> fullResidual=rhs;
            const std::vector<double> orthogonalAction=apply(
                domain,activeMask,system,result.pressure);
            for(std::size_t cell=0;cell<fullResidual.size();++cell)
                fullResidual[cell]-=orthogonalAction[cell];
            if(nonOrthogonalCorrections>0)
            {
                int applications=0,fallbacks=0;
                const std::vector<FaceAreaVector> gradients=
                    pressureGradientField(
                        domain,mesh,activeMask,inletPressureDirichlet,
                        result.pressure,true,&applications,&fallbacks);
                result.leastSquaresGradientApplications+=applications;
                result.leastSquaresGradientFallbacks+=fallbacks;
                const std::vector<double> nonOrthogonal=
                    nonOrthogonalEquationContribution(
                        domain,mesh,activeMask,system,gradients);
                for(std::size_t cell=0;cell<fullResidual.size();++cell)
                    fullResidual[cell]-=nonOrthogonal[cell];
            }
            return fullResidual;
        };
        std::vector<double> fullResidual=fullEquationResidual();
        result.finalResidualNorm=norm(activeMask,fullResidual);
        result.residualRatio=result.finalResidualNorm
            /result.initialResidualNorm;
        for(int correction=0;
            correction<nonOrthogonalCorrections
                &&result.residualRatio>relativeTolerance;
            ++correction)
        {
            int applications=0,fallbacks=0;
            const std::vector<FaceAreaVector> gradients=
                pressureGradientField(
                    domain,mesh,activeMask,inletPressureDirichlet,
                    result.pressure,true,&applications,&fallbacks);
            result.leastSquaresGradientApplications+=applications;
            result.leastSquaresGradientFallbacks+=fallbacks;
            const std::vector<double> nonOrthogonal=
                nonOrthogonalEquationContribution(
                    domain,mesh,activeMask,system,gradients);
            equationRhs=rhs;
            for(std::size_t cell=0;cell<equationRhs.size();++cell)
                equationRhs[cell]-=nonOrthogonal[cell];
            solveOrthogonal(equationRhs);
            ++result.nonOrthogonalCorrections;
            fullResidual=fullEquationResidual();
            result.finalResidualNorm=norm(activeMask,fullResidual);
            result.residualRatio=result.finalResidualNorm
                /result.initialResidualNorm;
        }
        result.nonOrthogonalResidualRatio=result.residualRatio;
        if(nonOrthogonalCorrections>0
            &&result.nonOrthogonalResidualRatio>relativeTolerance)
            throw std::runtime_error(
                "pressure Schur nonorthogonal correction did not converge");

        int finalGradientApplications=0,finalGradientFallbacks=0;
        const std::vector<FaceAreaVector> finalGradients=
            pressureGradientField(
                domain,mesh,activeMask,inletPressureDirichlet,
                result.pressure,useLeastSquaresGradient,
                &finalGradientApplications,&finalGradientFallbacks);
        result.leastSquaresGradientApplications+=finalGradientApplications;
        result.leastSquaresGradientFallbacks+=finalGradientFallbacks;
        for(int cell=0;cell<domain.N();++cell)
        {
            if(activeMask[static_cast<std::size_t>(cell)]==0U) continue;
            const int i=cell%nx;
            const int j=(cell/nx)%ny;
            const int k=cell/plane;
            const FaceAreaVector& gradient=finalGradients[
                static_cast<std::size_t>(cell)];
            const double factor=system.pseudoTime[
                static_cast<std::size_t>(cell)]/domain.Q(cell).rho;
            result.velocity[static_cast<std::size_t>(cell)]={
                factor*(residual[static_cast<std::size_t>(cell)].rhoU[0]
                    -gradient.x),
                factor*(residual[static_cast<std::size_t>(cell)].rhoU[1]
                    -gradient.y),
                factor*(residual[static_cast<std::size_t>(cell)].rhoU[2]
                    -gradient.z)};
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
                result.maximumVelocityMach,velocityMagnitude
                    /std::max(domain.Q(cell).acoustic.soundSpeed(),1.0));
        }
        return result;
    }
};

} // End namespace Foam

#endif
