#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "../src/framework/SectorStructuredMomentumInverse.H"

using namespace Foam;

namespace
{

int index(const Domain& domain,int i,int j,int k)
{ return k*domain.Nx()*domain.Ny()+j*domain.Nx()+i; }

double vectorMagnitude(const std::array<double,3>& value)
{
    return std::sqrt(value[0]*value[0]+value[1]*value[1]
        +value[2]*value[2]);
}

} // namespace

int main()
{
    SectorMeshMetrics mesh(vortexProfilePath(),10,8,6,120.0);
    Domain domain(10,8,6,mesh.length(),0.065,mesh.sectorAngle(),2);
    std::vector<unsigned char> active(
        static_cast<std::size_t>(domain.N()),0U);
    StructuredImplicitStencil stencil;
    stencil.ensureSize(static_cast<std::size_t>(domain.N()));
    const int nx=domain.Nx(),ny=domain.Ny();
    const int plane=nx*ny;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const int cell=index(domain,i,j,k);
        active[static_cast<std::size_t>(cell)]=1U;
        stencil.diagonal[static_cast<std::size_t>(cell)]=
            8.0+0.05*i+0.03*j;
        stencil.face[StructuredImplicitStencil::I_MINUS]
            [static_cast<std::size_t>(cell)]=0.9;
        stencil.face[StructuredImplicitStencil::I_PLUS]
            [static_cast<std::size_t>(cell)]=1.0;
        stencil.face[StructuredImplicitStencil::J_MINUS]
            [static_cast<std::size_t>(cell)]=0.7;
        stencil.face[StructuredImplicitStencil::J_PLUS]
            [static_cast<std::size_t>(cell)]=0.8;
        stencil.face[StructuredImplicitStencil::K_MINUS]
            [static_cast<std::size_t>(cell)]=0.6;
        stencil.face[StructuredImplicitStencil::K_PLUS]
            [static_cast<std::size_t>(cell)]=0.6;
    }

    using Vector=SectorStructuredMomentumInverse::Vector;
    std::vector<Vector> exact(
        static_cast<std::size_t>(domain.N()),Vector{0.0,0.0,0.0});
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const int cell=index(domain,i,j,k);
        const double theta=mesh.thetaCenter(k);
        // Cylindrical components are independent of theta, so the Cartesian
        // field is exactly rotationally periodic across the 120-degree seam.
        const double axial=0.4+0.03*i-0.02*j;
        const double radial=-0.2+0.01*i+0.04*j;
        const double tangential=0.15-0.02*i+0.01*j;
        exact[static_cast<std::size_t>(cell)]={
            axial,
            radial*std::cos(theta)-tangential*std::sin(theta),
            radial*std::sin(theta)+tangential*std::cos(theta)};
    }

    const double pseudoCfl=2.5;
    const std::vector<Vector> rightHandSide=
        SectorStructuredMomentumInverse::applyOperator(
            domain,mesh,active,stencil,pseudoCfl,exact);
    double maximumActionPeriodicCylindricalError=0.0;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    {
        const int lower=index(domain,i,j,1);
        const int upper=index(domain,i,j,domain.Nz()-2);
        const double lowerTheta=mesh.thetaCenter(1);
        const double upperTheta=mesh.thetaCenter(domain.Nz()-2);
        const Vector& lowerAction=rightHandSide[
            static_cast<std::size_t>(lower)];
        const Vector& upperAction=rightHandSide[
            static_cast<std::size_t>(upper)];
        const std::array<double,3> lowerCylindrical{
            lowerAction[0],
            lowerAction[1]*std::cos(lowerTheta)
                +lowerAction[2]*std::sin(lowerTheta),
            -lowerAction[1]*std::sin(lowerTheta)
                +lowerAction[2]*std::cos(lowerTheta)};
        const std::array<double,3> upperCylindrical{
            upperAction[0],
            upperAction[1]*std::cos(upperTheta)
                +upperAction[2]*std::sin(upperTheta),
            -upperAction[1]*std::sin(upperTheta)
                +upperAction[2]*std::cos(upperTheta)};
        for(int component=0;component<3;++component)
            maximumActionPeriodicCylindricalError=std::max(
                maximumActionPeriodicCylindricalError,
                std::abs(lowerCylindrical[
                    static_cast<std::size_t>(component)]
                    -upperCylindrical[
                        static_cast<std::size_t>(component)]));
    }
    const auto solved=SectorStructuredMomentumInverse::solve(
        domain,mesh,active,stencil,pseudoCfl,
        rightHandSide,160,1.0e-12);

    double maximumRelativeError=0.0;
    double maximumPeriodicCylindricalError=0.0;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    {
        const int reference=index(domain,i,j,1);
        const double referenceTheta=mesh.thetaCenter(1);
        const Vector& referenceValue=
            solved.correction[static_cast<std::size_t>(reference)];
        const std::array<double,3> referenceCylindrical{
            referenceValue[0],
            referenceValue[1]*std::cos(referenceTheta)
                +referenceValue[2]*std::sin(referenceTheta),
            -referenceValue[1]*std::sin(referenceTheta)
                +referenceValue[2]*std::cos(referenceTheta)};
        for(int k=1;k<domain.Nz()-1;++k)
        {
            const int cell=index(domain,i,j,k);
            const Vector& expected=exact[static_cast<std::size_t>(cell)];
            const Vector& actual=solved.correction[
                static_cast<std::size_t>(cell)];
            Vector error{
                actual[0]-expected[0],
                actual[1]-expected[1],
                actual[2]-expected[2]};
            maximumRelativeError=std::max(maximumRelativeError,
                vectorMagnitude(error)
                /std::max(vectorMagnitude(expected),1.0e-300));

            // At the two periodic edge planes, the cylindrical components
            // must match when the manufactured phase is sector-periodic.
            if(k==domain.Nz()-2)
            {
                const double theta=mesh.thetaCenter(k);
                const std::array<double,3> cylindrical{
                    actual[0],
                    actual[1]*std::cos(theta)
                        +actual[2]*std::sin(theta),
                    -actual[1]*std::sin(theta)
                        +actual[2]*std::cos(theta)};
                for(int component=0;component<3;++component)
                    maximumPeriodicCylindricalError=std::max(
                        maximumPeriodicCylindricalError,
                        std::abs(cylindrical[
                            static_cast<std::size_t>(component)]
                            -referenceCylindrical[
                                static_cast<std::size_t>(component)]));
            }
        }
    }

    if(solved.residualRatio>1.0e-11
        ||maximumRelativeError>1.0e-10
        ||maximumPeriodicCylindricalError>1.0e-10
        ||maximumActionPeriodicCylindricalError>1.0e-10
        ||!(solved.minimumDiagonalDominance>0.0)
        ||!(solved.maximumFaceCoupling>0.0))
    {
        std::cerr << "structured momentum inverse failed: sweeps/ratio/error/periodic="
                  << solved.sweeps << "/" << solved.residualRatio << "/"
                  << maximumRelativeError << "/"
                  << maximumPeriodicCylindricalError
                  << " action-periodic="
                  << maximumActionPeriodicCylindricalError
                  << " dominance=" << solved.minimumDiagonalDominance
                  << "\n";
        return 1;
    }

    std::cout << "sweeps/ratio/error/periodic="
              << solved.sweeps << "/" << solved.residualRatio << "/"
              << maximumRelativeError << "/"
              << maximumPeriodicCylindricalError
              << " action-periodic="
              << maximumActionPeriodicCylindricalError
              << " diagonal=" << solved.minimumDiagonal << "/"
              << solved.maximumDiagonal
              << " face/dominance=" << solved.maximumFaceCoupling << "/"
              << solved.minimumDiagonalDominance
              << "\nSECTOR_STRUCTURED_MOMENTUM_INVERSE_PASS\n";
    return 0;
}
