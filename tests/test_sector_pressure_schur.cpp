#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "../src/framework/SectorPressureSchurPreconditioner.H"

using namespace Foam;

int main()
{
    SectorMeshMetrics mesh(vortexProfilePath(),10,8,6,120.0);
    Domain domain(10,8,6,mesh.length(),0.065,mesh.sectorAngle(),2);
    domain.setSpecies({"O2","C3H8"});
    for(int i=0;i<domain.Nx();++i)
    for(int j=0;j<domain.Ny();++j)
    for(int k=0;k<domain.Nz();++k)
    {
        PrimitiveVars state(2);
        state.p=4.0e6;
        state.T=700.0;
        state.U={20.0,0.0,0.0};
        state.Y={1.0,0.0};
        domain.stagedW(i,j,k)=state;
    }
    if(!domain.primitiveToConserved())
    {
        std::cerr << "pressure Schur test Domain initialization failed\n";
        return 1;
    }

    PrimitiveVars pressureDiffusionLeft(2),pressureDiffusionRight(2);
    pressureDiffusionLeft.p=3.9e6;
    pressureDiffusionLeft.T=680.0;
    pressureDiffusionLeft.U={31.0,4.0,-2.0};
    pressureDiffusionRight.p=4.1e6;
    pressureDiffusionRight.T=720.0;
    pressureDiffusionRight.U={-17.0,-3.0,5.0};
    const std::vector<double> pressureDiffusionNormal={1.0,0.0,0.0};
    const double pressureDiffusionForward=
        ausmPressureDiffusionMassMobility(
            pressureDiffusionLeft,pressureDiffusionRight,
            pressureDiffusionNormal,8.0,420.0,3.0,510.0);
    const std::vector<double> pressureDiffusionReverseNormal={-1.0,0.0,0.0};
    const double pressureDiffusionReverse=
        ausmPressureDiffusionMassMobility(
            pressureDiffusionRight,pressureDiffusionLeft,
            pressureDiffusionReverseNormal,3.0,510.0,8.0,420.0);
    if(!std::isfinite(pressureDiffusionForward)
        ||pressureDiffusionForward<=0.0
        ||std::abs(pressureDiffusionForward-pressureDiffusionReverse)
            >1.0e-14*std::max(pressureDiffusionForward,1.0e-300))
    {
        std::cerr << "AUSM pressure-diffusion reciprocity failed: "
                  << pressureDiffusionForward << "/"
                  << pressureDiffusionReverse << "\n";
        return 1;
    }

    PrimitiveVars finiteDifferenceLeft(2),finiteDifferenceRight(2);
    finiteDifferenceLeft.p=finiteDifferenceRight.p=4.0e6;
    finiteDifferenceLeft.T=finiteDifferenceRight.T=700.0;
    finiteDifferenceLeft.U=finiteDifferenceRight.U={18.0,0.0,0.0};
    finiteDifferenceLeft.Y=finiteDifferenceRight.Y={1.0,0.0};
    const double finiteDifferenceRho=7.0;
    const double finiteDifferenceSound=440.0;
    const double finiteDifferenceEnthalpy=8.0e5;
    const double expectedPressureMobility=
        ausmPressureDiffusionMassMobility(
            finiteDifferenceLeft,finiteDifferenceRight,
            pressureDiffusionNormal,
            finiteDifferenceRho,finiteDifferenceSound,
            finiteDifferenceRho,finiteDifferenceSound);
    const double pressurePerturbation=100.0;
    PrimitiveVars positiveLeft=finiteDifferenceLeft;
    PrimitiveVars positiveRight=finiteDifferenceRight;
    positiveLeft.p-=0.5*pressurePerturbation;
    positiveRight.p+=0.5*pressurePerturbation;
    PrimitiveVars negativeLeft=finiteDifferenceLeft;
    PrimitiveVars negativeRight=finiteDifferenceRight;
    negativeLeft.p+=0.5*pressurePerturbation;
    negativeRight.p-=0.5*pressurePerturbation;
    AUSMPlusUpFlux ausm;
    ConservedVars positiveFlux(2),negativeFlux(2);
    ausm.computeFluxCached(
        positiveLeft,positiveRight,pressureDiffusionNormal,positiveFlux,
        finiteDifferenceRho,finiteDifferenceSound,finiteDifferenceEnthalpy,
        finiteDifferenceRho,finiteDifferenceSound,finiteDifferenceEnthalpy);
    ausm.computeFluxCached(
        negativeLeft,negativeRight,pressureDiffusionNormal,negativeFlux,
        finiteDifferenceRho,finiteDifferenceSound,finiteDifferenceEnthalpy,
        finiteDifferenceRho,finiteDifferenceSound,finiteDifferenceEnthalpy);
    const double measuredPressureDerivative=
        (positiveFlux.rho-negativeFlux.rho)
        /(2.0*pressurePerturbation);
    const double pressureDerivativeRelativeError=std::abs(
        measuredPressureDerivative+expectedPressureMobility)
        /std::max(std::abs(expectedPressureMobility),1.0e-300);
    if(!std::isfinite(pressureDerivativeRelativeError)
        ||pressureDerivativeRelativeError>1.0e-10)
    {
        std::cerr << "AUSM pressure-diffusion finite-difference contract failed: "
                  << measuredPressureDerivative << "/"
                  << -expectedPressureMobility << "/"
                  << pressureDerivativeRelativeError << "\n";
        return 1;
    }

    PrimitiveVars variableDensityLeft=finiteDifferenceLeft;
    PrimitiveVars variableDensityRight=finiteDifferenceRight;
    variableDensityLeft.U[0]=45.0;
    variableDensityRight.U[0]=20.0;
    const double variableDensityLeftRho=9.0;
    const double variableDensityRightRho=2.0;
    const double variableExpectedMobility=
        ausmPressureDiffusionMassMobility(
            variableDensityLeft,variableDensityRight,
            pressureDiffusionNormal,
            variableDensityLeftRho,finiteDifferenceSound,
            variableDensityRightRho,finiteDifferenceSound);
    positiveLeft=variableDensityLeft;
    positiveRight=variableDensityRight;
    positiveLeft.p-=0.5*pressurePerturbation;
    positiveRight.p+=0.5*pressurePerturbation;
    negativeLeft=variableDensityLeft;
    negativeRight=variableDensityRight;
    negativeLeft.p+=0.5*pressurePerturbation;
    negativeRight.p-=0.5*pressurePerturbation;
    ausm.computeFluxCached(
        positiveLeft,positiveRight,pressureDiffusionNormal,positiveFlux,
        variableDensityLeftRho,finiteDifferenceSound,
        finiteDifferenceEnthalpy,variableDensityRightRho,
        finiteDifferenceSound,finiteDifferenceEnthalpy);
    ausm.computeFluxCached(
        negativeLeft,negativeRight,pressureDiffusionNormal,negativeFlux,
        variableDensityLeftRho,finiteDifferenceSound,
        finiteDifferenceEnthalpy,variableDensityRightRho,
        finiteDifferenceSound,finiteDifferenceEnthalpy);
    const double variableMeasuredDerivative=
        (positiveFlux.rho-negativeFlux.rho)
        /(2.0*pressurePerturbation);
    const double variableDerivativeRelativeError=std::abs(
        variableMeasuredDerivative+variableExpectedMobility)
        /std::max(std::abs(variableExpectedMobility),1.0e-300);
    if(!std::isfinite(variableDerivativeRelativeError)
        ||variableDerivativeRelativeError>1.0e-10)
    {
        std::cerr << "variable-density AUSM pressure derivative failed: "
                  << variableMeasuredDerivative << "/"
                  << -variableExpectedMobility << "/"
                  << variableDerivativeRelativeError << "\n";
        return 1;
    }

    std::vector<unsigned char> active(
        static_cast<std::size_t>(domain.N()),0U);
    std::vector<ConservedVars> residual(
        static_cast<std::size_t>(domain.N()),ConservedVars(2));
    StructuredImplicitStencil spectral;
    spectral.ensureSize(static_cast<std::size_t>(domain.N()));
    const int nx=domain.Nx(),ny=domain.Ny();
    const int plane=nx*ny;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const int cell=k*plane+j*nx+i;
        active[static_cast<std::size_t>(cell)]=1U;
        spectral.diagonal[static_cast<std::size_t>(cell)]=1.0e8;
        const double source=1.0e3*(
            std::sin(0.6*i)+0.3*std::cos(0.7*j));
        residual[static_cast<std::size_t>(cell)].rho=source;
    }

    std::vector<unsigned char> allDirichletInlet(
        static_cast<std::size_t>(domain.Ny()*domain.Nz()),0U);
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
        allDirichletInlet[static_cast<std::size_t>(k*domain.Ny()+j)]=1U;
    std::vector<unsigned char> singleActive(
        static_cast<std::size_t>(domain.N()),0U);
    std::vector<ConservedVars> singleResidual(
        static_cast<std::size_t>(domain.N()),ConservedVars(2));
    StructuredImplicitStencil singleSpectral;
    singleSpectral.ensureSize(static_cast<std::size_t>(domain.N()));
    std::vector<unsigned char> singleInlet(
        static_cast<std::size_t>(domain.Ny()*domain.Nz()),0U);
    const int singleI=1,singleJ=3,singleK=2;
    const int singleCell=singleK*plane+singleJ*nx+singleI;
    singleActive[static_cast<std::size_t>(singleCell)]=1U;
    singleResidual[static_cast<std::size_t>(singleCell)].rho=2.5e3;
    singleSpectral.diagonal[static_cast<std::size_t>(singleCell)]=1.0e8;
    singleInlet[static_cast<std::size_t>(
        singleK*domain.Ny()+singleJ)]=1U;
    const double singlePseudoTime=5.0e-8;
    const double singleSound=domain.Q(singleCell).acoustic.soundSpeed();
    const double singleVolume=mesh.volume(singleI,singleJ,singleK);
    const double singleArea=mesh.axialArea(singleI,singleJ,singleK);
    const double singleDistance=mesh.cellCenter(
        singleI,singleJ,singleK).x-mesh.xFace(singleI);
    const double singleDiagonal=1.0/(
        singleSound*singleSound*singlePseudoTime)
        +singlePseudoTime*singleArea/(singleDistance*singleVolume);
    const double expectedSinglePressure=2.5e3/singleDiagonal;
    const double expectedSingleVelocity=-singlePseudoTime
        *expectedSinglePressure*singleArea
        /(domain.Q(singleCell).rho*singleVolume);
    const auto singleResult=SectorPressureSchurPreconditioner::build(
        domain,mesh,singleResidual,singleActive,singleSpectral,
        singleInlet,5.0,40,1.0e-12);
    const auto singleLeastSquaresResult=
        SectorPressureSchurPreconditioner::build(
            domain,mesh,singleResidual,singleActive,singleSpectral,
            singleInlet,5.0,40,1.0e-12,true);
    const auto singleControlledInletResult=
        SectorPressureSchurPreconditioner::build(
            domain,mesh,singleResidual,singleActive,singleSpectral,
            singleInlet,5.0,40,1.0e-12,false,0,true);
    const double expectedControlledPressure=2.5e3
        *singleSound*singleSound*singlePseudoTime;
    const double expectedControlledVelocity=-singlePseudoTime
        *expectedControlledPressure*singleArea
        /(domain.Q(singleCell).rho*singleVolume);
    const double controlledPressureRelativeError=std::abs(
        singleControlledInletResult.pressure[
            static_cast<std::size_t>(singleCell)]
            -expectedControlledPressure)
        /std::max(std::abs(expectedControlledPressure),1.0e-300);
    const double controlledVelocityRelativeError=std::abs(
        singleControlledInletResult.velocity[
            static_cast<std::size_t>(singleCell)][0]
            -expectedControlledVelocity)
        /std::max(std::abs(expectedControlledVelocity),1.0e-300);
    const double singlePressureRelativeError=std::abs(
        singleResult.pressure[static_cast<std::size_t>(singleCell)]
            -expectedSinglePressure)
        /std::max(std::abs(expectedSinglePressure),1.0e-300);
    const double singleVelocityRelativeError=std::abs(
        singleResult.velocity[static_cast<std::size_t>(singleCell)][0]
            -expectedSingleVelocity)
        /std::max(std::abs(expectedSingleVelocity),1.0e-300);
    if(singlePressureRelativeError>1.0e-12
        ||singleVelocityRelativeError>1.0e-12
        ||singleResult.inletPressureDirichletFaces!=1
        ||singleResult.inletDirichletFaces!=1
        ||singleResult.inletNeumannFaces!=0
        ||controlledPressureRelativeError>1.0e-12
        ||controlledVelocityRelativeError>1.0e-12
        ||singleControlledInletResult.inletPressureDirichletFaces!=1
        ||singleControlledInletResult.inletDirichletFaces!=0
        ||singleControlledInletResult.inletNeumannFaces!=1
        ||singleLeastSquaresResult.leastSquaresGradientApplications!=0
        ||singleLeastSquaresResult.leastSquaresGradientFallbacks!=1
        ||std::abs(singleLeastSquaresResult.velocity[
            static_cast<std::size_t>(singleCell)][0]
            -singleResult.velocity[static_cast<std::size_t>(singleCell)][0])
                >1.0e-12*std::max(std::abs(expectedSingleVelocity),1.0))
    {
        std::cerr << "inlet-face Dirichlet Schur contract failed: pressure="
                  << singlePressureRelativeError << " velocity="
                  << singleVelocityRelativeError
                  << " controlled=" << controlledPressureRelativeError
                  << "/" << controlledVelocityRelativeError
                  << " pressure-D/continuity-D/N="
                  << singleControlledInletResult.inletPressureDirichletFaces
                  << "/"
                  << singleControlledInletResult.inletDirichletFaces << "/"
                  << singleControlledInletResult.inletNeumannFaces << "\n";
        return 1;
    }

    std::vector<double> linearPressure(
        static_cast<std::size_t>(domain.N()),0.0);
    const std::array<double,3> exactGradient{2.3e5,-1.1e5,0.7e5};
    for(int i=0;i<domain.Nx();++i)
    for(int j=0;j<domain.Ny();++j)
    for(int k=0;k<domain.Nz();++k)
    {
        const int cell=k*plane+j*nx+i;
        const FaceAreaVector center=mesh.cellCenter(i,j,k);
        linearPressure[static_cast<std::size_t>(cell)]=5.0e6
            +exactGradient[0]*center.x+exactGradient[1]*center.y
            +exactGradient[2]*center.z;
    }
    const int gradientCell=2*plane+3*nx+4;
    bool leastSquaresApplied=false;
    const FaceAreaVector reconstructedGradient=
        SectorPressureSchurPreconditioner::reconstructPressureGradient(
            domain,mesh,active,allDirichletInlet,
            linearPressure,gradientCell,true,&leastSquaresApplied);
    const double exactGradientMagnitude=std::sqrt(
        exactGradient[0]*exactGradient[0]
        +exactGradient[1]*exactGradient[1]
        +exactGradient[2]*exactGradient[2]);
    const double leastSquaresGradientRelativeError=std::sqrt(
        (reconstructedGradient.x-exactGradient[0])
            *(reconstructedGradient.x-exactGradient[0])
        +(reconstructedGradient.y-exactGradient[1])
            *(reconstructedGradient.y-exactGradient[1])
        +(reconstructedGradient.z-exactGradient[2])
            *(reconstructedGradient.z-exactGradient[2]))
        /exactGradientMagnitude;
    if(!leastSquaresApplied||leastSquaresGradientRelativeError>1.0e-10)
    {
        std::cerr << "least-squares pressure-gradient contract failed: "
                  << leastSquaresApplied << "/"
                  << leastSquaresGradientRelativeError
                  << " reconstructed=" << reconstructedGradient.x << "/"
                  << reconstructedGradient.y << "/"
                  << reconstructedGradient.z << " exact="
                  << exactGradient[0] << "/" << exactGradient[1] << "/"
                  << exactGradient[2] << "\n";
        return 1;
    }

    const auto result=SectorPressureSchurPreconditioner::build(
        domain,mesh,residual,active,spectral,allDirichletInlet,
        5.0,400,1.0e-10);
    if(!(result.residualRatio<1.0e-8)
        ||!std::isfinite(result.maximumRelativePressure)
        ||!std::isfinite(result.maximumVelocityMach))
    {
        std::cerr << "pressure Schur solve did not converge: ratio="
                  << result.residualRatio << "\n";
        return 1;
    }

    double pressureSourceDot=0.0;
    double maximumPeriodicDifference=0.0;
    double maximumPressure=0.0;
    double maximumVelocity=0.0;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    {
        const int reference=(1*plane+j*nx+i);
        for(int k=1;k<domain.Nz()-1;++k)
        {
            const int cell=k*plane+j*nx+i;
            pressureSourceDot+=result.pressure[static_cast<std::size_t>(cell)]
                *residual[static_cast<std::size_t>(cell)].rho;
            maximumPressure=std::max(maximumPressure,
                std::abs(result.pressure[static_cast<std::size_t>(cell)]));
            maximumPeriodicDifference=std::max(maximumPeriodicDifference,
                std::abs(result.pressure[static_cast<std::size_t>(cell)]
                    -result.pressure[static_cast<std::size_t>(reference)]));
            const auto& velocity=result.velocity[static_cast<std::size_t>(cell)];
            maximumVelocity=std::max(maximumVelocity,std::sqrt(
                velocity[0]*velocity[0]+velocity[1]*velocity[1]
                    +velocity[2]*velocity[2]));
        }
    }
    if(!(pressureSourceDot>0.0)||!(maximumPressure>0.0)
        ||!(maximumVelocity>0.0)
        ||maximumPeriodicDifference>1.0e-10*maximumPressure)
    {
        std::cerr << "pressure Schur sign/periodic contract failed: dot="
                  << pressureSourceDot << " pressure=" << maximumPressure
                  << " velocity=" << maximumVelocity
                  << " periodic=" << maximumPeriodicDifference << "\n";
        return 1;
    }

    StructuredImplicitStencil variableSpectral=spectral;
    std::vector<ConservedVars> variableResidual=residual;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const int cell=k*plane+j*nx+i;
        variableSpectral.diagonal[static_cast<std::size_t>(cell)]=
            2.5e7*std::pow(2.0,static_cast<double>((i+j)%5));
        variableResidual[static_cast<std::size_t>(cell)].rhoU[0]=
            2.0e5*std::sin(0.4*i+0.2*j);
        variableResidual[static_cast<std::size_t>(cell)].rhoU[1]=
            -1.0e5*std::cos(0.3*i-0.5*j);
    }
    const auto variableResult=SectorPressureSchurPreconditioner::build(
        domain,mesh,variableResidual,active,variableSpectral,
        allDirichletInlet,5.0,400,1.0e-10);
    const auto variableLeastSquaresResult=
        SectorPressureSchurPreconditioner::build(
            domain,mesh,variableResidual,active,variableSpectral,
            allDirichletInlet,5.0,400,1.0e-10,true);
    const auto variableNonOrthogonalResult=
        SectorPressureSchurPreconditioner::build(
            domain,mesh,variableResidual,active,variableSpectral,
            allDirichletInlet,5.0,400,1.0e-10,true,2);
    if(!(variableResult.residualRatio<1.0e-8)
        ||!(variableResult.minimumFacePseudoTime>0.0)
        ||!(variableResult.maximumFacePseudoTime
            >4.0*variableResult.minimumFacePseudoTime)
        ||variableResult.maximumFaceReciprocityError>1.0e-12
        ||!(variableResult.maximumVelocityMach>0.0)
        ||!(variableResult.maximumFaceNonOrthogonalityDegrees>20.0)
        ||variableLeastSquaresResult.leastSquaresGradientApplications<=0
        ||variableLeastSquaresResult.leastSquaresGradientApplications
            +variableLeastSquaresResult.leastSquaresGradientFallbacks
                !=(domain.Nx()-2)*(domain.Ny()-2)*(domain.Nz()-2)
        ||variableNonOrthogonalResult.nonOrthogonalCorrections<1
        ||variableNonOrthogonalResult.nonOrthogonalCorrections>2
        ||!std::isfinite(
            variableNonOrthogonalResult.nonOrthogonalResidualRatio)
        ||variableNonOrthogonalResult.nonOrthogonalResidualRatio>1.0e-8
        ||variableNonOrthogonalResult.leastSquaresGradientFallbacks!=0)
    {
        std::cerr << "variable-mobility Schur contract failed: ratio="
                  << variableResult.residualRatio
                  << " dtau=" << variableResult.minimumFacePseudoTime
                  << "/" << variableResult.maximumFacePseudoTime
                  << " reciprocity="
                  << variableResult.maximumFaceReciprocityError
                  << " velocity=" << variableResult.maximumVelocityMach
                  << " nonorth="
                  << variableResult.maximumFaceNonOrthogonalityDegrees
                  << " LS="
                  << variableLeastSquaresResult.leastSquaresGradientApplications
                  << "/"
                  << variableLeastSquaresResult.leastSquaresGradientFallbacks
                  << " nonorth_ratio/corrections="
                  << variableNonOrthogonalResult.nonOrthogonalResidualRatio
                  << "/"
                  << variableNonOrthogonalResult.nonOrthogonalCorrections
                  << "\n";
        return 1;
    }

    std::vector<ConservedVars> rotationalResidual(
        static_cast<std::size_t>(domain.N()),ConservedVars(2));
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const int cell=k*plane+j*nx+i;
        const double theta=mesh.thetaCenter(k);
        const double axial=1.0e5*std::sin(0.5*i);
        const double radial=2.0e5*(0.2+0.05*j);
        const double tangential=-7.5e4*std::cos(0.3*i+0.2*j);
        rotationalResidual[static_cast<std::size_t>(cell)].rho=
            5.0e2*(std::sin(0.6*i)+0.2*std::cos(0.4*j));
        rotationalResidual[static_cast<std::size_t>(cell)].rhoU={
            axial,
            radial*std::cos(theta)-tangential*std::sin(theta),
            radial*std::sin(theta)+tangential*std::cos(theta)};
    }
    const auto rotationalResult=SectorPressureSchurPreconditioner::build(
        domain,mesh,rotationalResidual,active,variableSpectral,
        allDirichletInlet,5.0,400,1.0e-10);
    double rotationalPressureScale=0.0;
    double rotationalVelocityScale=0.0;
    double rotationalPressureError=0.0;
    double rotationalVelocityError=0.0;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    {
        const int reference=1*plane+j*nx+i;
        const double referenceTheta=mesh.thetaCenter(1);
        const auto& referenceVelocity=rotationalResult.velocity[
            static_cast<std::size_t>(reference)];
        const std::array<double,3> referenceCylindrical{
            referenceVelocity[0],
            referenceVelocity[1]*std::cos(referenceTheta)
                +referenceVelocity[2]*std::sin(referenceTheta),
            -referenceVelocity[1]*std::sin(referenceTheta)
                +referenceVelocity[2]*std::cos(referenceTheta)};
        for(int k=1;k<domain.Nz()-1;++k)
        {
            const int cell=k*plane+j*nx+i;
            const double theta=mesh.thetaCenter(k);
            const auto& velocity=rotationalResult.velocity[
                static_cast<std::size_t>(cell)];
            const std::array<double,3> cylindrical{
                velocity[0],
                velocity[1]*std::cos(theta)+velocity[2]*std::sin(theta),
                -velocity[1]*std::sin(theta)+velocity[2]*std::cos(theta)};
            rotationalPressureScale=std::max(rotationalPressureScale,
                std::abs(rotationalResult.pressure[
                    static_cast<std::size_t>(cell)]));
            rotationalPressureError=std::max(rotationalPressureError,
                std::abs(rotationalResult.pressure[
                    static_cast<std::size_t>(cell)]
                    -rotationalResult.pressure[
                        static_cast<std::size_t>(reference)]));
            for(int direction=0;direction<3;++direction)
            {
                rotationalVelocityScale=std::max(rotationalVelocityScale,
                    std::abs(cylindrical[static_cast<std::size_t>(direction)]));
                rotationalVelocityError=std::max(rotationalVelocityError,
                    std::abs(cylindrical[static_cast<std::size_t>(direction)]
                        -referenceCylindrical[
                            static_cast<std::size_t>(direction)]));
            }
        }
    }
    const double rotationalPressureRelativeError=rotationalPressureError
        /std::max(rotationalPressureScale,1.0e-300);
    const double rotationalVelocityRelativeError=rotationalVelocityError
        /std::max(rotationalVelocityScale,1.0e-300);
    if(!(rotationalResult.residualRatio<1.0e-8)
        ||rotationalPressureRelativeError>1.0e-10
        ||rotationalVelocityRelativeError>1.0e-10)
    {
        std::cerr << "rotational-periodic Schur contract failed: ratio="
                  << rotationalResult.residualRatio
                  << " pressure=" << rotationalPressureRelativeError
                  << " velocity=" << rotationalVelocityRelativeError
                  << "\n";
        return 1;
    }

    std::vector<unsigned char> mixedInlet=allDirichletInlet;
    int expectedDirichlet=0,expectedNeumann=0;
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const bool dirichlet=(j+k)%2==0;
        mixedInlet[static_cast<std::size_t>(k*domain.Ny()+j)]=
            dirichlet?1U:0U;
        if(dirichlet) ++expectedDirichlet;
        else ++expectedNeumann;
    }
    const auto mixedResult=SectorPressureSchurPreconditioner::build(
        domain,mesh,residual,active,spectral,mixedInlet,
        5.0,400,1.0e-10);
    const auto mixedControlledResult=
        SectorPressureSchurPreconditioner::build(
            domain,mesh,residual,active,spectral,mixedInlet,
            5.0,400,1.0e-10,false,0,true);
    const int activeInletFaces=(domain.Ny()-2)*(domain.Nz()-2);
    if(!(mixedResult.residualRatio<1.0e-8)
        ||mixedResult.inletPressureDirichletFaces!=expectedDirichlet
        ||mixedResult.inletDirichletFaces!=expectedDirichlet
        ||mixedResult.inletNeumannFaces!=expectedNeumann
        ||mixedControlledResult.inletPressureDirichletFaces
            !=expectedDirichlet
        ||mixedControlledResult.inletDirichletFaces!=0
        ||mixedControlledResult.inletNeumannFaces!=activeInletFaces
        ||result.inletPressureDirichletFaces!=activeInletFaces
        ||result.inletDirichletFaces!=activeInletFaces
        ||result.inletNeumannFaces!=0)
    {
        std::cerr << "mixed-inlet Schur contract failed: ratio="
                  << mixedResult.residualRatio
                  << " pressure-D/continuity-D/N="
                  << mixedControlledResult.inletPressureDirichletFaces
                  << "/" << mixedControlledResult.inletDirichletFaces
                  << "/" << mixedControlledResult.inletNeumannFaces
                  << " expected=" << expectedDirichlet
                  << "/0/" << activeInletFaces << "\n";
        return 1;
    }

    std::cout << "single_D_pressure/velocity_error="
              << singlePressureRelativeError << "/"
              << singleVelocityRelativeError
              << " controlled_pressure/velocity_error="
              << controlledPressureRelativeError << "/"
              << controlledVelocityRelativeError
              << " linear_LS_error="
              << leastSquaresGradientRelativeError
              << " sweeps/ratio=" << result.sweeps << "/"
              << result.residualRatio
              << " max_dp_rel/max_du_mach="
              << result.maximumRelativePressure << "/"
              << result.maximumVelocityMach
              << " periodic=" << maximumPeriodicDifference
              << " variable_sweeps/ratio/dtau/reciprocity="
              << variableResult.sweeps << "/"
              << variableResult.residualRatio << "/"
              << variableResult.minimumFacePseudoTime << "/"
              << variableResult.maximumFacePseudoTime << "/"
              << variableResult.maximumFaceReciprocityError
              << " nonorth/LS/fallback="
              << variableResult.maximumFaceNonOrthogonalityDegrees << "/"
              << variableLeastSquaresResult.leastSquaresGradientApplications
              << "/"
              << variableLeastSquaresResult.leastSquaresGradientFallbacks
              << " nonorth_ratio/corrections="
              << variableNonOrthogonalResult.nonOrthogonalResidualRatio
              << "/"
              << variableNonOrthogonalResult.nonOrthogonalCorrections
              << " rotational_ratio/pressure/velocity="
              << rotationalResult.residualRatio << "/"
              << rotationalPressureRelativeError << "/"
              << rotationalVelocityRelativeError
              << " mixed_ratio/D/N=" << mixedResult.residualRatio << "/"
              << mixedControlledResult.inletPressureDirichletFaces << "/"
              << mixedControlledResult.inletDirichletFaces << "/"
              << mixedControlledResult.inletNeumannFaces
              << "\nSECTOR_PRESSURE_SCHUR_PASS\n";
    return 0;
}
