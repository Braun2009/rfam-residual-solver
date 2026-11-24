#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "../src/framework/SectorStructuredPressureSchurPreconditioner.H"

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
        std::cerr << "structured pressure Schur Domain initialization failed\n";
        return 1;
    }

    std::vector<unsigned char> active(
        static_cast<std::size_t>(domain.N()),0U);
    StructuredImplicitStencil stencil;
    stencil.ensureSize(static_cast<std::size_t>(domain.N()));
    std::vector<unsigned char> inletPressureDirichlet(
        static_cast<std::size_t>(domain.Ny()*domain.Nz()),0U);
    const int nx=domain.Nx(),ny=domain.Ny();
    const int plane=nx*ny;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const int cell=index(domain,i,j,k);
        active[static_cast<std::size_t>(cell)]=1U;
        stencil.diagonal[static_cast<std::size_t>(cell)]=
            8.0+0.04*i+0.03*j;
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
        if(i==1)
            inletPressureDirichlet[static_cast<std::size_t>(
                k*domain.Ny()+j)]=(j%3)!=0?1U:0U;
    }

    std::vector<double> exactPressure(
        static_cast<std::size_t>(domain.N()),0.0);
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const int cell=index(domain,i,j,k);
        exactPressure[static_cast<std::size_t>(cell)]=1.0e3*(
            0.3+0.05*i-0.02*j+0.01*i*j);
    }

    const double pseudoCfl=2.5;
    SectorStructuredPressureSchurPreconditioner::Audit operatorAudit;
    const std::vector<double> rightHandSide=
        SectorStructuredPressureSchurPreconditioner::applyPressureOperator(
            domain,mesh,active,stencil,inletPressureDirichlet,
            pseudoCfl,exactPressure,operatorAudit);
    std::vector<ConservedVars> residual(
        static_cast<std::size_t>(domain.N()),ConservedVars(2));
    for(int cell=0;cell<domain.N();++cell)
        residual[static_cast<std::size_t>(cell)].rho=
            rightHandSide[static_cast<std::size_t>(cell)];

    const auto solved=SectorStructuredPressureSchurPreconditioner::build(
        domain,mesh,residual,active,stencil,inletPressureDirichlet,
        pseudoCfl,120,1.0e-10);

    double pressureErrorSquared=0.0;
    double pressureSquared=0.0;
    for(int cell=0;cell<domain.N();++cell)
    {
        if(active[static_cast<std::size_t>(cell)]==0U) continue;
        const double expected=exactPressure[static_cast<std::size_t>(cell)];
        const double error=solved.pressure[static_cast<std::size_t>(cell)]
            -expected;
        pressureErrorSquared+=error*error;
        pressureSquared+=expected*expected;
    }
    const double pressureRelativeError=std::sqrt(
        pressureErrorSquared/std::max(pressureSquared,1.0e-300));

    // Independently reconstruct the expected final momentum correction for
    // F_m=0: delta_m=-M_m^-1 G(delta_p).
    using Vector=SectorStructuredMomentumInverse::Vector;
    std::vector<Vector> pressureGradient(
        static_cast<std::size_t>(domain.N()),Vector{0.0,0.0,0.0});
    for(int cell=0;cell<domain.N();++cell)
    {
        if(active[static_cast<std::size_t>(cell)]==0U) continue;
        bool leastSquaresApplied=false;
        const FaceAreaVector gradient=
            SectorPressureSchurPreconditioner::reconstructPressureGradient(
                domain,mesh,active,inletPressureDirichlet,
                exactPressure,cell,true,&leastSquaresApplied);
        if(!leastSquaresApplied)
        {
            std::cerr << "structured Schur manufactured gradient fell back\n";
            return 1;
        }
        pressureGradient[static_cast<std::size_t>(cell)]={
            gradient.x,gradient.y,gradient.z};
    }
    const auto expectedMomentum=SectorStructuredMomentumInverse::solve(
        domain,mesh,active,stencil,pseudoCfl,pressureGradient,
        24,1.0e-14,24);
    double velocityErrorSquared=0.0;
    double velocitySquared=0.0;
    double maximumPeriodicCylindricalError=0.0;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    {
        const int lower=index(domain,i,j,1);
        const int upper=index(domain,i,j,domain.Nz()-2);
        const double lowerTheta=mesh.thetaCenter(1);
        const double upperTheta=mesh.thetaCenter(domain.Nz()-2);
        const auto& lowerVelocity=solved.velocity[
            static_cast<std::size_t>(lower)];
        const auto& upperVelocity=solved.velocity[
            static_cast<std::size_t>(upper)];
        const std::array<double,3> lowerCylindrical{
            lowerVelocity[0],
            lowerVelocity[1]*std::cos(lowerTheta)
                +lowerVelocity[2]*std::sin(lowerTheta),
            -lowerVelocity[1]*std::sin(lowerTheta)
                +lowerVelocity[2]*std::cos(lowerTheta)};
        const std::array<double,3> upperCylindrical{
            upperVelocity[0],
            upperVelocity[1]*std::cos(upperTheta)
                +upperVelocity[2]*std::sin(upperTheta),
            -upperVelocity[1]*std::sin(upperTheta)
                +upperVelocity[2]*std::cos(upperTheta)};
        for(int component=0;component<3;++component)
            maximumPeriodicCylindricalError=std::max(
                maximumPeriodicCylindricalError,
                std::abs(lowerCylindrical[
                    static_cast<std::size_t>(component)]
                    -upperCylindrical[
                        static_cast<std::size_t>(component)]));
    }
    for(int cell=0;cell<domain.N();++cell)
    {
        if(active[static_cast<std::size_t>(cell)]==0U) continue;
        const double rho=domain.Q(cell).rho;
        const Vector expected{
            -expectedMomentum.correction[static_cast<std::size_t>(cell)][0]
                /rho,
            -expectedMomentum.correction[static_cast<std::size_t>(cell)][1]
                /rho,
            -expectedMomentum.correction[static_cast<std::size_t>(cell)][2]
                /rho};
        const Vector& actual=solved.velocity[static_cast<std::size_t>(cell)];
        Vector error{actual[0]-expected[0],actual[1]-expected[1],
            actual[2]-expected[2]};
        velocityErrorSquared+=vectorMagnitude(error)*vectorMagnitude(error);
        velocitySquared+=vectorMagnitude(expected)*vectorMagnitude(expected);
    }
    const double velocityRelativeError=std::sqrt(
        velocityErrorSquared/std::max(velocitySquared,1.0e-300));
    int expectedPressureDiffusionFaces=0;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        expectedPressureDiffusionFaces+=2; // periodic K-/K+ neighbors
        if(i>1) ++expectedPressureDiffusionFaces;
        if(i<domain.Nx()-2) ++expectedPressureDiffusionFaces;
        if(j>1) ++expectedPressureDiffusionFaces;
        if(j<domain.Ny()-2) ++expectedPressureDiffusionFaces;
    }

    if(solved.residualRatio>1.0e-10
        ||pressureRelativeError>1.0e-9
        ||velocityRelativeError>1.0e-9
        ||maximumPeriodicCylindricalError>1.0e-9
        ||solved.audit.leastSquaresFallbacks!=0
        ||solved.audit.momentumSolves<=0
        ||solved.audit.maximumMomentumResidualRatio>1.0e-8
        ||solved.audit.acousticPreconditionerSolves<=0
        ||!std::isfinite(
            solved.audit.maximumAcousticPreconditionerResidualRatio)
        ||solved.audit.maximumAcousticPreconditionerResidualRatio>=1.0
        ||solved.audit.pressureDiffusionFaces
            !=expectedPressureDiffusionFaces
        ||!std::isfinite(solved.audit.minimumPressureDiffusionMobility)
        ||!(solved.audit.minimumPressureDiffusionMobility>0.0)
        ||!(solved.audit.maximumPressureDiffusionMobility
            >=solved.audit.minimumPressureDiffusionMobility)
        ||!(solved.audit.maximumPressureDiffusionCoefficient>0.0)
        ||solved.audit.maximumPressureDiffusionReciprocityError>1.0e-12)
    {
        std::cerr << "structured pressure Schur failed: krylov/ratio/p/u/periodic="
                  << solved.krylovIterations << "/"
                  << solved.residualRatio << "/"
                  << pressureRelativeError << "/"
                  << velocityRelativeError << "/"
                  << maximumPeriodicCylindricalError
                  << " momentum=" << solved.audit.momentumSolves << "/"
                  << solved.audit.maximumMomentumResidualRatio
                  << " acoustic="
                  << solved.audit.acousticPreconditionerSolves << "/"
                  << solved.audit.maximumAcousticPreconditionerResidualRatio
                  << " LSfallback=" << solved.audit.leastSquaresFallbacks
                  << " pressure_diffusion/faces/min/max/coefficient/reciprocity="
                  << solved.audit.pressureDiffusionFaces << "/"
                  << solved.audit.minimumPressureDiffusionMobility << "/"
                  << solved.audit.maximumPressureDiffusionMobility << "/"
                  << solved.audit.maximumPressureDiffusionCoefficient << "/"
                  << solved.audit.maximumPressureDiffusionReciprocityError
                  << "\n";
        return 1;
    }

    std::cout << "krylov/ratio/p/u/periodic="
              << solved.krylovIterations << "/"
              << solved.residualRatio << "/"
              << pressureRelativeError << "/"
              << velocityRelativeError << "/"
              << maximumPeriodicCylindricalError
              << " momentum_solves/sweeps/maxratio="
              << solved.audit.momentumSolves << "/"
              << solved.audit.momentumSweeps << "/"
              << solved.audit.maximumMomentumResidualRatio
              << " acoustic_solves/sweeps/maxratio="
              << solved.audit.acousticPreconditionerSolves << "/"
              << solved.audit.acousticPreconditionerSweeps << "/"
              << solved.audit.maximumAcousticPreconditionerResidualRatio
              << " LS/fallback=" << solved.audit.leastSquaresApplications
              << "/" << solved.audit.leastSquaresFallbacks
              << " pressure_diffusion/faces/min/max/coefficient/reciprocity="
              << solved.audit.pressureDiffusionFaces << "/"
              << solved.audit.minimumPressureDiffusionMobility << "/"
              << solved.audit.maximumPressureDiffusionMobility << "/"
              << solved.audit.maximumPressureDiffusionCoefficient << "/"
              << solved.audit.maximumPressureDiffusionReciprocityError
              << "\nSECTOR_STRUCTURED_PRESSURE_SCHUR_PASS\n";
    return 0;
}
