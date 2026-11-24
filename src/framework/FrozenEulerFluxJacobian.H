#ifndef FrozenEulerFluxJacobian_H
#define FrozenEulerFluxJacobian_H

/*---------------------------------------------------------------------------*\
    FrozenEulerFluxJacobian - analytic directional Euler flux linearisation.

    This helper applies dF_n/dU to a conserved perturbation at fixed species
    thermodynamic derivatives. It consumes the already committed equilibrium
    payload and performs no finite-difference residual or EOS evaluations.

    Pressure derivatives are transformed from (rho,T) to (rho,e):

      p_e|rho = p_T|rho / e_T|rho
      p_rho|e = p_rho|T - p_T|rho e_rho|T / e_T|rho.

    Composition perturbations are transported, but their direct contribution
    to pressure is deliberately frozen. This is a robust approximate block for
    LU-SGS, not a claim of a complete reacting-mixture Jacobian.
\*---------------------------------------------------------------------------*/

#include "Domain.H"
#include <algorithm>
#include <array>
#include <cmath>

namespace Foam
{

class FrozenEulerFluxJacobian
{
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

public:
    static bool apply(
        const PrimitiveVars& primitive,
        const thermo::EquilibriumState& equilibrium,
        const ConservedVars& perturbation,
        const std::array<double,3>& unitNormal,
        ConservedVars& product)
    {
        if(!equilibrium.valid()||!equilibrium.thermodynamics.has_value()
            ||!equilibrium.thermodynamics->valid()
            ||perturbation.rhoU.size()!=3
            ||perturbation.rhoY.size()!=primitive.Y.size()
            ||product.rhoU.size()!=3
            ||product.rhoY.size()!=primitive.Y.size())
            return false;
        const auto& thermo=*equilibrium.thermodynamics;
        const double rho=equilibrium.rho;
        if(!std::isfinite(rho)||rho<=0.0||thermo.dedT_rho<=0.0)
            return false;

        const double normalVelocity=
            primitive.U[0]*unitNormal[0]
            +primitive.U[1]*unitNormal[1]
            +primitive.U[2]*unitNormal[2];
        const double momentumNormal=
            perturbation.rhoU[0]*unitNormal[0]
            +perturbation.rhoU[1]*unitNormal[1]
            +perturbation.rhoU[2]*unitNormal[2];
        const double velocitySquared=
            primitive.U[0]*primitive.U[0]
            +primitive.U[1]*primitive.U[1]
            +primitive.U[2]*primitive.U[2];
        double velocityMomentum=0.0;
        for(int direction=0;direction<3;++direction)
            velocityMomentum+=primitive.U[direction]
                *perturbation.rhoU[static_cast<std::size_t>(direction)];
        const double deltaNormalVelocity=(momentumNormal
            -normalVelocity*perturbation.rho)/rho;
        const double deltaInternalEnergy=(perturbation.rhoE
            -velocityMomentum
            +(0.5*velocitySquared-equilibrium.e)*perturbation.rho)/rho;
        const double pressureEnergyDerivative=
            thermo.dpdT_rho/thermo.dedT_rho;
        const double pressureDensityDerivative=thermo.dpdrho_T
            -thermo.dpdT_rho*thermo.dedrho_T/thermo.dedT_rho;
        const double deltaPressure=pressureDensityDerivative*perturbation.rho
            +pressureEnergyDerivative*deltaInternalEnergy;
        if(!std::isfinite(deltaPressure)
            ||!std::isfinite(deltaNormalVelocity)) return false;

        product.rho=momentumNormal;
        for(int direction=0;direction<3;++direction)
            product.rhoU[static_cast<std::size_t>(direction)]=
                normalVelocity*perturbation.rhoU[
                    static_cast<std::size_t>(direction)]
                +rho*primitive.U[direction]*deltaNormalVelocity
                +unitNormal[static_cast<std::size_t>(direction)]
                    *deltaPressure;
        const double totalEnergyDensity=rho
            *(equilibrium.e+0.5*velocitySquared);
        product.rhoE=normalVelocity*(perturbation.rhoE+deltaPressure)
            +(totalEnergyDensity+equilibrium.p)*deltaNormalVelocity;
        for(std::size_t species=0;species<product.rhoY.size();++species)
            product.rhoY[species]=normalVelocity*perturbation.rhoY[species]
                +rho*primitive.Y[species]*deltaNormalVelocity;
        product.rhoZ=normalVelocity*perturbation.rhoZ
            +rho*primitive.Z*deltaNormalVelocity;
        product.rhoC=normalVelocity*perturbation.rhoC
            +rho*primitive.C*deltaNormalVelocity;

        auto finite=[](const ConservedVars& value)
        {
            if(!std::isfinite(value.rho)||!std::isfinite(value.rhoE)
                ||!std::isfinite(value.rhoZ)||!std::isfinite(value.rhoC))
                return false;
            for(double component:value.rhoU)
                if(!std::isfinite(component)) return false;
            for(double component:value.rhoY)
                if(!std::isfinite(component)) return false;
            return true;
        };
        return finite(product);
    }

    // Apply |A_n| without numerical eigensystems or EOS differencing. The
    // frozen Euler block has eigenvalues u_n-a, u_n, and u_n+a, so |A_n| is
    // the quadratic interpolation polynomial through those three values.
    static bool applyAbsolute(
        const PrimitiveVars& primitive,
        const thermo::EquilibriumState& equilibrium,
        const ConservedVars& perturbation,
        const std::array<double,3>& unitNormal,
        ConservedVars& product)
    {
        const double normalVelocity=
            primitive.U[0]*unitNormal[0]
            +primitive.U[1]*unitNormal[1]
            +primitive.U[2]*unitNormal[2];
        const double sound=equilibrium.acoustic.soundSpeed();
        if(!std::isfinite(normalVelocity)||!std::isfinite(sound)
            ||sound<=0.0) return false;
        ConservedVars first(perturbation.rhoY.size());
        ConservedVars second(perturbation.rhoY.size());
        if(!apply(primitive,equilibrium,perturbation,unitNormal,first)
            ||!apply(primitive,equilibrium,first,unitNormal,second))
            return false;

        const double eigenvalue[3]={
            normalVelocity-sound,normalVelocity,normalVelocity+sound};
        double coefficient[3]={0.0,0.0,0.0};
        for(int index=0;index<3;++index)
        {
            const int firstOther=(index+1)%3;
            const int secondOther=(index+2)%3;
            const double denominator=
                (eigenvalue[index]-eigenvalue[firstOther])
                *(eigenvalue[index]-eigenvalue[secondOther]);
            if(!std::isfinite(denominator)||std::abs(denominator)<1.0e-30)
                return false;
            const double weight=std::abs(eigenvalue[index])/denominator;
            coefficient[2]+=weight;
            coefficient[1]-=weight*(eigenvalue[firstOther]
                +eigenvalue[secondOther]);
            coefficient[0]+=weight*eigenvalue[firstOther]
                *eigenvalue[secondOther];
        }
        product=ConservedVars(perturbation.rhoY.size());
        addScaled(product,perturbation,coefficient[0]);
        addScaled(product,first,coefficient[1]);
        addScaled(product,second,coefficient[2]);
        auto finite=[](const ConservedVars& value)
        {
            if(!std::isfinite(value.rho)||!std::isfinite(value.rhoE)
                ||!std::isfinite(value.rhoZ)||!std::isfinite(value.rhoC))
                return false;
            for(double component:value.rhoU)
                if(!std::isfinite(component)) return false;
            for(double component:value.rhoY)
                if(!std::isfinite(component)) return false;
            return true;
        };
        return finite(product);
    }
};

} // End namespace Foam

#endif
