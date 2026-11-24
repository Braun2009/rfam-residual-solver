// Copyright 2026 Braun
// SPDX-License-Identifier: Apache-2.0
#ifndef Residual_H
#define Residual_H

/*---------------------------------------------------------------------------*\
    Residual - the dt-free (rate) form of a conserved update.

    SIGN CONVENTION (read this before adding an operator)
    ----------------------------------------------------
    R is the spatial term on the LEFT-HAND SIDE of the conservation law:

        dU/dt + R(U) = 0,       R = sum_i R_i

    so the update is

        Unew = U - dt * R.

    Therefore, for an operator whose additive contribution is c = dU/dt:
        flux divergence operator :  R += +div(F)      (c = -div F)
        volumetric source        :  R += -S           (c = +S)

    This is the sign that makes `(M/dt + dR/dU) dU = -R` and the preconditioned
    form `P dU/dtau = -R` read the standard way. It is deliberately NOT the same
    sign as `apply()`, and test_residual_integrator catches a sign slip.

    Every operator contributes its own R_i. No operator assigns another
    operator's contribution, so the chain is order-independent up to
    floating-point summation order, and dt is applied exactly once, by the
    integrator:

        Unew = U - dt * R

    RELATION TO THE ADDITIVE PATH
    -----------------------------
    The legacy additive path (`Operator::apply`) keeps its update arithmetic.
    Both paths use the same face-flux routines and the same halo treatment.

        apply:     Unew = U - dt/dx * (Fr - Fl)
        residual:  R   += (Fr - Fl)/dx

    tests/test_residual_integrator.cpp is the gate that keeps the two paths
    agreeing to round-off, and that proves order independence on the residual
    path where the additive path has a hard ordering contract.
\*---------------------------------------------------------------------------*/

#include "ConservedVars.H"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace Foam
{

// Per-cell spatial balance div(F) - S. Units: U per second.
using ResidualField = std::vector<ConservedVars>;

inline void zeroResidual(ResidualField& R, int cells, int nSpecies)
{
    R.assign(static_cast<std::size_t>(cells), ConservedVars(nSpecies));
}

// target += scale * source  (component-wise, species-layout checked by caller)
inline void addScaledRate(
    ConservedVars& target, const ConservedVars& source, double scale)
{
    target.rho  += scale*source.rho;
    target.rhoE += scale*source.rhoE;
    for (std::size_t k = 0; k < target.rhoU.size(); ++k)
        target.rhoU[k] += scale*source.rhoU[k];
    for (std::size_t k = 0; k < target.rhoY.size(); ++k)
        target.rhoY[k] += scale*source.rhoY[k];
    target.rhoZ += scale*source.rhoZ;
    target.rhoC += scale*source.rhoC;
}

// target -= dt * rate
inline void subtractScaledRate(
    ConservedVars& target, const ConservedVars& rate, double dt)
{
    target.rho  -= dt*rate.rho;
    target.rhoE -= dt*rate.rhoE;
    for (std::size_t k = 0; k < target.rhoU.size(); ++k)
        target.rhoU[k] -= dt*rate.rhoU[k];
    for (std::size_t k = 0; k < target.rhoY.size(); ++k)
        target.rhoY[k] -= dt*rate.rhoY[k];
    target.rhoZ -= dt*rate.rhoZ;
    target.rhoC -= dt*rate.rhoC;
}

// Add the divergence using face differences before scaling. This preserves
// exact cancellation for equal face fluxes even when R already has a source.
inline void addFluxDivergence(ConservedVars& r, const ConservedVars& right,
                              const ConservedVars& left, double invWidth)
{
    r.rho += invWidth*(right.rho-left.rho);
    r.rhoE += invWidth*(right.rhoE-left.rhoE);
    for (std::size_t k = 0; k < r.rhoU.size(); ++k)
        r.rhoU[k] += invWidth*(right.rhoU[k]-left.rhoU[k]);
    for (std::size_t k = 0; k < r.rhoY.size(); ++k)
        r.rhoY[k] += invWidth*(right.rhoY[k]-left.rhoY[k]);
    r.rhoZ += invWidth*(right.rhoZ-left.rhoZ);
    r.rhoC += invWidth*(right.rhoC-left.rhoC);
}

} // namespace Foam

#endif
