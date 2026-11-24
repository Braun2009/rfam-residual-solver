#ifndef ModelCapabilities_H
#define ModelCapabilities_H

namespace Foam
{

// Production-readiness manifest. A flag is true only after its implementation
// and validation matrix are both present. The reacting sector gate depends on
// the conjunction, so no single local switch can bypass an unfinished closure.
//
// ---------------------------------------------------------------------------
// Where RANS and CHT enter this manifest, and why they are absent today.
//
// The two closures exist (src/turbulence, src/cht) but neither is reachable
// from the tetrahedral production driver: solver_3d_tetra.cpp includes neither
// header, and its file comment states "No combustion, diffusion or turbulence
// model is enabled here." A flag here is a claim about the SOLVER, not about
// the existence of a header, so adding ransClosure/chtClosure before the
// driver can consume them would be exactly the kind of unearned capability
// claim this manifest was created to prevent.
//
// The intended decomposition, so the flags can be added one at a time:
//
//  RANS (turbulence)
//    ransClosure                 k-omega SST + WALE available and tested in
//                                isolation (largely true, sector path only)
//    ransWallDistance            a real wall-distance FIELD, not the
//                                first-cell half-height placeholder. Required
//                                before any F1/F2-dependent quantity is quoted.
//    ransRealFluidCoupling       density-consistent (non-Boussinesq) stress
//                                closure for the cryogenic interface
//    ransTetraIntegration        operator wired into TetraFiniteVolume's
//                                molecular-stress path with the same step
//                                transaction as U/W
//
//  CHT (conjugate heat transfer)
//    chtConductionSolid          SolidEnergyEqn / CHTCoupler validated
//    chtMultiMaterialWall        per-material rho*cp/k, not one averaged set
//    chtWallFluxCoupling         fluid-side q at the wall from a validated
//                                closure, with the FixedFlux sign convention
//                                enforced by construction
//    chtMonolithicResidual       one global nonlinear residual across
//                                fluid+solid+coolant, re-assembled per outer
//                                Newton step (the CHTCoupler system is
//                                monolithic only for FIXED wall data)
//
// Until those are true, turbulent or CHT numbers must be reported as
// diagnostics, and CHTFluidProductionReady()/ransTetraProductionReady() below
// must remain false. Keep this list in sync with the source comments in
// src/turbulence/KOmegaSST.H and src/cht/CHTCoupler.H.
// ---------------------------------------------------------------------------
struct ModelCapabilities
{
    static constexpr bool twoPhaseFlash = false;
    static constexpr bool phaseAwareEnergyInversion = false;
    static constexpr bool latentHeatClosure = false;
    // Two-phase acoustics are a DIAGNOSTIC PLACEHOLDER (self-consistent
    // ideal-gas-like closure, NOT a converged HEM derivative).  It must
    // never gate production fluxes/CFL.  Reverted 2026-08-04 per audit:
    // the placeholder sound speed was previously admitted to AUSM+/CFL.
    static constexpr bool phaseAwareSoundSpeed = false;
    static constexpr bool phaseAwareBoundaryStates = false;
    static constexpr bool phaseAwareCharacteristicBoundary = false;
    static constexpr bool denseCryogenicTransport = false;
    static constexpr bool phaseAwareChemistry = false;
    static constexpr bool condensedCarbonTransport = false;
    static constexpr bool condensedCarbonKinetics = false;
    static constexpr bool condensedCarbonWallDeposition = false;

    // Turbulence (RANS/LES). All false: the closures are implemented but the
    // tetrahedral production driver does not consume them. See the block
    // comment above for what each flag would require.
    static constexpr bool ransClosure = false;
    static constexpr bool ransWallDistance = false;
    static constexpr bool ransRealFluidCoupling = false;
    static constexpr bool ransTetraIntegration = false;

    // Conjugate heat transfer. All false: conduction-only 1D coupling exists;
    // a monolithic fluid-solid-coolant residual does not.
    static constexpr bool chtConductionSolid = false;
    static constexpr bool chtMultiMaterialWall = false;
    static constexpr bool chtWallFluxCoupling = false;
    static constexpr bool chtMonolithicResidual = false;

    static constexpr bool nonReactingSectorProductionReady()
    {
        return twoPhaseFlash
            && phaseAwareEnergyInversion
            && latentHeatClosure
            && phaseAwareSoundSpeed
            && phaseAwareBoundaryStates
            && phaseAwareCharacteristicBoundary
            && denseCryogenicTransport;
    }

    static constexpr bool reactingSectorProductionReady()
    {
        return nonReactingSectorProductionReady()
            && phaseAwareChemistry;
    }

    static constexpr bool condensedCarbonProductionReady()
    {
        return condensedCarbonTransport
            && condensedCarbonKinetics
            && condensedCarbonWallDeposition;
    }

    // Turbulence production readiness. Note this is deliberately NOT part of
    // nonReactingSectorProductionReady(): a solver must be able to run
    // laminar-with-provenance before a turbulence model can be layered on it.
    static constexpr bool ransTetraProductionReady()
    {
        return ransClosure
            && ransWallDistance
            && ransRealFluidCoupling
            && ransTetraIntegration;
    }

    // CHT production readiness. The monolithic residual is the last gate: the
    // first three flags describe a solver that can exchange wall data, which
    // is necessary but not the strong coupling the architecture targets.
    static constexpr bool chtProductionReady()
    {
        return chtConductionSolid
            && chtMultiMaterialWall
            && chtWallFluxSignCoupled()
            && chtMonolithicResidual;
    }

    // Kept separate so the sign/closure pair can be asserted together: a wall
    // flux without the correct into-solid convention is a wrong answer that
    // passes every magnitude-only check.
    static constexpr bool chtWallFluxSignCoupled()
    {
        return chtWallFluxCoupling;
    }
};

} // End namespace Foam

#endif
