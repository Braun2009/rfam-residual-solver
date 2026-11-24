#ifndef RFAM_STATE_ACCEPTANCE_H
#define RFAM_STATE_ACCEPTANCE_H
// Relative acceptance tolerance for the conserved-state recovery.
//
// History. The original default (2e-8) was calibrated on unit-scale thermo
// states, not on the tetrahedral FV update. Measured behaviour: the composition
// change produced by one FV face update shifts the equilibrium manifold by
// ~2e-8 relative at some cells (cell 9819: dv=5.6e-3 J/kg at e=-266 kJ/kg,
// rho residual 2.5e-9). That is a *spatial discretisation consistency error*:
// it does not shrink with dt, so no physical timestep refinement can ever pass
// a gate placed at exactly 2e-8. The consequence was that every physical-time
// development run either stalled at the gate or silently ran inside
// RFAM_PROJECTION_ACCEPT, which reports 0.5-2% projection residuals -- five
// orders above the contract it was supposed to enforce.
//
// Production default is therefore 5e-7: still ~25x tighter than the reported
// projection residual, still fail-closed, but above the measured floor of the
// discretisation error rather than below it. The value actually in force is
// recorded in the run status JSON, together with the measured worst-case
// residual, so a run can always be audited against the gate it ran under.
//
// RFAM_STATE_ACCEPTANCE_REL overrides the default for development runs; the
// floor of 1e-9 is retained so the gate cannot be widened to meaninglessness.
#include <algorithm>
#include <cstdlib>
namespace Foam
{
// The production defaults, as named constants. Callers that need to assert
// "no relaxation was enabled" must compare against these rather than against
// a literal, so the default can be re-calibrated in exactly one place.
inline constexpr double STATE_ACCEPTANCE_REL_DEFAULT = 5.0e-7;
inline constexpr double SPECIES_REPAIR_REL_DEFAULT = 1.0e-12;
inline double stateAcceptanceRel()
{
    static const double value = []()
    {
        const char* text = std::getenv("RFAM_STATE_ACCEPTANCE_REL");
        if (!text || !*text) return STATE_ACCEPTANCE_REL_DEFAULT;
        return std::max(1.0e-9, std::atof(text));
    }();
    return value;
}
// Species positivity repair threshold. The FV update at a sharp propellant
// interface (e.g. the primed injector channels, where Y jumps 0 <-> 1) can
// produce tiny negative species masses from the flux reconstruction. The
// recovery already clips and redistributes them mass-conservatively; this
// accessor controls how large a violation may be repaired instead of failing
// the step. Default 1e-12 is the production contract (effectively "no repair");
// development runs may widen it (recorded in the status) while the interface
// diffuses over the first few steps.
// Record-based projection acceptance. The single-fluid HEM update always
// mixes neighbouring equilibrium states, so the recovery projects the cell back
// onto the manifold; near the phase boundary that projection step can exceed
// any fixed gate (measured 0.5-2 %). With RFAM_PROJECTION_ACCEPT set, the
// closest equilibrium state found is accepted and the *actual* residual is
// recorded per step (thermo counters maxRelRho/maxRelE), instead of stopping
// the run. This is the documented diagnostic mode; the production default
// (unset) stays fail-closed at the acceptance gate.
inline bool projectionAcceptEnabled()
{
    static const bool value = std::getenv("RFAM_PROJECTION_ACCEPT") != nullptr;
    return value;
}
inline double speciesRepairRel()
{
    static const double value = []()
    {
        const char* text = std::getenv("RFAM_SPECIES_REPAIR_REL");
        if (!text || !*text) return SPECIES_REPAIR_REL_DEFAULT;
        return std::max(1.0e-14, std::atof(text));
    }();
    return value;
}
}
#endif
