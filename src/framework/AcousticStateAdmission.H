#ifndef AcousticStateAdmission_H
#define AcousticStateAdmission_H

#include "EquilibriumThermoModel.H"
#include "ModelCapabilities.H"

namespace Foam
{

inline bool solverConsumerAdmitsAcousticState(
    const thermo::EquilibriumState& state) noexcept
{
    if (!state.valid()) return false;

    switch (state.acoustic.closure)
    {
        case thermo::AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE:
            return true;
        case thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM:
#ifdef RFAM_ALLOW_DIAGNOSTIC_HEM_ACOUSTICS
            // EXPLICIT diagnostic override (2026-08-04 audit): the HEM
            // acoustic closure is a self-consistent placeholder, NOT a
            // converged derivative.  Only diagnostic solvers (e.g.
            // solver_mayer_o2) compile this in; production builds reject
            // the state below.
            return state.acoustic.evidence
                != thermo::AcousticEvidence::UNVERIFIED
                || state.valid();
#else
            return ModelCapabilities::phaseAwareSoundSpeed
                && state.acoustic.evidence
                    == thermo::AcousticEvidence::CONVERGED_RHOE_ISENTROPE;
#endif
        case thermo::AcousticClosure::UNKNOWN:
            return false;
    }
    return false;
}

} // End namespace Foam

#endif
