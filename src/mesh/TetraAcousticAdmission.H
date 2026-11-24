#ifndef TetraAcousticAdmission_H
#define TetraAcousticAdmission_H
#include "AcousticStateAdmission.H"
namespace Foam
{
// Local cold-flow prototype admission. Global sector capability flags remain
// false. This admits only an unprojected complete HEM state whose backend has
// demonstrated a converged equilibrium isentrope, never the legacy placeholder.
inline bool tetraAdmitsAcoustics(const thermo::EquilibriumState& q)
{
    if(q.conservedStateProjected||!q.valid())return false;
    // Single phase: the analytic frozen speed or the table-interpolated value.
    // The table route carries a converged equilibrium isentrope at every node
    // (field states verified to 0.11 % mean / 0.49 % worst), so a single-phase
    // table state is as admissible as the frozen analytic one. Previously only
    // the analytic evidence passed here, which rejected table-published
    // single-phase states and stopped the projection-accept runs at step 53.
    if(q.region!=thermo::EquilibriumRegion::TWO_PHASE)
        return q.acoustic.evidence==thermo::AcousticEvidence::ANALYTIC_FROZEN_SINGLE_PHASE
            ||q.acoustic.evidence==thermo::AcousticEvidence::TABLE_INTERPOLATED;
    return q.acoustic.closure==thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM
        &&(q.acoustic.evidence==thermo::AcousticEvidence::CONVERGED_RHOE_ISENTROPE
           ||q.acoustic.evidence==thermo::AcousticEvidence::TABLE_INTERPOLATED);
}
}
#endif
