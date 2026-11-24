#ifndef RFAM_TETRA_MOLECULAR_ADMISSION_H
#define RFAM_TETRA_MOLECULAR_ADMISSION_H

#include "EquilibriumThermoModel.H"
#include <cmath>

namespace Foam
{
// This development path has no admitted two-phase transport closure.
inline bool tetraAdmitsMolecularTransport(
    thermo::EquilibriumRegion region,double flashQuality)
{
    return (region==thermo::EquilibriumRegion::SINGLE_LIQUID
            ||region==thermo::EquilibriumRegion::SINGLE_VAPOR
            ||region==thermo::EquilibriumRegion::SUPERCRITICAL
            ||region==thermo::EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE)
        &&!(std::isfinite(flashQuality)&&flashQuality>=0.&&flashQuality<=1.);
}
}
#endif
