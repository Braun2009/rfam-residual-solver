#ifndef CoolPropThermoModel_H
#define CoolPropThermoModel_H
#include "EquilibriumThermoModel.H"
namespace Foam { namespace thermo {
// Native CoolProp HEOS. Pure-fluid rho-e flash uses CoolProp directly;
// mixture rho-e inversion closes against equilibrium PT flash (not imposed
// gas/liquid phases). HEM sound speed comes from independently converged
// equilibrium isentropes. Unsupported flash states return failure.
// Sensible h/e subtract mass-weighted ideal h(298.15 K) from the HEOS values.
std::shared_ptr<const EquilibriumThermoModel> makeCoolPropThermoModel();
bool coolPropThermoAvailable();
}}
#endif
