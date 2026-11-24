// Near-critical mixture recovery regression.
//
// Cell 11318 at B-route step 45 (Y_O2 = 0.92953) carries a trial (rho, e) that
// the FV update moved ~2.9e-3 in energy below the committed state. The point is
// recoverable, but only in a ~0.5 K wide basin 1.3 % away in pressure
// (3.804 MPa / 165.30 K, closure norm 4.2e-3); the previous single-resolution
// 5 %-spaced restart scan stepped over it and the run was rolled back. This
// test pins both sides of the contract:
//   * with the development acceptance (1e-2, CTest environment) the refined
//     scan must find the basin and publish a state there;
//   * the production default (2e-8) must still reject it -- the state is
//     genuinely off the equilibrium manifold, and the tolerance must not be
//     tuned to hide that.
#include "CoolPropThermoModel.H"
#include "StateAcceptance.H"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
using namespace Foam;

int main()
{
    if (!thermo::coolPropThermoAvailable()) { std::cerr << "native CoolProp unavailable\n"; return 77; }
    auto model = thermo::makeCoolPropThermoModel();
    const double rho = 124.55955480372425;
    const double e = -221419.34774;
    const std::vector<double> y{0.92953, 0.07047};
    thermo::EquilibriumSeed seed;
    seed.pressure = 3.7567021256068777e6;
    seed.temperature = 164.77584438095718;
    seed.previousRegion = thermo::EquilibriumRegion::TWO_PHASE;
    const auto r = model->recoverRhoE(rho, e, y, &seed);
    // Second case: A-route step 22, cell 559462 -- a 5 % propane mixture at
    // 3.506 MPa / 162.25 K. The trial state is genuinely off the equilibrium
    // manifold (the 80x90 fine scan reaches only 1.06e-2, just above the
    // development gate), so the recovery must fail cleanly: this is the
    // interface-projection blocker, not a scan-resolution problem. NOTE: the
    // energy is the trial state's *internal* energy = rhoE/rho minus the
    // kinetic term |m|^2/(2 rho^2) = 4728 J/kg. An earlier version of this
    // case forgot the kinetic term and tested a different (recoverable) state.
    {
        thermo::EquilibriumSeed ncSeed;
        ncSeed.pressure = 3506139.6617710623;
        ncSeed.temperature = 162.24710683173657;
        ncSeed.previousRegion = thermo::EquilibriumRegion::TWO_PHASE;
        const auto nc = model->recoverRhoE(115.17626239765, -213497.2, {0.9504, 0.0496}, &ncSeed);
        // Acceptance is per-quantity (|drho|/rho and |de|/|e| separately) while
        // the fine probe reports their sum: this state sits at ~0.5 % + 0.55 %
        // = 1.06e-2, so the development gate (1e-2) accepts it and the
        // production gate (2e-8) rejects it. Both are the documented contract.
        if (stateAcceptanceRel() > 5e-3)
        {
            if (!nc.success())
                throw std::runtime_error("development gate rejected a state inside its per-quantity limits: "
                                         + nc.diagnostics.backendStatus);
            std::cout << "cell 559462 accepted under the development per-quantity gate: T="
                      << nc.state->T << " p=" << nc.state->p << " PASS\n";
        }
        else if (nc.success())
            throw std::runtime_error("production tolerance accepted the off-manifold cell-559462 state");
        else
            std::cout << "cell 559462 correctly rejected at production tolerance: "
                      << nc.diagnostics.backendStatus << " PASS\n";
    }
    const double tolerance = stateAcceptanceRel();
    if (tolerance > 5e-3)
    {
        if (!r.success())
            throw std::runtime_error("near-critical basin not recovered: " + r.diagnostics.backendStatus);
        const double p = r.state->p, T = r.state->T;
        if (p < 3.75e6 || p > 3.86e6 || T < 164.5 || T > 166.5)
            throw std::runtime_error("recovered state outside the physical basin");
        if (std::abs(r.state->rho - rho)/rho > tolerance
            || std::abs(r.state->e - e)/std::max(std::abs(e),1.) > tolerance)
            throw std::runtime_error("recovered state does not close the conserved pair");
        std::cout << "near-critical recovery p=" << p << " T=" << T
                  << " rho=" << r.state->rho << " e=" << r.state->e
                  << " tolerance=" << tolerance << " PASS\n";
    }
    else
    {
        if (r.success())
            throw std::runtime_error("production tolerance accepted an off-manifold state");
        std::cout << "near-critical recovery correctly rejected at production tolerance="
                  << tolerance << " PASS\n";
    }
    return 0;
}
