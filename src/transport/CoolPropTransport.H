#ifndef RFAM_COOLPROP_TRANSPORT_H
#define RFAM_COOLPROP_TRANSPORT_H
// Molecular transport (dynamic viscosity and thermal conductivity) for the
// viscous/turbulence stage.
//
// Single phase: CoolProp 7.2.0 HEOS transport for the fluid/mixture at (p,T)
// (pure-fluid ECS/Chung-style models; validated against the supply table at the
// injector design point: LOX 130 K / 6 MPa -> mu = 8.577e-5 Pa s vs the
// authoritative 8.58e-5; propane 450 K / 4.5 MPa -> 1.4038e-5 vs 1.40e-5).
//
// Two phase: there is no standard mixture rule for an equilibrium two-phase
// cell, so the rule is defined here explicitly and validated by test:
//   * void fraction from the no-slip (homogeneous) lever rule,
//     alpha = (x/rho_v) / (x/rho_v + (1-x)/rho_l);
//   * viscosity: volume-weighted (the single velocity field shares the shear),
//     mu = alpha*mu_v + (1-alpha)*mu_l;
//   * conductivity: harmonic (phase boundaries act in series for conduction),
//     1/lambda = alpha/lambda_v + (1-alpha)/lambda_l.
// Limits are exact: alpha=0 gives the saturated liquid, alpha=1 the saturated
// vapour. Saturated phase properties come from the PQ flash at the cell
// pressure, so the phase compositions follow the mixture equilibrium.
//
// Liquid-branch viscosity: CoolProp's HEOS mixture viscosity returns NaN for
// the O2-rich saturated liquid at every injector pressure (measured 2-5 MPa;
// the pure-fluid and vapour branches are fine). The saturated liquid transport
// therefore falls back to the mass-weighted pure-component viscosities at the
// same (p, T_sat). That is exact in the pure limit and small in the O2-rich
// liquid (the liquid phase carries <=4 % propane by mass), and the fallback is
// reported by the test rather than silently blended.
#include "AbstractState.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace Foam
{
namespace transport
{
struct TransportSample
{
    double mu = 0.;        // dynamic viscosity [Pa s]
    double lambda = 0.;    // thermal conductivity [W/(m K)]
};

// A state already updated to a single-phase PT state.
inline TransportSample fromState(CoolProp::AbstractState& s)
{
    TransportSample out;
    out.mu = s.viscosity();
    out.lambda = s.conductivity();
    if (!(out.mu > 0.) || !std::isfinite(out.mu)
        || !(out.lambda > 0.) || !std::isfinite(out.lambda))
        throw std::runtime_error("COOLPROP_TRANSPORT_INVALID");
    return out;
}

// Per-thread pure-component state used by the liquid-viscosity fallback.
class PureComponentTransport
{
    CoolProp::AbstractState* state_;
public:
    explicit PureComponentTransport(const std::string& fluid)
      : state_(CoolProp::AbstractState::factory("HEOS", fluid)) {}
    ~PureComponentTransport() { delete state_; }
    PureComponentTransport(const PureComponentTransport&) = delete;
    PureComponentTransport& operator=(const PureComponentTransport&) = delete;
    double viscosity(double p, double T)
    {
        state_->update(CoolProp::PT_INPUTS, p, T);
        const double mu = state_->viscosity();
        if (!(mu > 0.) || !std::isfinite(mu))
            throw std::runtime_error("COOLPROP_TRANSPORT_PURE_FALLBACK_FAILED");
        return mu;
    }
    // Saturation pressure at T (0 when the fluid is supercritical there).
    double saturationPressure(double T)
    {
        try
        {
            state_->update(CoolProp::QT_INPUTS, 0., T);
            return state_->p();
        }
        catch (const std::exception&) { return 0.; }
    }
};

// Saturated liquid transport at pressure p on `state` (composition pre-set).
// `state` is left at the Q=0 state.
inline TransportSample saturatedLiquid(CoolProp::AbstractState& state, double p)
{
    state.update(CoolProp::PQ_INPUTS, p, 0.);
    TransportSample out;
    out.lambda = state.conductivity();
    double mu = 0.;
    bool mixtureOk = true;
    try { mu = state.viscosity(); }
    catch (const std::exception&) { mixtureOk = false; }
    if (!mixtureOk || !(mu > 0.) || !std::isfinite(mu))
    {
        // Documented fallback: mass-weighted pure components evaluated in their
        // own stable liquid state. The liquid phase of the O2-rich mixture
        // exists at (p, T_bubble) by Raoult's law; a pure component at the same
        // pressure may be superheated there (e.g. pure O2 at 4.11 MPa is vapour
        // at the mixture's 151.9 K bubble point, whose own saturation pressure
        // is 4.54 MPa). Each component is therefore evaluated at
        // max(p, p_sat,i(T)), which keeps it liquid and reduces exactly to the
        // pure-fluid value when the component is the only one present.
        state.update(CoolProp::PQ_INPUTS, p, 0.);
        const std::vector<double> y = state.get_mass_fractions();
        const double T = state.T();
        static thread_local PureComponentTransport oxygen("Oxygen");
        static thread_local PureComponentTransport propane("Propane");
        const double pO2 = oxygen.saturationPressure(T);
        const double pC3 = propane.saturationPressure(T);
        // 2 % above the saturation pressure: CoolProp refuses a flash exactly
        // on the saturation line, and the property change at constant T over
        // 2 % of pressure is negligible.
        const double muO2 = oxygen.viscosity(std::max(p, 1.02*pO2), T);
        const double muC3 = propane.viscosity(std::max(p, 1.02*pC3), T);
        mu = 0.;
        if (y.size() == 2) mu = y[0]*muO2 + y[1]*muC3;
        else mu = muO2;
        if (std::getenv("RFAM_TRANSPORT_DEBUG"))
            std::fprintf(stderr, "LIQFALLBACK p=%.8g T=%.8g yO2=%.8g yC3=%.8g "
                         "pSatO2=%.8g pSatC3=%.8g muO2=%.8g muC3=%.8g -> mu=%.8g\n",
                         p, T, y.empty()?0.:y[0], y.size()>1?y[1]:0., pO2, pC3, muO2, muC3, mu);
        if (!(mu > 0.) || !std::isfinite(mu))
            throw std::runtime_error("COOLPROP_TRANSPORT_LIQUID_FALLBACK_FAILED");
    }
    out.mu = mu;
    if (!(out.lambda > 0.) || !std::isfinite(out.lambda))
        throw std::runtime_error("COOLPROP_TRANSPORT_INVALID");
    return out;
}

// Saturated vapour transport at pressure p on `state` (composition pre-set).
inline TransportSample saturatedVapour(CoolProp::AbstractState& state, double p)
{
    state.update(CoolProp::PQ_INPUTS, p, 1.);
    return fromState(state);
}

// Two-phase equilibrium blend at the cell pressure. `state` is reused for the
// flashes (the caller owns it and sets the composition).
inline TransportSample twoPhaseBlend(CoolProp::AbstractState& state, double p, double x)
{
    x = std::clamp(x, 0., 1.);
    if (x <= 0.) return saturatedLiquid(state, p);
    if (x >= 1.) return saturatedVapour(state, p);
    const auto liquid = saturatedLiquid(state, p);
    const double rhoL = state.rhomass();
    const auto vapour = saturatedVapour(state, p);
    const double rhoV = state.rhomass();
    const double vL = (1. - x)/rhoL, vV = x/rhoV;
    const double v = vL + vV;
    const double alpha = v > 0. ? vV/v : 0.;
    TransportSample out;
    out.mu = alpha*vapour.mu + (1. - alpha)*liquid.mu;
    out.lambda = 1./(alpha/vapour.lambda + (1. - alpha)/liquid.lambda);
    return out;
}
}
}
#endif
