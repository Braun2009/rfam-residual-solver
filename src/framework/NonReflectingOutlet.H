#ifndef NonReflectingOutlet_H
#define NonReflectingOutlet_H

/*---------------------------------------------------------------------------*\
    NonReflectingOutlet - LODI-based subsonic outflow boundary (2D, +x face).

    A plain zero-gradient outlet (W[Nx-1]=W[Nx-2]) reflects outgoing acoustic
    waves back into the domain, which corrupts transient combustion / plume /
    thermoacoustic behaviour. This applies the Local One-Dimensional Inviscid
    (LODI) relations of Poinsot & Lele (1992) at the outlet face.

    For a subsonic OUTFLOW (0 < u < a) the x-normal characteristic waves are:
        L1 : (u - a)  incoming  (the ONE wave entering from outside)
        L2 : (u)      outgoing  (entropy)
        L3 : (u)      outgoing  (transverse velocity v)
        L5 : (u + a)  outgoing  (downstream acoustic)
    The three outgoing waves are computed from interior one-sided differences;
    the single incoming wave L1 is NOT extrapolated (that is what reflects) but
    set by a partial-reflection relation that relaxes pressure to p_inf:
        L1 = K (p - p_inf),   K = sigma (1 - M^2) a / L
    (Rudy & Strikwerda relaxation; sigma ~ 0.25). K=0 is perfectly
    non-reflecting; a small K prevents slow drift of the mean outlet pressure.

    Header-only, pure std C++. Operates directly on the Domain halo column via
    setCell, exactly like the hand-written BCs in solver2d_frame.
\*---------------------------------------------------------------------------*/

#include "Domain.H"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <utility>

namespace Foam
{

class NonReflectingOutlet
{
private:
    double p_inf_;     // far-field (target) static pressure [Pa]
    double sigma_;     // relaxation coefficient (0 = perfectly non-reflecting)
    double L_;         // characteristic domain length [m] (for the K scaling)
    std::vector<std::string> species_;   // EOS species (empty = legacy path)

    static bool supportsCharacteristicState(
        const thermo::EquilibriumState& state)
    {
        return state.valid()
            && state.acoustic.closure
                == thermo::AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE;
    }

public:
    NonReflectingOutlet(double p_inf, double L,
                        const std::vector<std::string>& species = {},
                        double sigma = 0.25)
        : p_inf_(p_inf), sigma_(sigma), L_(L), species_(species) {}

    // Apply to the +x outlet column (i = Nx-1), using interior columns
    // i=Nx-2 (adjacent) and i=Nx-3 (for the one-sided x-derivative).
    // dt is the current timestep. Rows j span the full column.
    void apply(Domain& d, double dt) const
    {
        const int Nx = d.Nx(), Ny = d.Ny();
        if (Nx < 3) return;                     // need 2 interior columns
        if (!d.stateReady() || d.tripped() || d.stepTransactionActive())
            throw std::logic_error(
                "NonReflectingOutlet requires a ready Domain outside a step transaction");
        if (!species_.empty() && species_ != d.speciesNames())
            throw std::invalid_argument(
                "NonReflectingOutlet species layout does not match the Domain");
        if (!d.thermoModel().supportsPrimitivePTYInitialization())
            throw std::runtime_error(
                "NonReflectingOutlet phase-aware characteristic closure is not implemented");
        const double dx = d.dx();

        struct BoundaryCandidate
        {
            int row = 0;
            PrimitiveVars primitive;
            thermo::EquilibriumState equilibrium;
        };
        std::vector<BoundaryCandidate> candidates;
        candidates.reserve(static_cast<std::size_t>(Ny));

        auto initializeState = [&](const PrimitiveVars& primitive)
        {
            thermo::EquilibriumResult result =
                d.thermoModel().initializePTY(
                    primitive.p,primitive.T,primitive.Y);
            if (!result.success() || !supportsCharacteristicState(*result.state))
                throw std::runtime_error(
                    "NonReflectingOutlet PTY recovery failed or returned an "
                    "unsupported acoustic closure");
            return std::move(*result.state);
        };

        for (int j = 0; j < Ny; ++j)
        {
            // NSCBC ghost-node LODI: evolve the boundary node's OWN persistent
            // state W(Nx-1) using the one-sided gradient AT the boundary
            // (W(Nx-1) - W(Nx-2))/dx. This lets the outgoing wave leave cleanly;
            // reconstructing the boundary FROM the interior each step instead
            // just re-injects the outgoing wave (= reflection).
            const PrimitiveVars& wb = d.W(Nx-1, j);   // boundary node (persistent)
            const PrimitiveVars& wi = d.W(Nx-2, j);   // first interior cell
            const thermo::EquilibriumState& qb = d.Q(Nx-1,j);
            const thermo::EquilibriumState& qi = d.Q(Nx-2,j);
            if (!supportsCharacteristicState(qb)
                || !supportsCharacteristicState(qi))
            {
                // Two-phase boundary state (HEM projection, e.g. the Mayer
                // transcritical jet): NSCBC LODI characteristic analysis is
                // frozen-composition single-phase only, so degrade to the
                // zero-gradient (extrapolated) outflow path with an
                // explicit audit print (2026-08-04, Gate 5).
#ifdef RFAM_THERMO_DIAGNOSTICS
                std::cerr << "[OUTLET-DIAG] two-phase boundary state -> "
                             "zero-gradient outflow fallback\n";
#endif
                candidates.push_back({j, wi, qi});
                continue;
            }

            double rho = qb.rho;
            double a   = qb.acoustic.soundSpeed();
            double u   = wb.U[0];

            // Supersonic outflow: everything leaves, pure extrapolation.
            if (u >= a)
            {
                candidates.push_back({j,wi,qi});
                continue;
            }
            // Reverse flow at an outlet is out of scope: zero-gradient fallback.
            if (u <= 0.0)
            {
                candidates.push_back({j,wi,qi});
                continue;
            }

            // One-sided x-derivatives at the boundary (backward toward interior).
            double dpdx   = (wb.p    - wi.p)    / dx;
            double dudx   = (wb.U[0] - wi.U[0]) / dx;
            double drhodx = (rho - qi.rho) / dx;

            // Outgoing waves from interior data: L5 (u+a acoustic), L2 (entropy).
            double L5 = (u + a) * (dpdx + rho * a * dudx);
            double L2 = u * (a*a*drhodx - dpdx);

            // Incoming acoustic wave L1 (u-a): NOT extrapolated. LODI relaxation
            // to p_inf (Rudy-Strikwerda); K=0 is perfectly non-reflecting.
            double M = u / a;
            double K = sigma_ * (1.0 - M*M) * a / L_;
            double L1 = K * (wb.p - p_inf_);

            // LODI time-advance of the boundary node's own primitives:
            //   dp/dt   = -1/2 (L5 + L1)
            //   du/dt   = -1/(2 rho a) (L5 - L1)
            //   drho/dt = -1/a^2 [ L2 + 1/2 (L5 + L1) ]
            double dpdt   = -0.5 * (L5 + L1);
            double dudt   = -(L5 - L1) / (2.0 * rho * a);
            double drhodt = -(L2 + 0.5 * (L5 + L1)) / (a*a);

            double p_new   = wb.p    + dt * dpdt;
            double u_new   = wb.U[0] + dt * dudt;
            double rho_new = rho     + dt * drhodt;

            if (p_new   < 1e3)  p_new   = 1e3;
            if (rho_new < 1e-6) rho_new = 1e-6;

            PrimitiveVars wout = wb;
            wout.p    = p_new;
            wout.U[0] = u_new;
            // Recover T so the EOS density matches rho_new (bisection on T).
            {
                double Tlo = 50.0, Thi = 6000.0;
                auto rhoErr = [&](double Tq){
                    PrimitiveVars t = wout; t.T = Tq;
                    return initializeState(t).rho - rho_new;
                };
                double flo = rhoErr(Tlo), fhi = rhoErr(Thi);
                if (flo * fhi < 0.0)
                {
                    double T = std::clamp(wb.T,Tlo,Thi);
                    for (int it = 0; it < 60; ++it)
                    {
                        double fm = rhoErr(T);
                        if (fm > 0.0) Tlo = T; else Thi = T;
                        double Tn = 0.5*(Tlo+Thi);
                        if (std::abs(Tn-T) < 1e-4) { T = Tn; break; }
                        T = Tn;
                    }
                    wout.T = T;
                }
                // else: keep wb.T (degrade gracefully)
            }
            candidates.push_back({j,wout,initializeState(wout)});
        }

        for (BoundaryCandidate& candidate : candidates)
        {
            if (!d.setCell(
                    Nx-1,candidate.row,candidate.primitive,
                    candidate.equilibrium))
                throw std::runtime_error(
                    "NonReflectingOutlet failed to publish a complete state: "
                    +d.failure().detail);
        }
    }

    void setPinf(double p) { p_inf_ = p; }
    double pinf() const { return p_inf_; }
};

} // End namespace Foam

#endif
