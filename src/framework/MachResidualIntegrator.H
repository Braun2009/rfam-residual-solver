// Copyright 2026 Braun
// SPDX-License-Identifier: Apache-2.0
#ifndef MachResidualIntegrator_H
#define MachResidualIntegrator_H

#include "TimeIntegrator.H"

#include <limits>

namespace Foam
{

// Single-phase, frozen-composition pressure-mode preconditioner. P has
// eigenvalue 1/beta^2 in the pressure mode, 1 in the complementary modes.
// Use (1 + sigma)*(I + sigma*P)^-1 so the pressure correction
// gets alpha = (1 + sigma)/(1 + sigma/beta^2), and other modes are unchanged.
// sigma = physical dt / reference acoustic pseudo dt. This includes the
// physical-time mass term: alpha -> 1 for small dt and -> beta^2 for large dt.
// In conserved coordinates the correction is rank one along the isentrope.
class MachPressurePreconditioner
{
    double floor_;

public:
    explicit MachPressurePreconditioner(double machFloor) : floor_(machFloor)
    {
        if (!(floor_ > 0.0 && floor_ <= 1.0))
            throw std::invalid_argument("Mach floor must be in (0, 1]");
    }

    void apply(const PrimitiveVars& w, const thermo::EquilibriumState& q,
               ConservedVars& correction, double sigma) const
    {
        if (floor_ == 1.0) return;
        if (!(sigma >= 0.0)) throw std::invalid_argument("invalid pseudo-time weight");
        if (sigma == 0.0) return;
        if (!solverConsumerAdmitsAcousticState(q))
            throw std::runtime_error("Mach preconditioning requires admitted acoustics");
        const double velocity2 = w.U[0]*w.U[0] + w.U[1]*w.U[1] + w.U[2]*w.U[2];
        const double mach = std::sqrt(velocity2)/q.acoustic.soundSpeed();
        if (mach >= 1.0) return;
        if (!q.thermodynamics || !q.thermodynamics->valid()
            || q.region == thermo::EquilibriumRegion::TWO_PHASE
            || (q.liquid && q.vapor) || q.condensed)
            throw std::runtime_error("Mach preconditioning requires single-phase EOS derivatives");

        const auto& t = *q.thermodynamics;
        const double dpde = t.dpdT_rho/t.dedT_rho;
        const double dpdrho_e = t.dpdrho_T - dpde*t.dedrho_T;
        const double a2 = dpdrho_e + dpde*(q.p/q.rho)/q.rho;
        if (!(a2 > 0.0 && std::isfinite(a2)))
            throw std::runtime_error("inconsistent isentropic EOS derivatives");

        double drhoe = correction.rhoE + (0.5*velocity2 - q.e)*correction.rho;
        for (std::size_t m = 0; m < 3; ++m) drhoe -= w.U[m]*correction.rhoU[m];
        const double dp = dpdrho_e*correction.rho + (dpde/q.rho)*drhoe;
        const double beta = std::max(mach, floor_);
        const double beta2 = beta*beta;
        // This form also has a well-defined steady limit at sigma = infinity.
        const double pseudoFraction = 1.0/(1.0 + 1.0/sigma);
        const double pressureScale = beta2/(beta2 + (1.0-beta2)*pseudoFraction);
        const double densityChange = (pressureScale - 1.0)*(dp/a2);
        correction.rho += densityChange;
        correction.rhoE += densityChange*(q.h + 0.5*velocity2);
        for (std::size_t m = 0; m < 3; ++m)
            correction.rhoU[m] += densityChange*w.U[m];
        for (std::size_t s = 0; s < w.Y.size(); ++s)
            correction.rhoY[s] += densityChange*w.Y[s];
        correction.rhoZ += densityChange*w.Z;
        correction.rhoC += densityChange*w.C;
    }
};

// Backward Euler: find U such that (U - U_start) + dt*R(U) == 0.
// Only inner corrections are preconditioned; the acceptance test is on the
// original, unpreconditioned equation. No change of physical time or flux.
class MachResidualIntegrator : public TimeIntegrator
{
    MachPressurePreconditioner preconditioner_;
    int maxIterations_;
    double tolerance_;
    double relaxation_;
    mutable int iterations_ = 0;
    mutable double defect_ = std::numeric_limits<double>::infinity();

public:
    explicit MachResidualIntegrator(double machFloor = 0.2,
        int maxIterations = 3000, double tolerance = 1e-12, double relaxation = 0.65)
        : preconditioner_(machFloor), maxIterations_(maxIterations),
          tolerance_(tolerance), relaxation_(relaxation)
    {
        if (maxIterations_ < 1 || !std::isfinite(tolerance_) || tolerance_ <= 0.0
            || !(relaxation_ > 0.0 && relaxation_ <= 1.0))
            throw std::invalid_argument("invalid inner iteration controls");
    }

    int iterations() const { return iterations_; }
    double defect() const { return defect_; }
    std::string name() const override { return "Mach-preconditioned residual backward Euler"; }

    bool step(Domain& d, const std::vector<std::unique_ptr<Operator>>& ops,
              double dt) const override
    {
        iterations_ = 0;
        defect_ = std::numeric_limits<double>::infinity();
        if (!detail::validResidualStep(d, dt, name())) return false;
        detail::StepTransaction transaction(d, ops);
        if (!transaction.active()) return false;
        std::vector<ConservedVars> old;
        old.reserve(d.N());
        for (int cell = 0; cell < d.N(); ++cell) old.push_back(d.U(cell));

        ResidualField work;
        for (;;)
        {
            if (!detail::assembleSpatialResidual(d, ops, work)) return false;
            defect_ = 0.0;
            for (int cell = 0; cell < d.N(); ++cell)
                formDefect(work[cell], d.U(cell), old[cell], dt);
            if (defect_ <= tolerance_)
            {
                transaction.accept();
                return true;
            }
            if (iterations_ == maxIterations_ || !std::isfinite(defect_))
            {
                d.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                    name(),"physical-step residual did not converge");
                return false;
            }

            d.beginStep();
            for (int cell = 0; cell < d.N(); ++cell)
            {
                const auto& w = d.W(cell);
                const auto& q = d.Q(cell);
                try { preconditioner_.apply(w, q, work[cell], acousticWeight(d,w,q,dt)); }
                catch (const std::exception& error)
                {
                    d.recordFailure(Domain::FailureCode::OPERATOR_EXCEPTION,cell,
                        name(),error.what());
                    return false;
                }
                subtractScaledRate(d.Unew(cell), work[cell], relaxation_);
            }
            if (!d.commit()) return false;
            ++iterations_;
        }
    }

private:
    // Directional acoustic rate on the same Cartesian grid as the residual.
    static double acousticWeight(const Domain& d, const PrimitiveVars& w,
                                 const thermo::EquilibriumState& q, double dt)
    {
        const double a = q.acoustic.soundSpeed();
        double rate = (std::abs(w.U[0]) + a)/d.dx();
        if (d.Ny() > 1) rate += (std::abs(w.U[1]) + a)/d.dy();
        if (d.Nz() > 1) rate += (std::abs(w.U[2]) + a)/d.dz();
        return dt*rate;
    }

    void formDefect(ConservedVars& r, const ConservedVars& u,
                    const ConservedVars& old, double dt) const
    {
        const auto component = [&](double& value, double current, double start)
        {
            value = (current - start) + dt*value;
            const double error = std::abs(value)/std::max({1.0, std::abs(current), std::abs(start)});
            defect_ = std::isfinite(error) ? std::max(defect_, error)
                                          : std::numeric_limits<double>::infinity();
        };
        component(r.rho, u.rho, old.rho);
        component(r.rhoE, u.rhoE, old.rhoE);
        component(r.rhoZ, u.rhoZ, old.rhoZ);
        component(r.rhoC, u.rhoC, old.rhoC);
        for (std::size_t m = 0; m < r.rhoU.size(); ++m)
            component(r.rhoU[m], u.rhoU[m], old.rhoU[m]);
        for (std::size_t s = 0; s < r.rhoY.size(); ++s)
            component(r.rhoY[s], u.rhoY[s], old.rhoY[s]);
    }
};

} // namespace Foam
#endif
