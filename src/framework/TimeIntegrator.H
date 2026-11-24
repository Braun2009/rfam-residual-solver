// Copyright 2026 Braun
// SPDX-License-Identifier: Apache-2.0
#ifndef TimeIntegrator_H
#define TimeIntegrator_H

/*---------------------------------------------------------------------------*\
    TimeIntegrator - time advance over an operator chain.

    OperatorSplittingIntegrator performs the ADDITIVE update that reproduces
    solver1d_fpv's single synchronous step:
        beginStep()  : Unew <- U
        apply ops    : each accumulates increments from the SAME t^n state
        commit()     : U <- Unew (+ NaN check + refresh primitives)
\*---------------------------------------------------------------------------*/

#include "Domain.H"
#include "Operator.H"
#include <cmath>
#include <exception>
#include <memory>
#include <string>
#include <vector>

namespace Foam
{

class TimeIntegrator
{
public:
    virtual ~TimeIntegrator() = default;
    virtual bool step(
        Domain& d,
        const std::vector<std::unique_ptr<Operator>>& ops,
        double dt
    ) const = 0;
    // The actual physical step accepted by step(). Most integrators accept
    // the requested value; adaptive integrators may accept a smaller retry.
    // SolverBase publishes this value to its clock and restart provenance.
    virtual double acceptedDt(double requestedDt) const
    {
        return requestedDt;
    }
    // Steady/pseudo-time integrators advance nonlinear iterations without
    // publishing pseudo time as physical clock time.
    virtual bool advancesPhysicalTime() const { return true; }
    virtual std::string name() const = 0;
};

namespace detail
{

class StepTransaction
{
    Domain& domain_;
    const std::vector<std::unique_ptr<Operator>>& operators_;
    std::size_t begun_ = 0;
    bool active_ = false;
    bool countRejectedOnRollback_ = true;

    void rollback() noexcept
    {
        if (!active_) return;
        for (std::size_t i = begun_; i > 0; --i)
            operators_[i-1]->rollbackStepTransaction();
        domain_.rollbackStepTransaction(countRejectedOnRollback_);
        active_ = false;
    }

public:
    StepTransaction(
        Domain& domain,
        const std::vector<std::unique_ptr<Operator>>& operators,
        bool countRejectedOnRollback=true)
        : domain_(domain), operators_(operators),
          countRejectedOnRollback_(countRejectedOnRollback)
    {
        if (!domain_.beginStepTransaction()) return;
        active_ = true;
        try
        {
            for (; begun_ < operators_.size(); ++begun_)
                operators_[begun_]->beginStepTransaction(domain_);
        }
        catch (const std::exception& error)
        {
            domain_.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                begun_ < operators_.size() ? operators_[begun_]->name()
                                            : "operator transaction",
                error.what());
            rollback();
            throw;
        }
        catch (...)
        {
            domain_.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                begun_ < operators_.size() ? operators_[begun_]->name()
                                            : "operator transaction",
                "unknown transaction setup failure");
            rollback();
            throw;
        }
    }

    ~StepTransaction() { rollback(); }

    bool active() const { return active_; }

    bool apply(double dt)
    {
        for (const auto& op : operators_)
        {
            try
            {
                op->apply(domain_,dt);
            }
            catch (const std::exception& error)
            {
                domain_.recordFailure(Domain::FailureCode::OPERATOR_EXCEPTION,-1,
                    op->name(),error.what());
                throw;
            }
            catch (...)
            {
                domain_.recordFailure(Domain::FailureCode::OPERATOR_EXCEPTION,-1,
                    op->name(),"unknown operator exception");
                throw;
            }
            if (domain_.tripped()) return false;
        }
        return true;
    }

    void averageWithStepStartState()
    {
        for (std::size_t i = 0; i < begun_; ++i)
        {
            try
            {
                operators_[i]->averageWithStepStartState();
            }
            catch (const std::exception& error)
            {
                domain_.recordFailure(
                    Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                    operators_[i]->name(),error.what());
                throw;
            }
            catch (...)
            {
                domain_.recordFailure(
                    Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                    operators_[i]->name(),
                    "unknown RK2 state averaging failure");
                throw;
            }
        }
    }

    void accept() noexcept
    {
        if (!active_) return;
        for (std::size_t i = 0; i < begun_; ++i)
            operators_[i]->acceptStepTransaction();
        domain_.acceptStepTransaction();
        active_ = false;
    }
};

} // namespace detail

class OperatorSplittingIntegrator : public TimeIntegrator
{
public:
    bool step(
        Domain& d,
        const std::vector<std::unique_ptr<Operator>>& ops,
        double dt
    ) const override
    {
        detail::StepTransaction transaction(d,ops);
        if (!transaction.active()) return false;
        d.beginStep();
        if (!transaction.apply(dt) || !d.commit()) return false;
        transaction.accept();
        return true;
    }
    std::string name() const override { return "OperatorSplitting (additive Euler)"; }
};

namespace detail
{
// All residual operators see the same thermodynamically recovered state.
inline bool assembleSpatialResidual(Domain& d,
    const std::vector<std::unique_ptr<Operator>>& ops, ResidualField& R)
{
    zeroResidual(R, d.N(), d.nSpecies());
    for (const auto& op : ops)
    {
        try
        {
            if (!op->supportsResidual())
                throw std::logic_error("Operator does not support residual assembly");
            op->residual(d, R);
        }
        catch (const std::exception& error)
        {
            d.recordFailure(Domain::FailureCode::OPERATOR_EXCEPTION,-1,
                op->name(),error.what());
            return false;
        }
        catch (...)
        {
            d.recordFailure(Domain::FailureCode::OPERATOR_EXCEPTION,-1,
                op->name(),"unknown residual exception");
            return false;
        }
        if (d.tripped()) return false;
    }
    return true;
}

inline bool validResidualStep(Domain& d, double dt, const std::string& name)
{
    if (std::isfinite(dt) && dt > 0.0) return true;
    d.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                    name,"dt must be finite and positive");
    return false;
}
} // namespace detail

// Explicit Euler. Assemble div(F) - S, then Unew = U - dt*R.
class ResidualIntegrator : public TimeIntegrator
{
public:
    bool step(Domain& d, const std::vector<std::unique_ptr<Operator>>& ops,
              double dt) const override
    {
        if (!detail::validResidualStep(d, dt, name())) return false;
        detail::StepTransaction transaction(d, ops);
        if (!transaction.active()) return false;
        ResidualField residual;
        if (!detail::assembleSpatialResidual(d, ops, residual)) return false;
        d.beginStep();
        for (int cell = 0; cell < d.N(); ++cell)
            subtractScaledRate(d.Unew(cell), residual[cell], dt);
        if (!d.commit()) return false;
        transaction.accept();
        return true;
    }
    std::string name() const override { return "Residual Euler"; }
};

/*---------------------------------------------------------------------------*\
    RK2Integrator - second-order Runge-Kutta (Heun's method).

    Two-stage predictor-corrector:
        U_pred = U^n + dt * f(U^n)            (stage 1)
        U^{n+1} = U^n + dt/2 * [f(U^n) + f(U_pred)]  (stage 2)

    Implementation:
      1. Save U^n
      2. beginStep() -> Unew = U
      3. apply ops with dt -> Unew += dt*f(U^n)  [predictor]
      4. commit() -> U = Unew (predictor becomes current state)
      5. beginStep() -> Unew = U_pred
      6. apply ops with dt -> Unew += dt*f(U_pred)
      7. Average: Unew = 0.5*(U_saved + Unew)
      8. commit() final
\*---------------------------------------------------------------------------*/
class RK2Integrator : public TimeIntegrator
{
public:
    bool step(
        Domain& d,
        const std::vector<std::unique_ptr<Operator>>& ops,
        double dt
    ) const override
    {
        detail::StepTransaction transaction(d,ops);
        if (!transaction.active()) return false;

        // Stage 1: save U^n, compute predictor = U^n + dt*f(U^n)
        d.beginStep();          // Unew <- U
        if (!transaction.apply(dt) || !d.commit()) return false;

        // Stage 2: compute f(U_pred), then average with saved U^n
        d.beginStep();          // Unew <- U_pred
        if (!transaction.apply(dt)) return false;
        // Now Unew = U_pred + dt*f(U_pred)
        // We need: U^{n+1} = U^n + dt/2*[f(U^n) + f(U_pred)]
        //                  = 0.5*[U^n + dt*f(U^n)] + 0.5*[U_pred + dt*f(U_pred)]
        // But we have: Unew_stage1 = U^n + dt*f(U^n) (now in Usaved_)
        //              Unew_stage2 = U_pred + dt*f(U_pred) (current Unew)
        // Wait, U_pred = U^n + dt*f(U^n), so:
        //   Unew_stage2 = (U^n + dt*f(U^n)) + dt*f(U_pred)
        //               = U^n + dt*[f(U^n) + f(U_pred)]
        // We need 0.5 of that plus 0.5*U^n = U^n + 0.5*dt*[f(U^n)+f(U_pred)]
        // Actually simpler: average Usaved (=U^n) with current Unew (=U^n+dt*[f1+f2])
        // gives (U^n + U^n + dt*[f1+f2])/2 = U^n + dt/2*[f1+f2]. Correct!
        d.averageWithStepStartState();
        if (d.tripped()) return false;
        transaction.averageWithStepStartState();
        if (d.tripped() || !d.commit()) return false;
        transaction.accept();
        return true;
    }
    std::string name() const override { return "RK2 (Heun)"; }
};

} // End namespace Foam

#endif
