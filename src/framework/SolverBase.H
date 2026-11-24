#ifndef SolverBase_H
#define SolverBase_H

/*---------------------------------------------------------------------------*\
    SolverBase - framework time-marching shell (Stage 0).

    Holds a Domain reference + an operator chain + a TimeIntegrator.
    Pure orchestration; physics lives in the operators/modules.

        computeTimestep(CFL) = min over operators' maxStableDt
        advanceToTime(tEnd, CFL): march until tEnd / maxIter / trip

    MVP uses this as a concrete composition class (no derivation needed);
    derived solvers come in Stage 1.
\*---------------------------------------------------------------------------*/

#include "Domain.H"
#include "Operator.H"
#include "TimeIntegrator.H"
#include <vector>
#include <memory>
#include <functional>
#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace Foam
{

class SolverBase
{
private:
    Domain& d_;
    std::vector<std::unique_ptr<Operator>> ops_;
    std::unique_ptr<TimeIntegrator> integ_;
    double t_;
    int    iter_;
    int    maxIter_;
    bool previousAcceptedDtPresent_;
    double previousAcceptedDt_;
    std::function<void(const SolverBase&)> monitor_;

public:
    explicit SolverBase(Domain& d)
        : d_(d), integ_(new OperatorSplittingIntegrator()),
          t_(0.0), iter_(0), maxIter_(50000),
          previousAcceptedDtPresent_(false), previousAcceptedDt_(0.0) {}

    void addOperator(std::unique_ptr<Operator> op) { ops_.push_back(std::move(op)); }
    void setIntegrator(std::unique_ptr<TimeIntegrator> ti) { integ_ = std::move(ti); }
    void setMaxIter(int m) { maxIter_ = m; }
    void setMonitor(std::function<void(const SolverBase&)> m) { monitor_ = std::move(m); }

    Domain& domain()    const { return d_; }
    double  time()      const { return t_; }
    int     iteration() const { return iter_; }
    bool hasPreviousAcceptedDt() const
    {
        return previousAcceptedDtPresent_;
    }

    std::optional<double> previousAcceptedDt() const
    {
        if (!previousAcceptedDtPresent_) return std::nullopt;
        return previousAcceptedDt_;
    }

    void setRestartState(
        double time,
        int iteration,
        std::optional<double> previousAcceptedDt = std::nullopt)
    {
        if(!std::isfinite(time)||time<0.0||std::signbit(time)||iteration<0)
            throw std::invalid_argument("invalid solver restart clock state");
        if (previousAcceptedDt.has_value()
            && (!std::isfinite(*previousAcceptedDt)
                || *previousAcceptedDt <= 0.0))
            throw std::invalid_argument(
                "invalid solver previous accepted time step");

        // Validate all fields before publishing any part of the restart
        // clock, so a malformed manifest cannot leave a mixed clock state.
        t_=time;
        iter_=iteration;
        previousAcceptedDtPresent_=previousAcceptedDt.has_value();
        previousAcceptedDt_=previousAcceptedDtPresent_
            ? *previousAcceptedDt : 0.0;
    }

    double computeTimestep(double CFL) const
    {
        double dt = 1e30;
        for (const auto& op : ops_) dt = std::min(dt, op->maxStableDt(d_, CFL));
        return dt;
    }

    // Single explicit step: beginStep + apply all operators + commit.
    // Returns false if the domain tripped.
    bool step(double dt)
    {
        if (!std::isfinite(dt) || dt <= 0.0)
            throw std::invalid_argument("invalid solver time step");
        if (!integ_->step(d_, ops_, dt)) return false;
        const double acceptedDt=integ_->acceptedDt(dt);
        if (!std::isfinite(acceptedDt) || acceptedDt <= 0.0
            || acceptedDt > dt)
            throw std::runtime_error(
                "time integrator returned invalid accepted step");
        if(integ_->advancesPhysicalTime())
        {
            t_ += acceptedDt;
            previousAcceptedDtPresent_=true;
            previousAcceptedDt_=acceptedDt;
        }
        iter_++;
        return true;
    }

    void advanceToTime(double t_end, double CFL)
    {
        while (t_ < t_end && iter_ < maxIter_ && !d_.tripped())
        {
            double dt = computeTimestep(CFL);
            if (!step(dt)) break;
            if (monitor_ && iter_ % 1000 == 0) monitor_(*this);
        }
    }
};

} // End namespace Foam

#endif
