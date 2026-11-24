#ifndef SteadyResidualEvaluator_H
#define SteadyResidualEvaluator_H

/*---------------------------------------------------------------------------*\
    SteadyResidualEvaluator - transactional fixed-point residual assembly.

    The evaluator forms the conservative steady residual

        F(U) = [U_candidate(U) - U] / dtReference

    with the production operator ordering. An optional boundary hook executes
    inside the same transaction before residual assembly, so inlet, wall,
    periodic, and steady-outlet ghost states are recomputed from the probed
    state rather than reused from an earlier evaluation.

    A successful evaluation is observational only: Domain U/Unew/W/Q and every
    operator-owned persistent state are rolled back exactly, and the physical
    accepted/rejected step ledger is not modified. A failed evaluation also
    rolls back state atomically while retaining the first Domain failure record.
\*---------------------------------------------------------------------------*/

#include "AcousticStateAdmission.H"
#include "TimeIntegrator.H"
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

namespace Foam
{

class SteadyResidualEvaluator
{
public:
    using BoundaryHook=std::function<bool(Domain&)>;

    struct Result
    {
        std::vector<ConservedVars> field;
        std::vector<unsigned char> active;
        double referenceStep=0.0;
    };

private:
    static void assignScaledDifference(
        ConservedVars& result,
        const ConservedVars& candidate,
        const ConservedVars& current,
        double inverseStep)
    {
        result.rho=(candidate.rho-current.rho)*inverseStep;
        result.rhoE=(candidate.rhoE-current.rhoE)*inverseStep;
        for(std::size_t component=0;component<result.rhoU.size();++component)
            result.rhoU[component]=
                (candidate.rhoU[component]-current.rhoU[component])
                *inverseStep;
        for(std::size_t component=0;component<result.rhoY.size();++component)
            result.rhoY[component]=
                (candidate.rhoY[component]-current.rhoY[component])
                *inverseStep;
        result.rhoZ=(candidate.rhoZ-current.rhoZ)*inverseStep;
        result.rhoC=(candidate.rhoC-current.rhoC)*inverseStep;
    }

public:
    static bool evaluate(
        Domain& domain,
        const std::vector<std::unique_ptr<Operator>>& operators,
        double referenceStep,
        Result& result,
        const BoundaryHook& boundaryHook={},
        bool excludeStructuredHalo=false)
    {
        result=Result{};
        if(!std::isfinite(referenceStep)||referenceStep<=0.0)
            throw std::invalid_argument(
                "steady residual reference step must be positive");
        if(!domain.stateReady()||domain.tripped()
            ||domain.stepTransactionActive())
        {
            if(!domain.tripped())
                domain.recordFailure(
                    Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                    "SteadyResidualEvaluator",
                    "residual evaluation requires a ready Domain outside "
                    "a step transaction");
            return false;
        }
        for(int cell=0;cell<domain.N();++cell)
            if(!solverConsumerAdmitsAcousticState(domain.Q(cell)))
            {
                domain.recordFailure(
                    Domain::FailureCode::STATE_RECOVERY_FAILED,cell,
                    "SteadyResidualEvaluator",
                    "solver acoustic state is not admitted");
                return false;
            }

        const int total=domain.N();
        const int nx=domain.Nx();
        const int ny=domain.Ny();
        const int nz=domain.Nz();
        const int plane=nx*ny;
        result.field.assign(
            static_cast<std::size_t>(total),
            ConservedVars(domain.nSpecies()));
        result.active.assign(static_cast<std::size_t>(total),1U);
        result.referenceStep=referenceStep;
        if(excludeStructuredHalo)
            for(int cell=0;cell<total;++cell)
            {
                const int i=cell%nx;
                const int j=(cell/nx)%ny;
                const int k=cell/plane;
                const bool active=i>0&&i<nx-1
                    &&(ny==1||(j>0&&j<ny-1))
                    &&(nz==1||(k>0&&k<nz-1));
                result.active[static_cast<std::size_t>(cell)]=active?1U:0U;
            }

        bool success=false;
        try
        {
            // A residual probe is not a physical or pseudo-time step. The
            // transaction restores Domain and operator state without adding a
            // rejected-step ledger entry.
            detail::StepTransaction transaction(domain,operators,false);
            if(!transaction.active()) return false;
            if(boundaryHook)
            {
                bool boundaryReady=false;
                try
                {
                    boundaryReady=boundaryHook(domain);
                }
                catch(const std::exception& error)
                {
                    if(!domain.tripped())
                        domain.recordFailure(
                            Domain::FailureCode::OPERATOR_EXCEPTION,-1,
                            "SteadyResidualEvaluator boundary",error.what());
                    return false;
                }
                catch(...)
                {
                    if(!domain.tripped())
                        domain.recordFailure(
                            Domain::FailureCode::OPERATOR_EXCEPTION,-1,
                            "SteadyResidualEvaluator boundary",
                            "unknown boundary exception");
                    return false;
                }
                if(!boundaryReady||domain.tripped())
                {
                    if(!domain.tripped())
                        domain.recordFailure(
                            Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                            "SteadyResidualEvaluator boundary",
                            "transactional boundary hook rejected the state");
                    return false;
                }
            }

            domain.beginStep();
            if(domain.tripped()||!transaction.apply(referenceStep))
                return false;
            const double inverseStep=1.0/referenceStep;
            for(int cell=0;cell<total;++cell)
            {
                if(result.active[static_cast<std::size_t>(cell)]==0U)
                    continue;
                assignScaledDifference(
                    result.field[static_cast<std::size_t>(cell)],
                    domain.Unew(cell),domain.U(cell),inverseStep);
            }
            success=true;
            // transaction destructor performs rollback(false).
        }
        catch(...)
        {
            // StepTransaction/Operator already publish deterministic failure
            // metadata where possible. Preserve it while rollback restores the
            // authoritative state.
            return false;
        }
        return success&&!domain.tripped();
    }
};

} // End namespace Foam

#endif
