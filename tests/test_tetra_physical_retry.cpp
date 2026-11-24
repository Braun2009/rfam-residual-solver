#include "TetraPhysicalStepRetry.H"
#include <iostream>
#include <stdexcept>

using namespace Foam;

int main()
{
    const std::string boundaryFailure=
        "inlet flow/area/backpressure violates supply entropy for fuel20";
    const auto trial=classifyPhysicalStepFailure(boundaryFailure,true);
    if(trial.action!=PhysicalStepFailureAction::RETRY_WITH_SMALLER_DT
       ||trial.boundaryInputRejected())
        throw std::runtime_error("trial-state boundary excursion was not retryable");

    const auto accepted=classifyPhysicalStepFailure(boundaryFailure,false);
    if(accepted.action!=PhysicalStepFailureAction::REJECT_BOUNDARY_INPUT
       ||!accepted.boundaryInputRejected())
        throw std::runtime_error("accepted-state boundary incompatibility was not terminal");

    const auto numerical=classifyPhysicalStepFailure(
        "physical-time inner LU-SGS did not converge",false);
    if(numerical.action!=PhysicalStepFailureAction::RETRY_WITH_SMALLER_DT)
        throw std::runtime_error("non-boundary numerical failure was not retryable");

    PhysicalStepRetryAudit audit;
    audit.beginStep();
    audit.recordRetry(true);
    audit.recordRetry(true);
    if(audit.stepRetries!=2||audit.boundaryTrialRetries!=2
       ||audit.totalRetries!=2||audit.totalBoundaryTrialRetries!=2)
        throw std::runtime_error("retry audit lost trial-boundary retry counts");
    audit.beginStep();
    if(audit.stepRetries!=0||audit.boundaryTrialRetries!=0
       ||audit.totalRetries!=2||audit.totalBoundaryTrialRetries!=2)
        throw std::runtime_error("new step did not reset only per-step retry counts");

    std::cout<<"TETRA_PHYSICAL_RETRY_PASS\n";
    return 0;
}
