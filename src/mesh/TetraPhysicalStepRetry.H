#ifndef TetraPhysicalStepRetry_H
#define TetraPhysicalStepRetry_H

#include <string>

namespace Foam
{
enum class PhysicalStepFailureAction
{
    RETRY_WITH_SMALLER_DT,
    REJECT_BOUNDARY_INPUT
};

struct PhysicalStepRetryDecision
{
    PhysicalStepFailureAction action;
    bool boundaryInputRejected() const
    {return action==PhysicalStepFailureAction::REJECT_BOUNDARY_INPUT;}
};

inline bool isNonRetryablePhysicalBoundaryFailure(const std::string& text)
{
    return text.find("inlet flow/area/backpressure violates supply entropy")!=std::string::npos
        ||text.find("inlet backpressure exceeds supply")!=std::string::npos
        ||text.find("inlet total enthalpy failed to close")!=std::string::npos
        ||text.find("static inlet flash failed")!=std::string::npos;
}

inline PhysicalStepRetryDecision classifyPhysicalStepFailure(
    const std::string& text,bool boundaryTrialRejected)
{
    if(isNonRetryablePhysicalBoundaryFailure(text)&&!boundaryTrialRejected)
        return {PhysicalStepFailureAction::REJECT_BOUNDARY_INPUT};
    return {PhysicalStepFailureAction::RETRY_WITH_SMALLER_DT};
}

struct PhysicalStepRetryAudit
{
    int stepRetries=0;
    int boundaryTrialRetries=0;
    int totalRetries=0;
    int totalBoundaryTrialRetries=0;

    void beginStep()
    {stepRetries=0;boundaryTrialRetries=0;}

    void recordRetry(bool boundaryTrial)
    {
        ++stepRetries;
        ++totalRetries;
        if(boundaryTrial)
        {
            ++boundaryTrialRetries;
            ++totalBoundaryTrialRetries;
        }
    }
};
}
#endif
